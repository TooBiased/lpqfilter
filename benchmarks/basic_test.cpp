/*******************************************************************************
 * benchmarks/basic_test.cpp
 *
 * This basic test fills a filter datastructure using p threads. Inputs:
 * 	Flags:  -growing -dynamic control if the table would need to grow
 *          -p       <# of threads>
 *          -n       <# of elements inserted/queried>
 *          -cap     <# initial capacity of the table>
 *          -min_cap <# starting point for growing tables>
 *          -it      <# number of outer repetitions>  // different seeds
 *          -sub_it  <# number of inner iterations> // same seed
 *
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "filter_select.hpp"
#include "benchmark_utils.hpp"

#include "utils/pin_thread.hpp"
#include "utils/command_line_parser.hpp"
#include "utils/thread_coordination.hpp"

#include "implementation/definitions.hpp"
#include "implementation/composite_filter/growing_qfilter.hpp"
#include "implementation/composite_filter/dynamic_qfilter.hpp"

#include <random>
#include <iostream>
#include <tuple>

namespace utm = utils_tm;
namespace ttm = utils_tm::thread_tm;
namespace otm = utils_tm::out_tm;

template<typename T>
inline void out(T t, size_t space)
{
    std::cout.width(space);
    std::cout << t << " " << std::flush;
}

void print_column_names()
{
    out("#i"        , 3);
    out("p"         , 3);
    out("bits"      , 7);
    out("n"         , 9);
    out("capacity"  , 10);
    out("t_insert"  , 10);
    out("t_cont_+"  , 10);
    out("t_cont_-"  , 10);
    out("fp_rate_%" , 10);
    out("error"     , 6);
    std::cout << std::endl;
}

void print_parameters(size_t it, size_t p, size_t bits, size_t n, size_t cap)
{
    out(it    , 3);
    out(p     , 3);
    out(bits  , 7);
    out(n     , 9);
    out(cap   , 10);
}


using key_type = size_t;

alignas(128)  static std::atomic_size_t current_block;
/* no need */ static uint64_t*          keys;
/* no need */ static uint64_t*          keys_succ;
/* no need */ static uint64_t*          keys_unsucc;
/* no need */ static std::atomic_size_t errors;
/* no need */ static std::atomic_size_t false_positives;

std::atomic_size_t  random_seed = 123;


/* STAGES OF THE MAIN TEST ****************************************************/
int generate_keys(size_t n)
{
    std::uniform_int_distribution<uint64_t> dis(1,1ull<<63);

    ttm::execute_blockwise_parallel(current_block, n,
        [&dis](size_t s, size_t e)
        {
            std::mt19937_64 re(s*random_seed);

            for (size_t i = s; i < e; i++)
            {
                keys[i] = dis(re);
            }
        });

    return 0;
}

int fill_succ_keys(size_t n)
{
	std::uniform_int_distribution<uint64_t> dis(0, n - 1);

	ttm::execute_blockwise_parallel(current_block, n,
							   [&dis](size_t s, size_t e)
	{
		std::mt19937_64 re(s*random_seed);

		for (size_t i = s; i < e; i++)
		{
			keys_succ[i] = keys[dis(re)];
		}
	});

	return 0;
}

template<class QuotientFilter>
int insert(QuotientFilter& filter, size_t end)
{
    ttm::execute_parallel(current_block, end,
        [&filter](size_t i)
        {
			const auto key = keys[i];

			if (!filter.insert(key))
			{
				//std::cerr << "i\n";
			}
        });

    return 0;
}

template<class QuotientFilter>
int contains_successful(QuotientFilter& filter, size_t end)
{
    size_t err = 0;
    ttm::execute_parallel(current_block, end,
        [&filter, &err](size_t i)
        {
            const auto key = keys[i];
            if (!filter.contains(key)) {
                err++;
                //std::cerr << "c";
            }
        });

    errors.fetch_add(err, std::memory_order_relaxed);
    return 0;
}

template<class QuotientFilter>
int contains_unsuccessful(QuotientFilter& filter, size_t end)
{
    size_t fp = 0;
    ttm::execute_parallel(current_block, end,
        [&filter, &fp](size_t i)
        {
            const auto key = keys_unsucc[i];
            if (filter.contains(key))
				fp++;
        });

    false_positives.fetch_add(fp, std::memory_order_relaxed);
    return 0;
}




/* MAIN TEST FUNCTION *********************************************************/
template <class FilterType>
struct test
{
    template <class ThreadType>
    struct filter_test
    {
        static int execute(ThreadType t, FilterType& filter, basic_benchmark_config config)
        {
            using HandleType = qf::HandleType<FilterType>;
            utm::pin_to_core(t.id);

            if (ThreadType::is_main)
            {
                keys        = new uint64_t[2 * config.n];
                keys_unsucc = &keys[config.n];
                keys_succ   = new uint64_t[config.n];
            }

            for (size_t i = 0; i < config.it * config.sub_it; ++i)
            {
                if (ThreadType::is_main)
                {
                    ++random_seed;
                    current_block.store(0);
                }

                if (i % config.sub_it == 0)
                {
                    t.synchronized(generate_keys, 2 * config.n);
                    if (ThreadType::is_main)
                        current_block.store(0);
                    t.synchronized(fill_succ_keys, config.n);
                }

                // STAGE 0.3 - Reinitialize the quotient filter + some synchronization
                if (ThreadType::is_main)
                {
                    FilterType newlocal = construct_quotient_filter<FilterType>(config);
                    filter = std::move(newlocal);
                }

                t.synchronize();

                auto handle = qf::create_handle(filter);
                auto start_cap = filter.capacity(); // use real cell count

                if (ThreadType::is_main)
                    print_parameters(i, config.p, config.bits, config.n, start_cap);

                t.synchronize();


                // STAGE 1 - Insert
                {
                    if (ThreadType::is_main)
                        current_block.store(0);

                    auto dur = t.synchronized(insert<HandleType>, handle, config.n);
                    t.out << otm::width(10) << dur.second / 1000000.;
                }

                // STAGE 2 - Successful Contains
                {
                    if (ThreadType::is_main)
                        current_block.store(0);

                    auto dur = t.synchronized(contains_successful<HandleType>, handle, config.n);
                    t.out << otm::width(10) << dur.second / 1000000.;
                }

                // STAGE 3 - Unsuccessful Contains
                {
                    if (ThreadType::is_main)
                        current_block.store(0);

                    auto dur = t.synchronized(contains_unsuccessful<HandleType>, handle, config.n);
                    t.out << otm::width(10) << dur.second / 1000000.
                          << otm::width(10) << double(false_positives.load()) / config.n * 100
                          << otm::width(6)  << errors.load() << std::endl;
                }

                if (ThreadType::is_main)
                {
                    errors.store(0);
                    false_positives.store(0);
                }

                t.synchronize();
            }

            if (ThreadType::is_main)
            {
                delete[] keys;
                delete[] keys_succ;
            }

            return 0;
        }
    };
};



template <class Type, bool instanciated>
class instanciation_handler
{
public:
    using type = Type;
};

template <class Type>
class instanciation_handler<Type,true>
{
public:
    using type = typename Type::template instanciated<qf::DEFAULT_REMAINDER_BITS>;
};

template <template <class> class Wrapper, class Type, bool growing_compatible>
class growing_compatibility_handler
{
public:
    using type = typename instanciation_handler<Type, Type::is_templated>::type;

};

template <template <class> class Wrapper, class Type>
class growing_compatibility_handler<Wrapper, Type, true>
{
public:
    using type = Wrapper<Type>;
};


/* MAIN FUNCTION: READS COMMANDLINE AND STARTS TEST THREADS *******************/
int main(int argn, char** argc)
{
	basic_benchmark_config config;
    utm::command_line_parser c{argn, argc};

	config.growing       = c.bool_arg("-growing");
	config.dynamic       = c.bool_arg("-dynamic");
	config.p             = c.int_arg("-p"         , config.p);
    config.n             = c.int_arg("-n"         , config.n);
	config.cap           = c.int_arg("-cap"       , 2 * config.n);
	config.min_cap       = c.int_arg("-min_cap"   , config.n / 4);
    config.it            = c.int_arg("-it"        , config.it);
	config.sub_it        = c.int_arg("-sub_it"    , config.sub_it);

    double fp_rate = c.double_arg("-fp_rate_%", 0.0);
    size_t rbits = c.int_arg("-r_bits", 0);

    if (rbits != 0 && fp_rate != 0.0)
    {
		std::cout << "Error: only one of \"-fp_rate_%\" and \"-r_bits\" can be used\n";
		return 1;
    }
    else if (rbits != 0)
    {
        config.set_rbits(rbits);
    }
    else if (fp_rate != 0.0)
    {
        config.set_fp_rate_percent(fp_rate);
    }

    if (!c.report()) return 1;

	if (config.growing && config.dynamic)
	{
		std::cout << "Error: only one of \"-growing\" or \"-dynamic\" can be used\n";
		return 1;
	}

	print_column_names();

	if (config.growing)
	{
        if constexpr (! QF_type<key_type>::is_growing_compatible)
        {
            std::cout << "Error: this quotient filter is not growing_compatible\n";
            return 1;
        }
        using qfilter_type = typename growing_compatibility_handler<
            qf::growing_quotient_filter, QF_type<key_type>,
            QF_type<key_type>::is_growing_compatible>::type;

        qfilter_type qfilter{};
        return ttm::start_threads<test<qfilter_type>::template filter_test>
            (config.p, std::ref(qfilter), config);
	}
	else if (config.dynamic)
	{
        if constexpr (! QF_type<key_type>::is_growing_compatible)
        {
            std::cout << "Error: this quotient filter is not growing_compatible\n";
            return 1;
        }
        using qfilter_type = typename growing_compatibility_handler<
            qf::dynamic_quotient_filter, QF_type<key_type>,
            QF_type<key_type>::is_growing_compatible>::type;

        qfilter_type qfilter{};
        return ttm::start_threads<test<qfilter_type>::template filter_test>
            (config.p, std::ref(qfilter), config);
	}
	else
	{
        using qfilter_type = typename instanciation_handler<
            QF_type<key_type>,
            QF_type<key_type>::is_templated>::type;

        qfilter_type qfilter{};
        return ttm::start_threads<test<qfilter_type>::template filter_test>
            (config.p, std::ref(qfilter), config);
	}

    return 0;
}
