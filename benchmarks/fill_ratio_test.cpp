#include "filter_select.hpp"
#include "benchmark_utils.hpp"

#include "utils/pin.h"
#include "utils/commandline.h"
#include "utils/test_coordination.h"

#include <random>
#include <iostream>
#include <tuple>

template<typename T>
inline void out(T t, size_t space)
{
    std::cout.width(space);
    std::cout << t << " " << std::flush;
}

void print_column_names()
{
    out("#i"       , 3);
    out("p"        , 3);
	out("bits"     , 7);
    out("n"        , 9);
    out("test_n"   , 6);
    out("f%"       , 4);
    out("t_cont_+" , 10);
    out("t_cont_-" , 10);
    out("t_insert" , 10);
    out("fp_rate_%", 10);
    out("error"    , 6);
    std::cout << std::endl;
}

void print_parameters(size_t it, size_t p, size_t bits, size_t n, size_t tsize)
{
    out(it    , 3);
    out(p     , 3);
    out(bits  , 7);
    out(n     , 9);
    out(tsize , 6);
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
	std::uniform_int_distribution<uint64_t> dis(1, 1ull << 63);

	execute_blockwise_parallel(current_block, n,
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

    execute_blockwise_parallel(current_block, n,
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
	execute_parallel(current_block, end,
	    [&filter](size_t i)
		{
			auto key = keys[i];

			if (!filter.insert(key))
			{
				std::cerr << "i\n";
			}
		});

	return 0;
}

template<class QuotientFilter>
int contains_successful(QuotientFilter& filter, uint64_t* lbuffer, size_t test_size)
{
    size_t err = 0;
    for (size_t i = 0; i < test_size; ++i)
    {
        auto key = lbuffer[i];
		if (!filter.contains(key))
		{
			err++;
			std::cerr << "c";
		}
    }

    errors.fetch_add(err, std::memory_order_relaxed);
    return 0;
}

template<class QuotientFilter>
int contains_unsuccessful(QuotientFilter& filter, uint64_t* lbuffer, size_t test_size)
{
    size_t fp = 0;
    for (size_t i = 0; i < test_size; ++i)
    {
        auto key = lbuffer[i];
        if (filter.contains(key)) fp++;
    }

    false_positives.fetch_add(fp, std::memory_order_relaxed);
    return 0;
}

template<class QuotientFilter>
int insert_test(QuotientFilter& filter, uint64_t* lbuffer, size_t test_size)
{
    for (size_t i = 0; i < test_size; ++i)
    {
        auto key = lbuffer[i];
        filter.insert(key);
    };

    return 0;
}

/* GENERATE KEYS FOR TESTS ****************************************************/

void fill_loc_succ(uint64_t* lbuffer, size_t id, size_t inserted, size_t test_size)
{
    std::mt19937_64 re((id + 1) * random_seed);
    std::uniform_int_distribution<uint64_t> dis(0,inserted - 1);

    for (size_t i = 0; i < test_size; ++i)
    {
        lbuffer[i] = keys[dis(re)];
    }
}

void fill_loc_random(uint64_t* lbuffer, size_t id, size_t test_size)
{
    std::mt19937_64 re((id + 1) * random_seed);
    std::uniform_int_distribution<uint64_t> dis(1,1ull<<63);

    for (size_t i = 0; i < test_size; ++i)
    {
        lbuffer[i] = dis(re);
    }
}

void fill_loc_insert(uint64_t* lbuffer, size_t id, size_t inserted, size_t test_size)
{
    for (size_t i = 0; i < test_size; ++i)
    {
        lbuffer[i] = keys[inserted + id*test_size + i];
    }
}


/* MAIN TEST FUNCTION *********************************************************/
template <class ThreadType, class FilterType>
int test(FilterType& filter, fill_benchmark_config config, size_t id)
{
    pin_to_core(id);

	const size_t input_n = config.n;
	config.n = qf::filter_capacity(config.n);

    size_t local_test_size = std::min(config.test_size / config.p, config.n / 10  / config.p);
    config.test_size = local_test_size * config.p;

    size_t stage = 0;
	const size_t max_steps = std::min(config.max_fill / 10, size_t{ 10 });


    if (ThreadType::is_main)
    {
        keys        = new uint64_t[2* config.n];
        keys_unsucc = &keys[config.n];
    }

    uint64_t *loc_buffer = new uint64_t[local_test_size];

	for (size_t i = 0; i < config.it * config.sub_it; ++i)
	{
		if (ThreadType::is_main)
		{
			++random_seed;
			current_block.store(0);
		}

		if (i % config.sub_it == 0)
			ThreadType::synchronized(generate_keys, ++stage, config.p-1, 2* config.n);

        // STAGE 0.3 - Initialize the quotient filter + some synchronization
		if (ThreadType::is_main)
			filter = construct_quotient_filter<FilterType>(config.bits, input_n);

		ThreadType::synchronized([] { return 0; }, ++stage, config.p - 1);

        size_t inserted = 0;
		//n = filter.capacity(); // use real cell count

        for ( size_t fill_step = 1; fill_step < max_steps; ++fill_step )
        {
            if (ThreadType::is_main) print_parameters(i, config.p, config.bits, config.n, config.test_size);
            auto next_test = size_t(fill_step * .1 * config.n);

            ThreadType::out(fill_step*10, 4);

            // STAGE 1 - Fills until next 10% bound
            {
                if (ThreadType::is_main) current_block.store(inserted);
                inserted = next_test;

                ThreadType::synchronized(insert<FilterType>, ++stage, config.p-1, filter, next_test);

                //filter.check_consistency();
            }

            // STAGE 2 - Contains for some randomly chosen inserted elements
            {
                //if (ThreadType::is_main) current_block.store(0);

                fill_loc_succ(loc_buffer, id, inserted, local_test_size);

                auto dur = ThreadType::synchronized(contains_successful<FilterType>, ++stage, config.p-1, filter, loc_buffer, local_test_size);

                ThreadType::out(dur.second/1000000., 10);
            }

            // STAGE 3 - Contains for non-inserted elements
            {
                fill_loc_random(loc_buffer, id, local_test_size);

                auto dur = ThreadType::synchronized(contains_unsuccessful<FilterType>, ++stage, config.p-1, filter, loc_buffer, local_test_size);

                ThreadType::out(dur.second/1000000., 10);
            }

            // STAGE 4 - Insert Test
            {
                fill_loc_insert(loc_buffer, id, inserted, local_test_size);

                auto dur = ThreadType::synchronized(insert_test<FilterType>, ++stage, config.p-1, filter, loc_buffer, local_test_size);

                ThreadType::out(dur.second/1000000., 10);
                inserted += config.test_size;
            }

            ThreadType::out(double(false_positives.load()) / config.test_size * 100, 10);
            ThreadType::out(errors.load(), 6);
            ThreadType() << std::endl;

            // End of Iteration Stuff + Some Synchronization
            if (ThreadType::is_main)
            {
                errors.store(0);
                false_positives.store(0);
            }

            ThreadType::synchronized([]{ return 0; }, ++stage, config.p-1);
        }
    }

    delete[] loc_buffer;

    if (ThreadType::is_main)
    {
        delete[] keys;
        reset_stages();
    }

    return 0;
}


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


/* MAIN FUNCTION: READS COMMANDLINE AND STARTS TEST THREADS *******************/
int main(int argn, char** argc)
{
	fill_benchmark_config config;
	CommandLine c{ argn, argc };

	config.p         = c.intArg("-p"         , config.p);
	config.n         = c.intArg("-n"         , config.n);
	config.test_size = c.intArg("-ts"        , config.test_size);
	config.it        = c.intArg("-it"        , config.it);
	config.sub_it    = c.intArg("-sub_it"    , config.sub_it);
	config.max_fill  = c.intArg("-max_fill_%", config.max_fill);
	config.set_fp_rate_percent(c.doubleArg("-fp_rate_%", config.fp_rate * 100));

	double fp_rate = c.doubleArg("-fp_rate_%", 0.0);
	size_t rbits = c.intArg("-r_bits", 0);

    // does not need to be specified anymore
	[[maybe_unused]] bool adjust_lp_fp_rate = c.boolArg("-adjust_lp_fp_rate");

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

    // // is now done automatically
	// if (QF_concurrency_variant == QF_concurrency_type::linear_probing && adjust_lp_fp_rate)
	// {
	// 	config.set_rbits(config.bits + 3);
	// }

    if (QF_type<key_type>::is_sequential) config.p = 1;

    print_column_names();

	auto start_test = [&](auto&& filter)
	{
		using QuotientFilterType = std::remove_reference_t<decltype(filter)>;
		start_threads(test<TimedMainThread, QuotientFilterType>, test<UnTimedSubThread, QuotientFilterType>, filter, config);
	};


    using qfilter_type = typename instanciation_handler<
        QF_type<key_type>,
        QF_type<key_type>::is_templated>::type;

    qfilter_type qfilter{};
    start_test(qfilter);

    return 0;
}
