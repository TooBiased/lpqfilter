#include "filter_select.hpp"
#include "benchmark_utils.hpp"

#include "utils/pin.h"
#include "utils/test_coordination.h"

#include "implementation/composite_filter/growing_qfilter.hpp"
#include "implementation/composite_filter/dynamic_qfilter.hpp"

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
    out("#i"              , 3);
    out("p"               , 3);
    out("n"               , 9);
	out("target_fp_rate_%", 17);
	out("min_cap"         , 10);
    out("sample"          , 7);
    out("t_insert"        , 10);
    out("t_cont_+"        , 10);
    out("t_cont_-"        , 10);
    out("fp_rate_%"       , 10);
    out("error"           , 6);
    std::cout << std::endl;
}

void print_parameters(size_t it, size_t p, double fp_rate, size_t n, size_t cap, size_t sample)
{
    out(it           , 3);
    out(p            , 3);
    out(n            , 9);
	out(fp_rate * 100, 17);
	out(cap          , 10);
	out(sample       , 7);
}

using key_type = size_t;

alignas(128)  static std::atomic_size_t  current_block;
/* no need */ static uint64_t*           keys;
/* no need */ static uint64_t*           keys_succ;
/* no need */ static uint64_t*           keys_unsucc;
/* no need */ static std::atomic_size_t  errors;
/* no need */ static std::atomic_size_t  false_positives;

std::atomic_size_t offset;
std::atomic_size_t random_seed = 123;


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

void fill_succ(size_t max_index, size_t count, size_t sample)
{
	std::uniform_int_distribution<size_t> dis(0, max_index - 1);
	std::mt19937_64 re(sample * random_seed + max_index);

	for (size_t i = 0; i < count; i++)
	{
		keys_succ[i] = keys[dis(re)];
	}
}

template<class QuotientFilter>
int insert(QuotientFilter& filter, size_t end)
{
	const size_t off = offset;
	execute_parallel(current_block, end,
					 [&filter, off](size_t i)
	{
		auto key = keys[off + i];
 		if (!filter.insert(key))
		{
			std::cerr << "i\n";
		}
	});

	return 0;
}

template<class QuotientFilter>
int contains_successful(QuotientFilter& filter, size_t end)
{
	const size_t off = offset;
	size_t err = 0;
	execute_parallel(current_block, end,
					 [&filter, &err, off](size_t i)
	{
		auto key = keys_succ[i];
		if (!filter.contains(key))
		{
			err++;
			std::cerr << "c";
		}
	});

	errors.fetch_add(err, std::memory_order_relaxed);
	return 0;
}

template<class QuotientFilter>
int contains_unsuccessful(QuotientFilter& filter, size_t end)
{
	const size_t off = offset;
	size_t fp = 0;
	execute_parallel(current_block, end,
					 [&filter, &fp, off](size_t i)
	{
		auto key = keys_unsucc[off + i];
		if (filter.contains(key))
			fp++;
	});

	false_positives.fetch_add(fp, std::memory_order_relaxed);
	return 0;
}




/* MAIN TEST FUNCTION *********************************************************/
template <class ThreadType, class FilterType>
int test(FilterType& filter, progression_benchmark_config config, size_t id)
{
	using HandleType = qf::HandleType<FilterType>;
	pin_to_core(id);

	size_t stage = 0;
	const size_t sample_size = config.n / config.sample_count;

	if (ThreadType::is_main)
	{
		keys        = new uint64_t[2 * config.n];
		keys_succ   = new uint64_t[sample_size];
		keys_unsucc = &keys[config.n];
	}

	for (size_t i = 0; i < config.it * config.sub_it; ++i)
	{
		if (ThreadType::is_main)
		{
			++random_seed;
			current_block.store(0);
		}

		if (i % config.sub_it == 0)
			ThreadType::synchronized(generate_keys, ++stage, config.p - 1, 2 * config.n);

		// STAGE 0.3 - Reinitialize the quotient filter + some synchronization
		if (ThreadType::is_main)
			filter = construct_quotient_filter<FilterType>(config);

		ThreadType::synchronized([] { return 0; }, ++stage, config.p - 1);

		auto handle = create_handle(filter);
		auto start_cap = filter.capacity(); // use real cell count

		for (size_t sample = 0; sample < config.sample_count; sample++)
		{
			if (ThreadType::is_main)
			{
				errors = 0;
				false_positives = 0;
				offset = sample * sample_size;
				print_parameters(i, config.p, config.fp_rate, config.n, start_cap, sample);
			}

			ThreadType::synchronized([] { return 0; }, ++stage, config.p - 1);

			// STAGE 1 - Insert
			{
				if (ThreadType::is_main)
					current_block.store(0);

				auto dur = ThreadType::synchronized(insert<HandleType>, ++stage, config.p - 1, handle, sample_size);
				ThreadType::out(dur.second / 1000000., 10);
			}

			// STAGE 2 - Successful Contains
			{
				if (ThreadType::is_main) {
					fill_succ(offset + sample_size, sample_size, sample);
					current_block.store(0);
				}

				ThreadType::synchronized([] { return 0; }, ++stage, config.p - 1);

				auto dur = ThreadType::synchronized(contains_successful<HandleType>, ++stage, config.p - 1, handle, sample_size);
				ThreadType::out(dur.second / 1000000., 10);
			}

			// STAGE 3 - Unsuccessful Contains
			{
				if (ThreadType::is_main)
					current_block.store(0);

				auto dur = ThreadType::synchronized(contains_unsuccessful<HandleType>, ++stage, config.p - 1, handle, sample_size);
				ThreadType::out(dur.second / 1000000., 10);
				ThreadType::out(double(false_positives.load()) / sample_size * 100, 10);
				ThreadType::out(errors.load(), 6);
			}

			// End of Iteration Stuff + Some Synchronization
			ThreadType() << std::endl;

			ThreadType::synchronized([] { return 0; }, ++stage, config.p - 1);
		}

		ThreadType::synchronized([] { return 0; }, ++stage, config.p - 1);
	}

	if (ThreadType::is_main)
	{
		delete[] keys;
		delete[] keys_succ;
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
	progression_benchmark_config config;
    CommandLine c{argn, argc};

	config.growing       = c.boolArg("-growing");
	config.dynamic       = c.boolArg("-dynamic");
	config.p             = c.intArg("-p"         , config.p);
    config.n             = c.intArg("-n"         , config.n);
	config.cap           = c.intArg("-cap"       , config.n * 2);
	config.min_cap       = c.intArg("-min_cap"   , config.n / 100);
    config.it            = c.intArg("-it"        , config.it);
	config.sub_it        = c.intArg("-sub_it"    , config.sub_it);
	config.sample_count  = c.intArg("-samples"   , config.sample_count);

	double fp_rate = c.doubleArg("-fp_rate_%", 0.0);
	size_t rbits = c.intArg("-r_bits", 0);
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

	// if (QF_concurrency_variant == QF_concurrency_type::linear_probing && adjust_lp_fp_rate)
	// {
	// 	config.set_rbits(config.bits + 3);
	// }

	if (config.growing && config.dynamic)
	{
		std::cout << "Error: only one of \"-growing\" or \"-dynamic\" can be used\n";
		return 1;
	}

	print_column_names();

	auto start_test = [&](auto&& filter)
	{
		using QuotientFilterType = std::remove_reference_t<decltype(filter)>;
		start_threads(test<TimedMainThread, QuotientFilterType>, test<UnTimedSubThread, QuotientFilterType>, filter, config);
	};

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
        start_test(qfilter);
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
        start_test(qfilter);
	}
	else
	{
        using qfilter_type = typename instanciation_handler<
            QF_type<key_type>,
            QF_type<key_type>::is_templated>::type;
            // typename std::conditional
            // <
            //     QF_type<key_type>::is_templated,
            //     typename QF_type<key_type>::template instanciated<qf::DEFAULT_REMAINDER_BITS> ,
            //     QF_type<key_type>
            // >::type;

        qfilter_type qfilter{};
        start_test(qfilter);
	}

    return 0;
}
