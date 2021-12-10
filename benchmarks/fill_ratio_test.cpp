/*******************************************************************************
 * benchmarks/fill_ratio_test.cpp
 *
 * This basic test fills a filter datastructure using p threads in regular
 * intervalls, the speed of the data structure is benchmarked. Inputs:
 *          -p       <# of threads>
 *          -n       <# of overall elements>
 *          -ts      <# test size>     // used to test at different fill degrees
 *          -max_fill_% <% where the algorithm stops>
 *          -fp_rate_%  <d% target false positive rate>
 *          -r_bits  <# of remainder bits>                //overwrites fp_rate_%
 *          -it      <# number of outer repetitions>          // different seeds
 *          -sub_it  <# number of inner iterations>                 // same seed
 *
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "benchmark_utils.hpp"
#include "filter_select.hpp"

#include "utils/command_line_parser.hpp"
#include "utils/pin_thread.hpp"
#include "utils/thread_coordination.hpp"

#include <iostream>
#include <random>
#include <tuple>


namespace utm = utils_tm;
namespace ttm = utils_tm::thread_tm;
namespace otm = utils_tm::out_tm;

template <typename T> inline void out(T t, size_t space)
{
    std::cout.width(space);
    std::cout << t << " " << std::flush;
}

void print_column_names()
{
    out("#i", 3);
    out("p", 3);
    out("bits", 7);
    out("n", 9);
    out("test_n", 6);
    out("f%", 4);
    out("t_cont_+", 10);
    out("t_cont_-", 10);
    out("t_insert", 10);
    out("fp_rate_%", 10);
    out("error", 6);
    std::cout << std::endl;
}

void print_parameters(size_t it, size_t p, size_t bits, size_t n, size_t tsize)
{
    out(it, 3);
    out(p, 3);
    out(bits, 7);
    out(n, 9);
    out(tsize, 6);
}


using key_type = size_t;

alignas(128) static std::atomic_size_t current_block;
/* no need */ static uint64_t*          keys;
/* no need */ static uint64_t*          keys_succ;
/* no need */ static uint64_t*          keys_unsucc;
/* no need */ static std::atomic_size_t errors;
/* no need */ static std::atomic_size_t false_positives;

std::atomic_size_t random_seed = 123;


/* STAGES OF THE MAIN TEST ****************************************************/
int generate_keys(size_t n)
{
    std::uniform_int_distribution<uint64_t> dis(1, 1ull << 63);

    ttm::execute_blockwise_parallel(current_block, n,
                                    [&dis](size_t s, size_t e) {
                                        std::mt19937_64 re(s * random_seed);

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
                                    [&dis](size_t s, size_t e) {
                                        std::mt19937_64 re(s * random_seed);

                                        for (size_t i = s; i < e; i++)
                                        {
                                            keys_succ[i] = keys[dis(re)];
                                        }
                                    });

    return 0;
}

template <class QuotientFilter> int insert(QuotientFilter& filter, size_t end)
{
    ttm::execute_parallel(current_block, end, [&filter](size_t i) {
        auto key = keys[i];

        if (!filter.insert(key)) { std::cerr << "i\n"; }
    });

    return 0;
}

template <class QuotientFilter>
int contains_successful(QuotientFilter& filter, uint64_t* lbuffer,
                        size_t test_size)
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

template <class QuotientFilter>
int contains_unsuccessful(QuotientFilter& filter, uint64_t* lbuffer,
                          size_t test_size)
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

template <class QuotientFilter>
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

void fill_loc_succ(uint64_t* lbuffer, size_t id, size_t inserted,
                   size_t test_size)
{
    std::mt19937_64                         re((id + 1) * random_seed);
    std::uniform_int_distribution<uint64_t> dis(0, inserted - 1);

    for (size_t i = 0; i < test_size; ++i) { lbuffer[i] = keys[dis(re)]; }
}

void fill_loc_random(uint64_t* lbuffer, size_t id, size_t test_size)
{
    std::mt19937_64                         re((id + 1) * random_seed);
    std::uniform_int_distribution<uint64_t> dis(1, 1ull << 63);

    for (size_t i = 0; i < test_size; ++i) { lbuffer[i] = dis(re); }
}

void fill_loc_insert(uint64_t* lbuffer, size_t id, size_t inserted,
                     size_t test_size)
{
    for (size_t i = 0; i < test_size; ++i)
    {
        lbuffer[i] = keys[inserted + id * test_size + i];
    }
}


/* MAIN TEST FUNCTION *********************************************************/
template <class FilterType> struct test
{
    template <class ThreadType> struct test_filter
    {
        static int
        execute(ThreadType t, FilterType& filter, fill_benchmark_config config)
        {
            utm::pin_to_core(t.id);

            const size_t input_n = config.n;
            config.n             = qf::filter_capacity(config.n);

            size_t local_test_size =
                std::min(config.test_size / config.p, config.n / 10 / config.p);
            config.test_size = local_test_size * config.p;

            const size_t max_steps = std::min(config.max_fill / 10, size_t{10});


            if (ThreadType::is_main)
            {
                keys        = new uint64_t[2 * config.n];
                keys_unsucc = &keys[config.n];
            }

            uint64_t* loc_buffer = new uint64_t[local_test_size];

            for (size_t i = 0; i < config.it * config.sub_it; ++i)
            {
                if (ThreadType::is_main)
                {
                    ++random_seed;
                    current_block.store(0);
                }

                if (i % config.sub_it == 0)
                    t.synchronized(generate_keys, 2 * config.n);

                // STAGE 0.3 - Initialize the quotient filter + some
                // synchronization
                if (ThreadType::is_main)
                    filter = construct_quotient_filter<FilterType>(config.bits,
                                                                   input_n);

                t.synchronize();

                size_t inserted = 0;
                // n = filter.capacity(); // use real cell count

                for (size_t fill_step = 1; fill_step < max_steps; ++fill_step)
                {
                    if (ThreadType::is_main)
                        print_parameters(i, config.p, config.bits, config.n,
                                         config.test_size);
                    auto next_test = size_t(fill_step * .1 * config.n);

                    t.out << otm::width(4) << fill_step * 10;

                    // STAGE 1 - Fills until next 10% bound
                    {
                        if (ThreadType::is_main) current_block.store(inserted);
                        inserted = next_test;

                        t.synchronized(insert<FilterType>, filter, next_test);

                        // filter.check_consistency();
                    }

                    // STAGE 2 - Contains for some randomly chosen inserted
                    // elements
                    {
                        // if (ThreadType::is_main) current_block.store(0);

                        fill_loc_succ(loc_buffer, t.id, inserted,
                                      local_test_size);

                        auto dur =
                            t.synchronized(contains_successful<FilterType>,
                                           filter, loc_buffer, local_test_size);

                        t.out << otm::width(10) << dur.second / 1000000.;
                    }

                    // STAGE 3 - Contains for non-inserted elements
                    {
                        fill_loc_random(loc_buffer, t.id, local_test_size);

                        auto dur =
                            t.synchronized(contains_unsuccessful<FilterType>,
                                           filter, loc_buffer, local_test_size);

                        t.out << otm::width(10) << dur.second / 1000000.;
                    }

                    // STAGE 4 - Insert Test
                    {
                        fill_loc_insert(loc_buffer, t.id, inserted,
                                        local_test_size);

                        auto dur =
                            t.synchronized(insert_test<FilterType>, filter,
                                           loc_buffer, local_test_size);

                        t.out << otm::width(10) << dur.second / 1000000.;
                        inserted += config.test_size;
                    }

                    t.out << otm::width(10)
                          << double(false_positives.load()) / config.test_size *
                                 100
                          << otm::width(6) << errors.load() << std::endl;

                    // End of Iteration Stuff + Some Synchronization
                    if (ThreadType::is_main)
                    {
                        errors.store(0);
                        false_positives.store(0);
                    }

                    t.synchronize();
                }
            }

            delete[] loc_buffer;

            if (ThreadType::is_main) { delete[] keys; }

            return 0;
        }
    };
};

template <class Type, bool instanciated> class instanciation_handler
{
  public:
    using type = Type;
};

template <class Type> class instanciation_handler<Type, true>
{
  public:
    using type =
        typename Type::template instanciated<qf::DEFAULT_REMAINDER_BITS>;
};


/* MAIN FUNCTION: READS COMMANDLINE AND STARTS TEST THREADS *******************/
int main(int argn, char** argc)
{
    fill_benchmark_config    config;
    utm::command_line_parser c{argn, argc};

    config.p         = c.int_arg("-p", config.p);
    config.n         = c.int_arg("-n", config.n);
    config.test_size = c.int_arg("-ts", config.test_size);
    config.it        = c.int_arg("-it", config.it);
    config.sub_it    = c.int_arg("-sub_it", config.sub_it);
    config.max_fill  = c.int_arg("-max_fill_%", config.max_fill);
    config.set_fp_rate_percent(
        c.double_arg("-fp_rate_%", config.fp_rate * 100));

    double fp_rate = c.double_arg("-fp_rate_%", 0.0);
    size_t rbits   = c.int_arg("-r_bits", 0);

    if (rbits != 0 && fp_rate != 0.0)
    {
        std::cout << "Error: only one of \"-fp_rate_%\" and \"-r_bits\" can be "
                     "used\n";
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

    if (QF_type<key_type>::is_sequential) config.p = 1;

    print_column_names();

    using qfilter_type =
        typename instanciation_handler<QF_type<key_type>,
                                       QF_type<key_type>::is_templated>::type;

    qfilter_type qfilter{};
    // start_test(qfilter);
    ttm::start_threads<test<qfilter_type>::template test_filter>(
        config.p, qfilter, config);

    return 0;
}
