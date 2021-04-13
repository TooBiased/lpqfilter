#pragma once
/*******************************************************************************
 * benchmarks/benchmark_utils.cpp
 *
 * This file contains the base configuration of all benchmarks
 * + some common functions
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include "implementation/utilities.hpp"
#include "utils/command_line_parser.hpp"

struct benchmark_config
{
	size_t p = 1;
	size_t it = 3;
	size_t sub_it = 3;
	size_t n = 1'000'000;
	size_t bits = 13;
};

struct basic_benchmark_config : public benchmark_config
{
	size_t cap = 2'000'000;
	size_t min_cap = 250'000;
	double fp_rate = 0.0001;
	bool growing = false;
	bool dynamic = false;

	void set_rbits(size_t rbits)
	{
		bits = rbits;
		fp_rate = qf::remainder_bits_to_false_positive_rate(bits);
	}

	void set_fp_rate(double fp_rate)
	{
		this->fp_rate = fp_rate;
		bits = qf::false_positive_rate_to_remainder_bits(fp_rate);
	}

	void set_fp_rate_percent(double fp_rate_percent)
	{
		set_fp_rate(fp_rate_percent / 100);
	}
};

struct fill_benchmark_config : public basic_benchmark_config
{
	size_t test_size = 200000;
	size_t max_fill = 80;
};

struct progression_benchmark_config : public basic_benchmark_config
{
	size_t sample_count = 20;
};


template<class FilterType>
FilterType construct_quotient_filter(size_t remainder_bits, size_t capacity)
{
	if constexpr (FilterType::is_templated)
        return FilterType{ capacity };
    else
        return FilterType{ capacity, remainder_bits };
}

template<class FilterType>
FilterType construct_quotient_filter(size_t remainder_bits, size_t capacity, size_t min_capacity, double fp_rate)
{
	if constexpr      (  FilterType::is_dynamic  )
	{
		return FilterType{ fp_rate, min_capacity };
	}
	else if constexpr (  FilterType::is_growing  )
	{
		return FilterType{ min_capacity, remainder_bits };
	}
    else if constexpr (  FilterType::is_templated  )
	{
		return FilterType{ capacity };
	}
    else if constexpr (! FilterType::is_templated  )
	{
		return FilterType{ capacity, remainder_bits };
	}

	return {};
}

template<class FilterType>
FilterType construct_quotient_filter(const basic_benchmark_config& config)
{
	return construct_quotient_filter<FilterType>(config.bits, config.cap, config.min_cap, config.fp_rate);
}
