#pragma once

#include "utils/utilities.h"
#include "utils/commandline.h"

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
	size_t cap = 2000000;
	size_t min_cap = 250000;
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

// struct mixed_benchmark_config : public basic_benchmark_config
// {
// 	size_t contains_ratio = 2;
// };

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







// template <class key_type, size_t... Bits, class Functor>
// void run_test_compact_qf_helper(size_t remainder_bits, Functor&& func, std::index_sequence<Bits...>)
// {
// 	((remainder_bits == Bits && (
// 		func(QF_type<key_type, Bits>{}),
// 		true)) || ...);
// }

// template<size_t First, size_t... Bits>
// std::index_sequence<Bits...> remove_first(std::index_sequence<First, Bits...>)
// {
// 	return {};
// };

// template<class key_type, class Functor>
// void run_test_compact_qf(size_t remainder_bits, Functor&& func)
// {
// 	if constexpr (QF_concurrency_variant == QF_concurrency_type::linear_probing)
// 	{
// 		run_test_compact_qf_helper<key_type>(remainder_bits,
//                                              std::forward<Functor>(func),
//                                              remove_first(qf::REMAINDER_SEQUENCE));
// 	}
// 	else
// 	{
// 		run_test_compact_qf_helper<key_type>(remainder_bits,
//                                              std::forward<Functor>(func),
//                                              qf::REMAINDER_SEQUENCE);
// 	}
// }
