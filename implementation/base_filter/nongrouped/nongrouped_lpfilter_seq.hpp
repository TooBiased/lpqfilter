#pragma once
/*******************************************************************************
 * implementation/base_filter/nongrouped/nongrouped_lpfilter_seq.hpp
 *
 * sequential linear probing filter with non-grouped slots (one slot per atomic)
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/

#include <iostream>
#include <memory>

#include "utils/default_hash.hpp"
namespace htm = utils_tm::hash_tm;

#include "implementation/utilities.hpp"
#include "implementation/definitions.hpp"

namespace qf
{

template < class Key,
           class Hash = htm::default_hash>
class nongrouped_lpfilter_seq
{
public:
	using key_type = Key;
	using hash_function_type = Hash;
private:
	using entry_type = std::uint64_t;
	using cell = entry_type;

	bool is_empty(size_t pos) const
	{
		return table[pos] == 0;
	}

	alignas(128) size_t                    table_capacity;
	alignas(128) std::unique_ptr<cell[]>   table;
	alignas(128) hash_function_type        hf;
	alignas(128) size_t                    remainder_bits;
	alignas(128) size_t                    quotient_bits;
	alignas(128) size_t                    modulo_mask;

	static constexpr entry_type            one = static_cast<entry_type>(1);

public:

	nongrouped_lpfilter_seq(size_t capacity = qf::DEFAULT_MIN_CAPACITY,
                            size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS,
                            const hash_function_type& hf = hash_function_type())
		: hf(hf),
          remainder_bits(remainder_bits+3),
          quotient_bits(qf::capacity_to_quotient_bits(capacity))
	{
		table_capacity = one << quotient_bits;
		table = std::make_unique<cell[]>(table_capacity);
		modulo_mask = table_capacity - 1;
	}

	nongrouped_lpfilter_seq(const nongrouped_lpfilter_seq& other)
		: table_capacity(other.table_capacity),
          table(new cell[table_capacity]),
          hf(other.hf),
          remainder_bits(other.remainder_bits),
          quotient_bits(other.quotient_bits),
          modulo_mask(other.modulo_mask)
	{
		for (size_t i = 0; i < table_capacity; i++)
		{
			table[i] = other.table[i];
		}
	}

	nongrouped_lpfilter_seq& operator=(const nongrouped_lpfilter_seq& other)
	{
		table_capacity = other.table_capacity;
		table = std::make_unique<cell[]>(other.table_capacity);
		hf = other.hf;
		remainder_bits = other.remainder_bits;
		quotient_bits = other.quotient_bits;
		modulo_mask = other.modulo_mask;

		for (size_t i = 0; i < table_capacity; i++)
		{
			table[i] = other.table[i];
		}

		return *this;
	}

	nongrouped_lpfilter_seq(nongrouped_lpfilter_seq&&) = default;
	nongrouped_lpfilter_seq& operator=(nongrouped_lpfilter_seq&&) = default;

	bool insert(const key_type& key)
	{
		const auto[q, r] = get_quotient_and_remainder(key);
		size_t pos = q;

		while (!is_empty(pos))
		{
			if (table[pos] == r)
				return true;

			pos = (pos + 1) & modulo_mask;

			if (pos == q)
				return false;
		}

		table[pos] = r;
		return true;
	}

	bool contains(const key_type& key) const
	{
		const auto[q, r] = get_quotient_and_remainder(key);
		size_t pos = q;

		while (!is_empty(pos))
		{
			if (table[pos] == r)
				return true;

			pos = (pos + 1) & modulo_mask;

			if (pos == q)
				return false;
		}

		return false;
	}

	std::pair<size_t, entry_type>
    get_quotient_and_remainder(const key_type& key) const
	{
		auto[q,r] = qf::get_quotient_and_remainder<entry_type>(hf(key),
                                                               quotient_bits,
                                                               remainder_bits);
		if (r == 0) r++; // 0 means empty slot and can't be used as remainder
		return { q,r };
	}

	size_t capacity() const
	{
		return ONE << quotient_bits;
	}


	size_t memory_usage_bytes() const
	{
		return table_capacity * sizeof(cell);
	}

	size_t unused_memory_bits() const
	{
		return table_capacity * (sizeof(cell) * 8 - remainder_bits);
	}

	double fill_level() const
	{
		size_t used_slots = 0;

		for (size_t pos = 0; pos < table_capacity; pos++)
		{
			if (!is_empty(pos))
				used_slots++;
		}

		return static_cast<double>(used_slots) / table_capacity;
	}
};

} // namespace qf
