#pragma once
#include <iostream>
#include <memory>
#include "utils/utilities.h"
#include "standard_lpfilter_cell.hpp"


namespace qf {

template <class Key,
          class Hash = default_hash>
class standard_lpfilter_seq
{
public:
    static constexpr bool is_growing_compatible = false;
    static constexpr bool is_templated  = false;
    static constexpr bool is_sequential = true;
    static constexpr bool is_growing    = false;
    static constexpr bool is_dynamic    = false;
    static constexpr bool uses_handle   = false;

	using key_type = Key;
	using hasher = Hash;

private:
	using cell_base_type = qf::grouped_cell_base_type;
	using cell_type = qf::standard_lpfilter_cell_non_atomic<cell_base_type>;
	using remainder_type = typename cell_type::entry_type;

    using entry_type     = typename cell_type::entry_type;
    using cell_stuff     = typename cell_type::cstuff;
    using entry_pos      = std::pair<size_t, size_t>;

    cell_stuff                    cc;
	size_t                        quotient_bits;
	size_t			              table_capacity;
	std::unique_ptr<cell_type[]>  table;
	hasher                        hf;

public:

	standard_lpfilter_seq(size_t capacity = qf::DEFAULT_MIN_CAPACITY,
                          size_t remainder = qf::DEFAULT_REMAINDER_BITS,
                          const hasher& hf = hasher())
        : cc(remainder),
          quotient_bits(qf::capacity_to_quotient_bits(capacity)),
          hf(hf)
	{
		table_capacity = (ONE << quotient_bits) / cc.capacity +qf::shift_buffer;
		table = std::make_unique<cell_type[]>(table_capacity);
	}

	standard_lpfilter_seq(const standard_lpfilter_seq& other)
		: table_capacity(other.table_capacity),
          table(new cell_type[table_capacity]), hf(other.hf),
          quotient_bits(other.quotient_bits)
	{
		for (size_t i = 0; i < table_capacity; i++)
		{
			table[i] = other.table[i];
		}
	}

	standard_lpfilter_seq& operator=(const standard_lpfilter_seq& other)
	{
		table_capacity = other.table_capacity;
		table = std::make_unique<cell_type[]>(other.table_capacity);
		hf = other.hf;
		quotient_bits = other.quotient_bits;

		for (size_t i = 0; i < table_capacity; i++)
		{
			table[i] = other.table[i];
		}

		return *this;
	}

	standard_lpfilter_seq(standard_lpfilter_seq&& other) = default;
	standard_lpfilter_seq& operator=(standard_lpfilter_seq&& other) = default;

	bool insert(const key_type& key)
	{
		const auto[q, r] = this->get_quotient_and_remainder(key);
		const auto q_pos = quotient_position(q);
		auto iter = q_pos;
		cell_type current_cell = table[iter.first];

		while (!current_cell.is_empty(cc, iter.second))
		{
			if (current_cell.remainder(cc, iter.second) == r)
				return true;

			increment(iter);

			if (iter.first >= table_capacity)
				return false;

			current_cell = table[iter.first];
		}

		table[iter.first].set_remainder(cc, r, iter.second);
		return true;
	}

	bool contains(const key_type& key) const
	{
		const auto[q, r] = this->get_quotient_and_remainder(key);
		const auto q_pos = quotient_position(q);
		auto iter = q_pos;
		cell_type current_cell = table[iter.first];

		while (!current_cell.is_empty(cc, iter.second))
		{
			if (current_cell.remainder(cc, iter.second) == r)
				return true;

			increment(iter);

			current_cell = table[iter.first];
		}

		return false;
	}

	std::pair<size_t, remainder_type>
    get_quotient_and_remainder(const key_type& key) const
	{
		auto[q, r] = qf::get_quotient_and_remainder<remainder_type>(hf(key),
                                                             quotient_bits,
                                                             cc.remainder_bits);
		if (r == 0) r++; // 0 means empty slot and can't be used as remainder
		return { q,r };
	}

	size_t capacity() const
	{
		return ONE << quotient_bits;
	}

	size_t memory_usage_bytes() const
	{
		return table_capacity * sizeof(cell_type);
	}

	size_t unused_memory_bits() const
	{
		return table_capacity * (sizeof(cell_type) * 8
                                 - cc.capacity * cc.remainder_bits);
	}

	double fill_level() const
	{
		size_t used_slots = 0;

		for (entry_pos pos{ 0,0 }; pos.first < table_capacity; increment(pos))
		{
			if (!table[pos.first].is_empty(cc, pos.second))
				used_slots++;
		}

		return static_cast<double>(used_slots) / table_capacity / cc.capacity;
	}

private:
    entry_pos quotient_position(const size_t quotient) const
    {
        return { quotient / cc.capacity, quotient % cc.capacity };
    }

    void increment(entry_pos& pos) const
    {
        pos.second++;
        if (pos.second == cc.capacity)
        {
            pos.second = 0;
            pos.first++;
        }
    }
};

} // namespace qf
