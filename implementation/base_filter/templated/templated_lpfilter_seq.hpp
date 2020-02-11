#pragma once
#include <iostream>
#include <memory>
#include "utils/utilities.h"
#include "templated_lpfilter_cell.hpp"


namespace qf {

template <class Key,
          size_t remainder_bits,
          class Hash = default_hash>
class templated_lp_filter_seq
{
public:
    static constexpr bool is_growing_compatible = false;
    static constexpr bool is_templated  = true;
    static constexpr bool is_sequential = true;
    static constexpr bool is_growing    = false;
    static constexpr bool is_dynamic    = false;
    static constexpr bool uses_handle   = false;
    template <size_t rem>
    using instanciated = templated_lp_filter_seq<Key, rem, Hash>;

	using key_type = Key;
	using hash_function_type = Hash;

private:
	using cell_base_type = qf::grouped_cell_base_type;
	using cell_type = qf::templated_lpfilter_cell_non_atomic<cell_base_type,
                                                           remainder_bits>;
	using remainder_type = typename cell_type::remainder_type;
	using remainder_pos = qf::entry_position<cell_type::capacity>;

	static constexpr remainder_pos quotient_position(size_t quotient)
	{
		return { static_cast<std::uint32_t>(quotient / cell_type::capacity),
                 static_cast<std::uint8_t> (quotient % cell_type::capacity) };
	}

	alignas(128) size_t				                table_capacity;
	alignas(128) std::unique_ptr<cell_type[]>       table;
	alignas(128) hash_function_type                 hf;
	alignas(128) size_t                             quotient_bits;

public:

	templated_lp_filter_seq(size_t capacity = qf::DEFAULT_MIN_CAPACITY,
                          const hash_function_type& hf = hash_function_type())
	 : hf(hf), quotient_bits(qf::capacity_to_quotient_bits(capacity))
	{
		table_capacity = (ONE << quotient_bits) / cell_type::capacity
                         + qf::shift_buffer;
		table = std::make_unique<cell_type[]>(table_capacity);
	}

	templated_lp_filter_seq(const templated_lp_filter_seq& other)
		: table_capacity(other.table_capacity),
          table(new cell_type[table_capacity]),
          hf(other.hf),
          quotient_bits(other.quotient_bits)
	{
		for (size_t i = 0; i < table_capacity; i++)
		{
			table[i] = other.table[i];
		}
	}

	templated_lp_filter_seq& operator=(const templated_lp_filter_seq& other)
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

	templated_lp_filter_seq(templated_lp_filter_seq&& other) = default;
	templated_lp_filter_seq& operator=(
        templated_lp_filter_seq&& other) = default;

	bool insert(const key_type& key)
	{
		const auto[q, r] = this->get_quotient_and_remainder(key);
		const auto q_pos = quotient_position(q);
		auto iter = q_pos;
		cell_type current_cell = table[iter.table];

		while (!current_cell.is_empty(iter.cell))
		{
			if (current_cell.remainder(iter.cell) == r)
				return true;

			++iter;
			iter.table %= table_capacity;

			if (iter == q_pos)
				return false;

			current_cell = table[iter.table];
		}

		table[iter.table].set_remainder(r, iter.cell);
		return true;
	}

	bool contains(const key_type& key) const
	{
		const auto[q, r] = this->get_quotient_and_remainder(key);
		const auto q_pos = quotient_position(q);
		auto iter = q_pos;
		cell_type current_cell = table[iter.table];

		while (!current_cell.is_empty(iter.cell))
		{
			if (current_cell.remainder(iter.cell) == r)
				return true;

			++iter;
			iter.table %= table_capacity;

			if (iter == q_pos)
				return false;

			current_cell = table[iter.table];
		}

		return false;
	}

	std::pair<size_t, remainder_type>
    get_quotient_and_remainder(const key_type& key) const
	{
		auto[q, r] = qf::get_quotient_and_remainder<remainder_type>(hf(key),
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
		return table_capacity * sizeof(cell_type);
	}

	size_t unused_memory_bits() const
	{
		return table_capacity * (sizeof(cell_type) * 8
                                 - cell_type::capacity * remainder_bits);
	}

	double fill_level() const
	{
		size_t used_slots = 0;

		for (remainder_pos pos{ 0,0 }; pos.table < table_capacity; pos++)
		{
			if (!table[pos.table].is_empty(pos.cell))
				used_slots++;
		}

		return static_cast<double>(used_slots) / table_capacity
                                               / cell_type::capacity;
	}
};

} // namespace qf
