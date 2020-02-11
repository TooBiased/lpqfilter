#pragma once
#include <iostream>
#include <memory>

#include "utils/utilities.h"
#include "utils/definitions.h"


namespace qf
{

template < class Key,
	class Hash = default_hash>
	class nongrouped_lpfilter_conc
{
public:
	using key_type = Key;
	using hash_function_type = Hash;
private:
	using remainder_type = std::uint64_t;
	using cell = remainder_type;
	using cell_atomic = std::atomic<cell>;

	bool table_CAS(size_t pos, cell& expected, const cell& desired)
	{
		return table[pos].compare_exchange_strong(expected, desired);
	}

	alignas(128) size_t                           table_capacity;
	alignas(128) std::unique_ptr<cell_atomic[]>   table;
	alignas(128) hash_function_type               hf;
	alignas(128) size_t                           remainder_bits;
	alignas(128) size_t                           quotient_bits;
	alignas(128) size_t                           modulo_mask;

	static constexpr remainder_type one = static_cast<remainder_type>(1);

public:

	nongrouped_lpfilter_conc(size_t capacity = qf::DEFAULT_MIN_CAPACITY,
                             size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS,
                             const hash_function_type& hf =hash_function_type())
		: hf(hf),
          remainder_bits(remainder_bits),
          quotient_bits(qf::capacity_to_quotient_bits(capacity))
	{
		table_capacity = one << quotient_bits;
		table = std::make_unique<cell_atomic[]>(table_capacity);
		modulo_mask = table_capacity - 1;
	}

	nongrouped_lpfilter_conc(const nongrouped_lpfilter_conc& other)
		: table_capacity(other.table_capacity),
          table(new cell_atomic[table_capacity]),
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

	nongrouped_lpfilter_conc& operator=(const nongrouped_lpfilter_conc& other)
	{
		table_capacity = other.table_capacity;
		table = std::make_unique<cell_atomic[]>(other.table_capacity);
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

	nongrouped_lpfilter_conc(nongrouped_lpfilter_conc&& other) = default;
	nongrouped_lpfilter_conc& operator=(nongrouped_lpfilter_conc&&) = default;

	bool insert(const key_type& key)
	{
		const auto[q, r] = get_quotient_and_remainder(key);
		size_t pos = q;
		cell current = table[pos];

		do
		{
			while (current != 0)
			{
				if (current == r)
					return true;

				pos = (pos + 1) & modulo_mask;

				if (pos == q)
					return false;

				current = table[pos];
			}

		} while (!table_CAS(pos, current, r));
		return true;
	}

	bool contains(const key_type& key) const
	{
		const auto[q, r] = get_quotient_and_remainder(key);
		size_t pos = q;
		cell current = table[pos];

		while (current != 0)
		{
			if (current == r)
				return true;

			pos = (pos + 1) & modulo_mask;

			if (pos == q)
				return false;

			current = table[pos];
		}

		return false;
	}

	std::pair<size_t, remainder_type>
    get_quotient_and_remainder(const key_type& key) const
	{
		auto[q,r] = qf::get_quotient_and_remainder<remainder_type>(hf(key),
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
			if (table[pos] != 0)
				used_slots++;
		}

		return static_cast<double>(used_slots) / table_capacity;
	}
};

} // namespace qf
