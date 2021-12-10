#pragma once
/*******************************************************************************
 * implementation/base_filter/standard/standard_lpfilter_conc.hpp
 *
 * concurrent linear probing filter with grouped slots
 * (multiple slots per data slot non-templated)
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
#include "standard_lpfilter_cell.hpp"


namespace qf
{

template <class Key, class Hash = htm::default_hash>
class standard_lpfilter_conc
{
  public:
    static constexpr bool is_growing_compatible = false;
    static constexpr bool is_templated          = false;
    static constexpr bool is_sequential         = false;
    static constexpr bool is_growing            = false;
    static constexpr bool is_dynamic            = false;
    static constexpr bool uses_handle           = false;

    using key_type           = Key;
    using hash_function_type = Hash;

  private:
    using this_type      = standard_lpfilter_conc;
    using cell_base_type = qf::grouped_cell_base_type;
    using cell_type   = qf::standard_lpfilter_cell_non_atomic<cell_base_type>;
    using cell_atomic = typename cell_type::atomic;

    using entry_type = typename cell_type::entry_type;
    using cell_stuff = typename cell_type::cstuff;
    using entry_pos  = std::pair<size_t, size_t>;

    cell_stuff cc;
    alignas(128) size_t quotient_bits;
    alignas(128) size_t table_capacity;
    alignas(128) std::unique_ptr<cell_atomic[]> table;
    alignas(128) hash_function_type hf;

  public:
    standard_lpfilter_conc(size_t capacity       = qf::DEFAULT_MIN_CAPACITY,
                           size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS,
                           const hash_function_type& hf = hash_function_type())
        : cc(remainder_bits + 3),
          quotient_bits(qf::capacity_to_quotient_bits(capacity)), hf(hf)
    {
        table_capacity =
            (ONE << quotient_bits) / cc.capacity + qf::shift_buffer;
        table = std::make_unique<cell_atomic[]>(table_capacity);
    }

    standard_lpfilter_conc(const standard_lpfilter_conc& other)
        : cc(other.cc), quotient_bits(other.quotient_bits),
          table_capacity(other.table_capacity),
          table(new cell_atomic[table_capacity]), hf(other.hf)
    {
        for (size_t i = 0; i < table_capacity; i++)
        {
            table[i] = other.table[i];
        }
    }

    standard_lpfilter_conc& operator=(const standard_lpfilter_conc& other)
    {
        if (&other == this) return *this;

        this->~this_type();
        new (this) this_type(other);
        return *this;
    }

    standard_lpfilter_conc(standard_lpfilter_conc&& other) = default;
    standard_lpfilter_conc& operator=(standard_lpfilter_conc&& other) = default;

    inline bool insert(const key_type& key)
    {
        const auto [q, r]     = this->get_quotient_and_remainder(key);
        const entry_pos q_pos = quotient_position(q);
        entry_pos       iter  = q_pos;

        cell_type current_cell = table[iter.first];
        cell_type inserted_cell;

        do {
            while (!current_cell.is_empty(cc, iter.second))
            {
                if (current_cell.remainder(cc, iter.second) == r) return true;

                increment(iter);

                if (iter.first >= table_capacity) return false;

                current_cell = table[iter.first];
            }
            inserted_cell = current_cell;
            inserted_cell.set_remainder(cc, r, iter.second);

        } while (!table[iter.first].CAS(current_cell, inserted_cell));

        return true;
    }

    inline bool contains(const key_type& key) const
    {
        const auto [q, r]       = this->get_quotient_and_remainder(key);
        const auto q_pos        = quotient_position(q);
        auto       iter         = q_pos;
        cell_type  current_cell = table[iter.first];

        while (!current_cell.is_empty(cc, iter.second))
        {
            if (current_cell.remainder(cc, iter.second) == r) return true;

            increment(iter);

            current_cell = table[iter.first];
        }

        return false;
    }

    inline std::pair<size_t, size_t>
    get_quotient_and_remainder(const key_type& key) const
    {
        auto [q, r] = qf::get_quotient_and_remainder<size_t>(
            hf(key), quotient_bits, cc.remainder_bits);
        if (r == 0) r++; // 0 means empty slot and can't be used as remainder
        return {q, r};
    }

    inline size_t capacity() const { return ONE << quotient_bits; }

    inline size_t memory_usage_bytes() const
    {
        return table_capacity * sizeof(cell_type);
    }

    inline size_t unused_memory_bits() const
    {
        return table_capacity *
               (sizeof(cell_type) * 8 - cc.capacity * cc.remainder_bits);
    }

    inline double fill_level() const
    {
        size_t used_slots = 0;

        for (entry_pos pos{0, 0}; pos.first < table_capacity; increment(pos))
        {
            if (!table[pos.first].is_empty(cc, pos.second)) used_slots++;
        }

        return static_cast<double>(used_slots) / table_capacity / cc.capacity;
    }

  private:
    inline entry_pos quotient_position(const size_t quotient) const
    {
        return {quotient / cc.capacity, quotient % cc.capacity};
    }

    inline void increment(entry_pos& pos) const
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
