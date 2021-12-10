#pragma once
/*******************************************************************************
 * implementation/base_filter/templated/templated_lpfilter_conc.hpp
 *
 * concurrent linear probing filter with templated grouped slots
 * (multiple slots per data element templated)
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
#include "templated_lpfilter_cell.hpp"

namespace qf
{

template <class Key, size_t RemainderBits, class Hash = htm::default_hash>
class templated_lpfilter_conc
{
  public:
    static constexpr bool   is_growing_compatible = false;
    static constexpr bool   is_templated          = true;
    static constexpr bool   is_sequential         = false;
    static constexpr bool   is_growing            = false;
    static constexpr bool   is_dynamic            = false;
    static constexpr bool   uses_handle           = false;
    static constexpr size_t remainder_bits        = RemainderBits + 3;

    template <size_t rem>
    using instanciated = templated_lpfilter_conc<Key, rem, Hash>;

    using key_type           = Key;
    using hash_function_type = Hash;

  private:
    using cell_base_type = qf::grouped_cell_base_type;
    using cell_type =
        qf::templated_lpfilter_cell_non_atomic<cell_base_type, remainder_bits>;
    using cell_atomic    = typename cell_type::atomic;
    using remainder_type = typename cell_type::remainder_type;
    using remainder_pos  = qf::entry_position<cell_type::capacity>;

    static constexpr remainder_pos quotient_position(size_t quotient)
    {
        return {static_cast<std::uint32_t>(quotient / cell_type::capacity),
                static_cast<std::uint8_t>(quotient % cell_type::capacity)};
    }

    alignas(128) size_t table_capacity;
    alignas(128) std::unique_ptr<cell_atomic[]> table;
    alignas(128) hash_function_type hf;
    alignas(128) size_t quotient_bits;

  public:
    templated_lpfilter_conc(size_t capacity = qf::DEFAULT_MIN_CAPACITY,
                            const hash_function_type& hf = hash_function_type())
        : hf(hf), quotient_bits(qf::capacity_to_quotient_bits(capacity))
    {
        table_capacity =
            (ONE << quotient_bits) / cell_type::capacity + qf::shift_buffer;
        table = std::make_unique<cell_atomic[]>(table_capacity);
    }

    templated_lpfilter_conc(const templated_lpfilter_conc& other)
        : table_capacity(other.table_capacity),
          table(new cell_atomic[table_capacity]), hf(other.hf),
          quotient_bits(other.quotient_bits)
    {
        for (size_t i = 0; i < table_capacity; i++)
        {
            table[i] = other.table[i];
        }
    }

    templated_lpfilter_conc& operator=(const templated_lpfilter_conc& other)
    {
        table_capacity = other.table_capacity;
        table          = std::make_unique<cell_atomic[]>(other.table_capacity);
        hf             = other.hf;
        quotient_bits  = other.quotient_bits;

        for (size_t i = 0; i < table_capacity; i++)
        {
            table[i] = other.table[i];
        }

        return *this;
    }

    templated_lpfilter_conc(templated_lpfilter_conc&& other) = default;
    templated_lpfilter_conc&
    operator=(templated_lpfilter_conc&& other) = default;

    bool insert(const key_type& key)
    {
        const auto [q, r] = this->get_quotient_and_remainder(key);
        const auto q_pos  = quotient_position(q);
        auto       iter   = q_pos;

        cell_type current_cell = table[iter.table];
        cell_type inserted_cell;

        do {
            while (!current_cell.is_empty(iter.cell))
            {
                if (current_cell.remainder(iter.cell) == r) return true;

                ++iter;
                iter.table %= table_capacity;

                if (iter == q_pos) return false;

                current_cell = table[iter.table];
            }

            inserted_cell = current_cell;
            inserted_cell.set_remainder(r, iter.cell);

        } while (!table[iter.table].CAS(current_cell, inserted_cell));

        return true;
    }

    bool contains(const key_type& key) const
    {
        const auto [q, r]       = this->get_quotient_and_remainder(key);
        const auto q_pos        = quotient_position(q);
        auto       iter         = q_pos;
        cell_type  current_cell = table[iter.table];

        while (!current_cell.is_empty(iter.cell))
        {
            if (current_cell.remainder(iter.cell) == r) return true;

            ++iter;
            iter.table %= table_capacity;

            if (iter == q_pos) return false;

            current_cell = table[iter.table];
        }

        return false;
    }

    std::pair<size_t, remainder_type>
    get_quotient_and_remainder(const key_type& key) const
    {
        auto [q, r] = qf::get_quotient_and_remainder<remainder_type>(
            hf(key), quotient_bits, remainder_bits);
        if (r == 0) r++; // 0 means empty slot and can't be used as remainder
        return {q, r};
    }

    size_t capacity() const { return ONE << quotient_bits; }

    size_t memory_usage_bytes() const
    {
        return table_capacity * sizeof(cell_type);
    }

    size_t unused_memory_bits() const
    {
        return table_capacity *
               (sizeof(cell_type) * 8 - cell_type::capacity * remainder_bits);
    }

    double fill_level() const
    {
        size_t used_slots = 0;

        for (remainder_pos pos{0, 0}; pos.table < table_capacity; pos++)
        {
            if (!table[pos.table].is_empty(pos.cell)) used_slots++;
        }

        return static_cast<double>(used_slots) / table_capacity /
               cell_type::capacity;
    }
};


} // namespace qf
