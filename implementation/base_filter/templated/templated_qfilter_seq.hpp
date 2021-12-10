#pragma once
/*******************************************************************************
 * implementation/base_filter/templated/templated_qfilter_seq.hpp
 *
 * sequential quotient filter with templated grouped slots
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
#include "templated_qfilter_cell.hpp"

namespace qf
{

template <class Key, size_t remainder_bits, class Hash>
class templated_qfilter_seq;

template <class Key, class Hash = htm::default_hash>
class templated_qfilter_seq_base
{
  public:
    using key_type           = Key;
    using hash_function_type = Hash;
    using hashed_type = decltype(std::declval<Hash>()(std::declval<Key>()));

    template <size_t new_remainder>
    using instanciated = templated_qfilter_seq<Key, new_remainder, Hash>;
    static constexpr bool is_growing_compatible = true;
    static constexpr bool is_templated          = true;
    static constexpr bool is_sequential         = true;
    static constexpr bool is_growing            = false;
    static constexpr bool is_dynamic            = false;
    static constexpr bool uses_handle           = false;


    virtual ~templated_qfilter_seq_base() = default;

    virtual void               prefetch(const hashed_type& hashed) const    = 0;
    virtual qf::InsertResult   insert(const key_type& key)                  = 0;
    virtual qf::InsertResult   insert_hash(const hashed_type& hashed)       = 0;
    virtual bool               quick_insert(const key_type& key)            = 0;
    virtual bool               quick_insert_hash(const hashed_type& hashed) = 0;
    virtual qf::ContainsResult contains(const key_type& key) const          = 0;
    virtual qf::ContainsResult
                   contains_hash(const hashed_type& hashed) const = 0;
    virtual size_t capacity() const                               = 0;
    virtual size_t memory_usage_bytes() const                     = 0;
    virtual size_t unused_memory_bits() const                     = 0;
    virtual double fill_level() const                             = 0;
    virtual bool   check_consistency()                            = 0;
    virtual templated_qfilter_seq_base* create_bigger_QF(
        const hash_function_type& hf = hash_function_type()) const        = 0;
    virtual void grow(templated_qfilter_seq_base<Key, Hash>& target_base) = 0;
};


template <class Key, size_t remainder_bits, class Hash = htm::default_hash>
class templated_qfilter_seq : public templated_qfilter_seq_base<Key, Hash>
{
  public:
    template <size_t new_remainder>
    using other_remainder = templated_qfilter_seq<Key, new_remainder, Hash>;
    static constexpr bool is_growing_compatible = true;
    static constexpr bool is_templated          = true;
    static constexpr bool is_sequential         = true;
    static constexpr bool is_growing            = false;
    static constexpr bool is_dynamic            = false;
    static constexpr bool uses_handle           = false;

    using key_type           = Key;
    using hash_function_type = Hash;
    using hashed_type = decltype(std::declval<Hash>()(std::declval<Key>()));
    using next_smaller_qf =
        templated_qfilter_seq<Key, remainder_bits + 1, Hash>;
    using next_bigger_qf = templated_qfilter_seq<Key, remainder_bits - 1, Hash>;
    using base_qf        = templated_qfilter_seq_base<Key, Hash>;

    friend base_qf;
    friend next_smaller_qf;

  private:
    using cell_base_type = qf::grouped_cell_base_type;
    using cell_type =
        qf::templated_cell_non_atomic<cell_base_type, remainder_bits>;
    using entry_type = typename cell_type::entry_type;
    using entry_pos  = qf::entry_position<cell_type::capacity>;

    entry_pos find_cluster_start(entry_pos pos) const
    {
        while (table[pos.table].is_shifted(pos.cell)) { --pos; }

        return pos;
    }

    entry_pos find_run_start(entry_pos cluster_start, entry_pos q) const
    {
        auto bucket = cluster_start;
        auto run    = bucket;

        while (bucket != q)
        {
            do {
                if (run.table == table_capacity - 1)
                {
                    return run; // invalid value
                }

                ++run;
            } while (table[run.table].is_continuation(run.cell));

            do {
                ++bucket;
            } while (!table[bucket.table].is_occupied(bucket.cell));
        }

        return run;
    }

    // assumes that pos is not empty
    bool has_space_for_shift(entry_pos pos) const
    {
        while (!table[pos.table].is_empty(pos.cell)) { ++pos; }

        return pos.table != table_capacity - 1;
    }

    // inserts entry at pos and shifts the following entries
    // cell at pos must have enough space for the shift i.e. an empty entry
    void insert_entry_not_full(entry_type entry, entry_pos insert_pos,
                               size_t empty_pos)
    {
        const cell_type cell = table[insert_pos.table];
        cell_type       shifted_cell;

        // write unshifted entries at beginning of cell
        for (size_t i = 0; i < insert_pos.cell; i++)
        {
            shifted_cell.set_entry(cell.entry(i), i);
        }

        // write new entry
        shifted_cell.set_entry(entry, insert_pos.cell);

        // write shifted entries
        for (size_t i = insert_pos.cell + 1; i <= empty_pos; i++)
        {
            shifted_cell.set_entry(cell.entry(i - 1), i);
            shifted_cell.set_shifted(true, i);
        }

        // write unshifted entries at end of cell
        for (size_t i = empty_pos + 1; i < cell_type::capacity; i++)
        {
            shifted_cell.set_entry(cell.entry(i), i);
        }

        shifted_cell.set_occupied_full(cell.occupied_mask_full());
        table[insert_pos.table] = shifted_cell;
    }

    // inserts entry at pos and shifts the following entries
    // the last entry which is lost due to the shift is written back to the
    // entry reference
    void insert_entry_full(entry_type& entry, entry_pos insert_pos)
    {
        const cell_type cell = table[insert_pos.table];
        cell_type       shifted_cell;

        // write unshifted entries at beginning of cell
        for (size_t i = 0; i < insert_pos.cell; i++)
        {
            shifted_cell.set_entry(cell.entry(i), i);
        }

        // write new entry
        shifted_cell.set_entry(entry, insert_pos.cell);

        // write shifted entries
        for (size_t i = insert_pos.cell + 1; i < cell_type::capacity; i++)
        {
            shifted_cell.set_entry(cell.entry(i - 1), i);
            shifted_cell.set_shifted(true, i);
        }

        shifted_cell.set_occupied_full(cell.occupied_mask_full());
        table[insert_pos.table] = shifted_cell;

        entry = cell.entry(cell_type::capacity - 1);
        cell_type::entry_set_shifted(entry, true);
    }

    void insert_and_shift(entry_type entry, entry_pos& insert_pos)
    {
        entry_pos end_pos = insert_pos;
        while (!table[end_pos.table].is_empty(end_pos.cell)) { ++end_pos; }

        while (insert_pos.table != end_pos.table)
        {
            insert_entry_full(entry, insert_pos);
            insert_pos.table++;
            insert_pos.cell = 0;
        };

        insert_entry_not_full(entry, insert_pos, end_pos.cell);
        insert_pos = end_pos;
    }

    static constexpr entry_pos quotient_position(size_t quotient)
    {
        return {static_cast<std::uint32_t>(quotient / cell_type::capacity),
                static_cast<std::uint8_t>(quotient % cell_type::capacity)};
    }

    entry_pos grow_cluster(next_bigger_qf& target, entry_pos cluster_start)
    {
        static constexpr size_t remainder_bits_trgt = remainder_bits - 1;
        using cell_type_trgt  = typename next_bigger_qf::cell_type;
        using entry_pos_trgt  = typename next_bigger_qf::entry_pos;
        using entry_type_trgt = typename next_bigger_qf::entry_type;

        static constexpr entry_type_trgt target_remainder_mask =
            ~((~0ull) << remainder_bits_trgt);
        auto run_start = cluster_start;
        auto run_iter  = cluster_start;

        entry_type      source_entry;
        entry_type_trgt target_entry = 0;
        bool            leading_bit, last_leading_bit;
        entry_pos_trgt  target_pos, last_target_pos = {0, 0};

        auto read = [&]() {
            source_entry = table[run_iter.table].entry(run_iter.cell);
            leading_bit =
                cell_type::entry_remainder(source_entry) >> remainder_bits_trgt;
            cell_type_trgt::entry_set_remainder(
                target_entry, cell_type::entry_remainder(source_entry) &
                                  target_remainder_mask);
        };

        auto advance = [&]() {
            // write
            target.table[target_pos.table].set_entry(target_entry,
                                                     target_pos.cell);

            // backup state
            last_leading_bit = leading_bit;
            last_target_pos  = target_pos;

            target_pos++;
            run_iter++;

            read();
        };

        read();
        target_pos = entry_pos_trgt(run_iter.offset() * 2 +
                                    static_cast<size_t>(leading_bit));

        bool second_pass = false;

        cell_type_trgt::entry_set_occupied(target_entry, false);

        while (!cell_type::entry_is_empty(source_entry))
        {
            // move run_iter start
            cell_type_trgt::entry_set_continuation(target_entry, false);

            if (cell_type::entry_is_cluster_start(source_entry))
            {
                cell_type_trgt::entry_set_shifted(target_entry, false);
            }
            else if (target_pos <= last_target_pos)
            {
                target_pos = last_target_pos;
                target_pos++;
                cell_type_trgt::entry_set_shifted(target_entry, true);
            }
            else
            {
                cell_type_trgt::entry_set_shifted(target_entry, false);
            }

            advance();

            // move continuation of run_iter

            if (cell_type::entry_is_continuation(source_entry))
            {
                cell_type_trgt::entry_set_continuation(target_entry, true);
                cell_type_trgt::entry_set_shifted(target_entry, true);

                while (cell_type::entry_is_continuation(source_entry) &&
                       leading_bit == last_leading_bit)
                {
                    advance();
                }

                if (!second_pass &&
                    cell_type::entry_is_continuation(source_entry))
                {
                    // leading bits differ -> new run_iter start
                    second_pass = true;
                    target_pos =
                        entry_pos_trgt(run_start.offset() * 2 +
                                       static_cast<size_t>(leading_bit));
                    continue;
                }
            }

            if (second_pass)
            {
                target_pos = entry_pos_trgt(run_start.offset() * 2);
                target.table[target_pos.table].set_occupied(true,
                                                            target_pos.cell);
                ++target_pos;
                target.table[target_pos.table].set_occupied(true,
                                                            target_pos.cell);
            }
            else
            {
                target_pos =
                    entry_pos_trgt(run_start.offset() * 2 +
                                   static_cast<size_t>(last_leading_bit));
                target.table[target_pos.table].set_occupied(true,
                                                            target_pos.cell);
            }

            second_pass = false;

            do {
                run_start++;
            } while (run_start.offset() < capacity() &&
                     !table[run_start.table].is_occupied(run_start.cell));

            target_pos = entry_pos_trgt(run_start.offset() * 2 +
                                        static_cast<size_t>(leading_bit));
        }

        return run_start;
    }


    alignas(128) size_t table_capacity;
    alignas(128) std::unique_ptr<cell_type[]> table;
    alignas(128) hash_function_type hf;
    alignas(128) size_t fingerprint_mask;
    alignas(128) size_t quotient_bits;

  public:
    templated_qfilter_seq(size_t capacity = qf::DEFAULT_MIN_CAPACITY,
                          const hash_function_type& hf = hash_function_type())
        : hf(hf), quotient_bits(qf::capacity_to_quotient_bits(capacity))
    {
        table_capacity =
            (ONE << quotient_bits) / cell_type::capacity + qf::shift_buffer;
        table            = std::make_unique<cell_type[]>(table_capacity);
        fingerprint_mask = (ONE << (remainder_bits + quotient_bits)) - 1;
    }

    templated_qfilter_seq(const templated_qfilter_seq& other)
        : table_capacity(other.table_capacity),
          table(new cell_type[table_capacity]), hf(other.hf),
          fingerprint_mask(other.fingerprint_mask),
          quotient_bits(other.quotient_bits)
    {
        for (size_t i = 0; i < table_capacity; i++)
        {
            table[i] = other.table[i];
        }
    }

    templated_qfilter_seq& operator=(const templated_qfilter_seq& other)
    {
        table_capacity   = other.table_capacity;
        table            = std::make_unique<cell_type[]>(other.table_capacity);
        hf               = other.hf;
        fingerprint_mask = other.fingerprint_mask;
        quotient_bits    = other.quotient_bits;

        for (size_t i = 0; i < table_capacity; i++)
        {
            table[i] = other.table[i];
        }

        return *this;
    }

    templated_qfilter_seq(templated_qfilter_seq&& other) = default;
    templated_qfilter_seq& operator=(templated_qfilter_seq&& other) = default;

    ~templated_qfilter_seq() = default;

    base_qf* create_bigger_QF(
        const hash_function_type& hash = hash_function_type()) const override
    {
        if constexpr (remainder_bits > 1)
        {
            return new next_bigger_qf(ONE << (quotient_bits + 1), hash);
        }
        else
        {
            return nullptr;
        }
    }

    void prefetch(const hashed_type& hashed) const override
    {
        const auto [q, r] = this->get_quotient_and_remainder(hashed);
        const auto q_pos  = quotient_position(q);
        __builtin_prefetch(&table[q_pos.first]);
    }

    qf::InsertResult insert(const key_type& key) override
    {
        return insert_hash(hf(key));
    }

    qf::InsertResult insert_hash(const hashed_type& hashed) override
    {
        const auto [q, r] = this->get_quotient_and_remainder(hashed);
        const auto q_pos  = quotient_position(q);

        if (table[q_pos.table].is_empty(q_pos.cell))
        {
            table[q_pos.table].set_remainder(r, q_pos.cell);
            table[q_pos.table].set_occupied(true, q_pos.cell);
            return 1;
        }

        const bool run_already_existing =
            table[q_pos.table].is_occupied(q_pos.cell);
        table[q_pos.table].set_occupied(true, q_pos.cell);

        const auto cluster_start = find_cluster_start(q_pos);
        const auto run           = find_run_start(cluster_start, q_pos);

        if (run.table == table_capacity - 1) { return qf::FAILED_OPERATION; }

        auto       insert_pos = run;
        entry_type new_entry  = 0;
        cell_type::entry_set_remainder(new_entry, r);

        if (run_already_existing)
        {
            // find position in existing run
            do {
                if (table[insert_pos.table].remainder(insert_pos.cell) == r)
                    // element already present
                    return entry_pos::distance(cluster_start, insert_pos);
                else if (table[insert_pos.table].remainder(insert_pos.cell) > r)
                    break;

                ++insert_pos;
            } while (table[insert_pos.table].is_continuation(insert_pos.cell));

            // element is not present so we have to shift
            if (!has_space_for_shift(insert_pos))
            {
                return qf::FAILED_OPERATION;
            }

            if (insert_pos == run)
            {
                // old run start will be shifted
                table[run.table].set_continuation(true, run.cell);
            }
            else
            {
                // new entry is not run start
                cell_type::entry_set_continuation(new_entry, true);
            }
        }
        else if (!has_space_for_shift(insert_pos))
        {
            table[q_pos.table].set_occupied(false, q_pos.cell);
            return qf::FAILED_OPERATION;
        }

        if (insert_pos != q_pos)
        {
            cell_type::entry_set_shifted(new_entry, true);
        }

        insert_and_shift(new_entry, insert_pos);
        return entry_pos::distance(cluster_start, insert_pos);
    }

    bool quick_insert(const key_type& key) override
    {
        return quick_insert_hash(hf(key));
    }

    bool quick_insert_hash(const hashed_type& hashed) override
    {
        const auto [q, r] = this->get_quotient_and_remainder(hashed);
        const auto q_pos  = quotient_position(q);

        if (table[q_pos.table].is_empty(q_pos.cell))
        {
            table[q_pos.table].set_remainder(r, q_pos.cell);
            table[q_pos.table].set_occupied(true, q_pos.cell);
            return true;
        }

        return contains_hash(hashed);
    }

    qf::ContainsResult contains(const key_type& key) const override
    {
        return contains_hash(hf(key));
    }

    qf::ContainsResult contains_hash(const hashed_type& hashed) const override
    {
        const auto [q, r]        = this->get_quotient_and_remainder(hashed);
        const auto       q_pos   = quotient_position(q);
        const entry_type q_entry = table[q_pos.table].entry(q_pos.cell);

        if (cell_type::entry_is_empty(q_entry))
            return ContainsResultEnum::CANONICAL_SLOT_EMPTY;
        if (!cell_type::entry_is_occupied(q_entry))
            return ContainsResultEnum::ENTRY_NOT_PRESENT;

        auto iter = find_run_start(find_cluster_start(q_pos), q_pos);

        do {
            if (table[iter.table].remainder(iter.cell) == r)
                return ContainsResultEnum::ENTRY_PRESENT;

            ++iter;
        } while (table[iter.table].is_continuation(iter.cell));

        return ContainsResultEnum::ENTRY_NOT_PRESENT;
    }

    std::pair<size_t, entry_type>
    get_quotient_and_remainder(const hashed_type& hashed) const
    {
        return qf::get_quotient_and_remainder<entry_type, remainder_bits>(
            hashed, quotient_bits);
    }

    size_t capacity() const override { return ONE << quotient_bits; }

    size_t memory_usage_bytes() const override
    {
        return table_capacity * sizeof(cell_type);
    }

    size_t unused_memory_bits() const override
    {
        return table_capacity *
               (sizeof(cell_type) * 8 -
                cell_type::capacity * (remainder_bits + qf::status_bits));
    }

    double fill_level() const override
    {
        size_t used_slots = 0;

        for (entry_pos pos{0, 0}; pos.table < table_capacity; pos++)
        {
            if (!table[pos.table].is_empty(pos.cell)) used_slots++;
        }

        return static_cast<double>(used_slots) / table_capacity /
               cell_type::capacity;
    }

    virtual void grow(base_qf& target_base) override
    {
        if constexpr (remainder_bits <= 1) { return; }
        else
        {
            auto& target = dynamic_cast<next_bigger_qf&>(target_base);

            entry_pos cluster_start{0, 0};

            while (table[cluster_start.table].is_empty(cluster_start.cell))
            {
                cluster_start++;
            }

            while (cluster_start.offset() < ONE << quotient_bits)
            {
                cluster_start = this->grow_cluster(target, cluster_start);
            }
        }
    }

    void print()
    {
        std::cout << "occu | cont | shif || qoutient | remainder\n";

        for (size_t table_pos = 0; table_pos < table_capacity; table_pos++)
        {
            for (size_t cell_pos = 0; cell_pos < cell_type::capacity;
                 cell_pos++)
            {
                print_entry({table_pos, cell_pos});
                std::cout << "   ";
            }
            std::cout << "\n";
        }

        std::cout << std::endl;
    }

    void print_entry(entry_pos pos)
    {
        print_entry(table[pos.table].entry(pos.cell), pos.offset());
    }

    void print_entry(entry_type entry, size_t index)
    {
        std::cout << cell_type::entry_is_occupied(entry) << "|"
                  << cell_type::entry_is_continuation(entry) << "|"
                  << cell_type::entry_is_shifted(entry) << " || "
                  << qf::bitstring(index, quotient_bits) << "|"
                  << qf::bitstring(cell_type::entry_remainder(entry),
                                   remainder_bits);
    }

    void print_entry(entry_type entry)
    {
        std::cout << cell_type::entry_is_occupied(entry) << "|"
                  << cell_type::entry_is_continuation(entry) << "|"
                  << cell_type::entry_is_shifted(entry) << " ||     "
                  << qf::bitstring(cell_type::entry_remainder(entry),
                                   remainder_bits);
    }

    bool check_consistency() override
    {
        bool       consistent = true;
        bool       sorted     = true;
        entry_type last       = 0;
        cell_type  current;

        for (entry_pos pos{0, 0}; pos.table < table_capacity; pos++)
        {
            current = table[pos.table];

            if (current.is_continuation(pos.cell) &&
                !current.is_shifted(pos.cell))
            {
                std::cerr << "error1: [" << pos << "] cont =/=> shift"
                          << std::endl;
                consistent = false;
            }
            else if (current.is_continuation(pos.cell) &&
                     !cell_type::entry_is_shifted(last) &&
                     !cell_type::entry_is_cluster_start(last))
            {
                std::cerr << "error2: [" << pos << "] cont =/=> shift"
                          << std::endl;
                consistent = false;
            }
            else if (current.is_continuation(pos.cell))
            {
                if (current.remainder(pos.cell) <
                    cell_type::entry_remainder(last))
                {
                    sorted     = false;
                    consistent = false;
                }
            }

            last = current.entry(pos.cell);
        }

        if (!sorted) std::cerr << "error:  remainder not sorted" << std::endl;

        return consistent;
    }
};

} // namespace qf
