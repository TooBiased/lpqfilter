#pragma once
/*******************************************************************************
 * implementation/base_filter/nongrouped/nongrouped_qfilter_seq.hpp
 *
 * sequential quotient filter with non-grouped slots (one slot per data slot)
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/



#include <iostream>
#include <memory>
#include <unordered_map>


#include "utils/default_hash.hpp"
namespace htm = utils_tm::hash_tm;
#include "implementation/definitions.hpp"
#include "implementation/utilities.hpp"
#include "nongrouped_qfilter_cell.hpp"


namespace qf
{

class qf_histograms
{
  public:
    using histogram = std::unordered_map<size_t, size_t>;

  private:
    histogram runs;
    histogram cluster;
    histogram super_cluster;
    size_t    max_super_cluster = 0;

    size_t get_count(histogram& hist, size_t val)
    {
        auto iter = hist.find(val);
        if (iter == hist.end()) { return 0; }
        else
        {
            return iter->second;
        }
    }

    void add_entry(histogram& hist, size_t& val)
    {
        if (val == 0) return;

        auto iter = hist.find(val);
        if (iter == hist.end()) { hist[val] = 1; }
        else
        {
            iter->second++;
        }

        val = 0;
    }

  public:
    size_t max_super_cluster_length() { return max_super_cluster; }

    size_t get_run_count(size_t length) { return get_count(runs, length); }

    size_t get_cluster_count(size_t length)
    {
        return get_count(cluster, length);
    }

    size_t get_super_cluster_count(size_t length)
    {
        return get_count(super_cluster, length);
    }

    void add_run(size_t& val) { add_entry(runs, val); }

    void add_cluster(size_t& val) { add_entry(cluster, val); }

    void add_super_cluster(size_t& val)
    {
        max_super_cluster = val > max_super_cluster ? val : max_super_cluster;
        add_entry(super_cluster, val);
    }
};

template <class Key, class Hash = htm::default_hash>
class nongrouped_qfilter_seq
{
  public:
    using key_type = Key;
    using hasher   = Hash;
    using HashVal  = decltype(std::declval<Hash>()(std::declval<Key>()));

  private:
    using entry_type = std::uint64_t;
    using cell       = qf::nongrouped_cell_non_atomic<entry_type>;


    size_t find_cluster_start(size_t q) const
    {
        while (table[q].is_shifted()) { q--; }

        return q;
    }

    size_t find_run_start(size_t cluster_start, size_t q) const
    {
        size_t bucket = cluster_start;
        size_t run    = bucket;

        while (bucket != q)
        {
            do {
                if (run == table_capacity - 1)
                {
                    return run; // invalid value
                }

                run++;
            } while (table[run].is_continuation());

            do {
                bucket++;
            } while (!table[bucket].is_occupied());
        }

        return run;
    }

    // insert a remainder in a slot that is not its canonical slot
    void insert_and_shift(cell entry, size_t& pos)
    {
        cell prev;
        cell curr = entry;
        bool empty;

        do {
            prev  = table[pos];
            empty = prev.is_empty();

            if (!empty)
            {
                prev.set_shifted(true);
                if (prev.is_occupied())
                {
                    prev.set_occupied(false);
                    curr.set_occupied(true);
                }
            }

            table[pos] = curr;
            curr       = prev;
            pos++;
        } while (!empty);
        pos--;
    }

    // assumes that pos is not empty
    bool has_space_for_shift(size_t pos) const
    {
        while (!table[pos].is_empty()) { pos++; }

        return pos != table_capacity - 1;
    }

    static size_t
    grow_cluster(nongrouped_qfilter_seq& source, nongrouped_qfilter_seq& target,
                 size_t cluster_start)
    {
        const entry_type target_remainder_mask =
            ~((~static_cast<entry_type>(0)) << target.remainder_bits);
        auto   src_table_size = source.table_capacity;
        size_t run_start      = cluster_start;
        size_t run_iter       = cluster_start;

        cell   source_entry, target_entry;
        bool   leading_bit, last_leading_bit;
        size_t target_index, last_target_index = 0;

        auto read = [&]() {
            source_entry = source.table[run_iter];
            leading_bit  = source_entry.remainder() >> target.remainder_bits;
            target_entry.set_remainder(source_entry.remainder() &
                                       target_remainder_mask);
        };

        auto advance = [&]() {
            // write
            target.table[target_index] = target_entry;

            // backup state
            last_leading_bit  = leading_bit;
            last_target_index = target_index;

            target_index++;
            run_iter++;

            read();
        };

        read();
        target_index = run_iter * 2 + static_cast<size_t>(leading_bit);

        bool second_pass = false;

        target_entry.set_occupied(false);

        while (!source_entry.is_empty())
        {
            // move run_iter start
            target_entry.set_continuation(false);

            if (source_entry.is_cluster_start())
            {
                target_entry.set_shifted(false);
            }
            else if (target_index < last_target_index + 1)
            {
                target_index = last_target_index + 1;
                target_entry.set_shifted(true);
            }
            else
            {
                target_entry.set_shifted(false);
            }

            advance();

            // move continuation of run_iter

            if (source_entry.is_continuation())
            {
                target_entry.set_continuation(true);
                target_entry.set_shifted(true);

                while (source_entry.is_continuation() &&
                       leading_bit == last_leading_bit)
                {
                    advance();
                }

                if (!second_pass && source_entry.is_continuation())
                {
                    // leading bits differ -> new run_iter start
                    second_pass = true;
                    target_index =
                        run_start * 2 + static_cast<size_t>(leading_bit);
                    continue;
                }
            }

            if (second_pass)
            {
                target.table[run_start * 2].set_occupied(true);
                target.table[run_start * 2 + 1].set_occupied(true);
            }
            else
            {
                target
                    .table[run_start * 2 +
                           static_cast<size_t>(last_leading_bit)]
                    .set_occupied(true);
            }

            second_pass = false;

            do {
                run_start++;
            } while (run_start < src_table_size &&
                     !source.table[run_start].is_occupied());

            target_index = run_start * 2 + static_cast<size_t>(leading_bit);
        }

        return run_start;
    }

    alignas(128) size_t table_capacity;
    alignas(128) std::unique_ptr<cell[]> table;
    alignas(128) hasher hf;
    alignas(128) size_t remainder_bits;
    alignas(128) size_t quotient_bits;

    static constexpr entry_type one = static_cast<entry_type>(1);

  public:
    nongrouped_qfilter_seq(size_t capacity       = qf::DEFAULT_MIN_CAPACITY,
                           size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS,
                           const hasher& hf      = hasher())
        : hf(hf), remainder_bits(remainder_bits),
          quotient_bits(qf::capacity_to_quotient_bits(capacity))
    {
        table_capacity = (one << quotient_bits) + qf::shift_buffer;
        table          = std::make_unique<cell[]>(table_capacity);
    }

    nongrouped_qfilter_seq(const nongrouped_qfilter_seq& other)
        : table_capacity(other.table_capacity), table(new cell[table_capacity]),
          hf(other.hf), remainder_bits(other.remainder_bits),
          quotient_bits(other.quotient_bits)
    {
        for (size_t i = 0; i < table_capacity; i++)
        {
            table[i].entry = other.table[i].entry;
        }
    }

    nongrouped_qfilter_seq& operator=(const nongrouped_qfilter_seq& other)
    {
        table_capacity = other.table_capacity;
        table          = std::make_unique<cell[]>(other.table_capacity);
        hf             = other.hf;
        remainder_bits = other.remainder_bits;
        quotient_bits  = other.quotient_bits;

        for (size_t i = 0; i < table_capacity; i++)
        {
            table[i].entry = other.table[i].entry;
        }

        return *this;
    }

    nongrouped_qfilter_seq(nongrouped_qfilter_seq&& other) = default;
    nongrouped_qfilter_seq& operator=(nongrouped_qfilter_seq&& other) = default;

    void prefetch(const HashVal& hashed) const
    {
        auto [q, r] = get_quotient_and_remainder(hashed);
        __builtin_prefetch(&table[q]);
    }

    qf::InsertResult insert(const key_type& key)
    {
        return insert_hash(hf(key));
    }

    qf::InsertResult insert_hash(const HashVal& hashed)
    {
        auto [q, r] = get_quotient_and_remainder(hashed);

        if (table[q].is_empty())
        {
            table[q].set_remainder(r);
            table[q].set_occupied(true);
            return 1;
        }

        bool run_already_existing = table[q].is_occupied();
        table[q].set_occupied(true);

        const size_t cluster_start = find_cluster_start(q);
        const size_t run           = find_run_start(cluster_start, q);

        if (run == table_capacity - 1) { return qf::FAILED_OPERATION; }

        size_t insert_pos = run;

        cell new_entry;
        new_entry.set_remainder(r);

        if (run_already_existing)
        {
            // find position in existing run_iter
            do {
                if (table[insert_pos].remainder() == r)
                    // element already present
                    return insert_pos - cluster_start;
                else if (table[insert_pos].remainder() > r)
                    break;

                insert_pos++;
            } while (table[insert_pos].is_continuation());

            // element is not present so we have to shift
            if (!has_space_for_shift(insert_pos))
            {
                return qf::FAILED_OPERATION;
            }

            if (insert_pos == run)
            {
                // old run_iter start will be shifted
                table[run].set_continuation(true);
            }
            else
            {
                // new entry is not run_iter start
                new_entry.set_continuation(true);
            }
        }
        else if (!has_space_for_shift(insert_pos))
        {
            table[q].set_occupied(false);
            return qf::FAILED_OPERATION;
        }

        if (insert_pos != q) { new_entry.set_shifted(true); }

        insert_and_shift(new_entry, insert_pos);
        return insert_pos - cluster_start;
    }

    bool quick_insert(const key_type& key)
    {
        return quick_insert_hash(hf(key));
    }

    bool quick_insert_hash(const HashVal& hashed)
    {
        auto [q, r] = get_quotient_and_remainder(hashed);

        if (table[q].is_empty())
        {
            table[q].set_remainder(r);
            table[q].set_occupied(true);
            return true;
        }

        return contains_hash(hashed);
    }

    qf::ContainsResult contains(const key_type& key) const
    {
        return contains_hash(hf(key));
    }

    qf::ContainsResult contains_hash(const HashVal& hashed) const
    {
        auto [q, r] = get_quotient_and_remainder(hashed);

        if (table[q].is_empty())
            return ContainsResultEnum::CANONICAL_SLOT_EMPTY;
        if (!table[q].is_occupied())
            return ContainsResultEnum::ENTRY_NOT_PRESENT;

        size_t iter = find_run_start(find_cluster_start(q), q);

        do {
            if (table[iter].remainder() == r)
                return ContainsResultEnum::ENTRY_PRESENT;

            iter++;
        } while (table[iter].is_continuation());

        return ContainsResultEnum::ENTRY_NOT_PRESENT;
    }

    std::pair<size_t, entry_type>
    get_quotient_and_remainder(const HashVal& hashed) const
    {
        return qf::get_quotient_and_remainder<entry_type>(hashed, quotient_bits,
                                                          remainder_bits);
    }

    size_t capacity() const { return ONE << quotient_bits; }


    size_t memory_usage_bytes() const { return table_capacity * sizeof(cell); }

    size_t unused_memory_bits() const
    {
        return table_capacity *
               (sizeof(cell) * 8 - remainder_bits - qf::status_bits);
    }

    double fill_level() const
    {
        size_t used_slots = 0;

        for (size_t pos = 0; pos < table_capacity; pos++)
        {
            if (!table[pos].is_empty()) used_slots++;
        }

        return static_cast<double>(used_slots) / table_capacity;
    }

    void print()
    {
        std::cout << "occu | cont | shif || qoutient | remainder\n";

        for (size_t i = 0; i < table_capacity; i++)
        {
            print_entry(i);
            std::cout << "\n";
        }

        std::cout << std::endl;
    }

    void print_entry(size_t index) { print_entry(table[index], index); }

    void print_entry(cell entry, size_t index)
    {
        std::cout << entry.is_occupied() << "|" << entry.is_continuation()
                  << "|" << entry.is_shifted() << " || "
                  << qf::bitstring(index, quotient_bits) << "|"
                  << qf::bitstring(entry.remainder(), remainder_bits);
    }

    void grow(nongrouped_qfilter_seq& target)
    {
        size_t cluster_start = 0;

        while (table[cluster_start].is_empty()) { cluster_start++; }

        while (cluster_start < table_capacity)
        {
            cluster_start = grow_cluster(*this, target, cluster_start);
        }
    }

    bool check_consistency()
    {
        bool consistent = true;
        bool sorted     = true;
        cell last;
        cell current;

        for (size_t pos = 0; pos < table_capacity; pos++)
        {
            current = table[pos];

            if (current.is_continuation() && !current.is_shifted())
            {
                std::cerr << "error1: [" << pos << "] cont =/=> shift"
                          << std::endl;
                consistent = false;
            }
            else if (current.is_continuation() && !last.is_shifted() &&
                     !last.is_cluster_start())
            {
                std::cerr << "error2: [" << pos << "] cont =/=> shift"
                          << std::endl;
                consistent = false;
            }
            else if (current.is_continuation())
            {
                if (current.remainder() < last.remainder())
                {
                    sorted     = false;
                    consistent = false;
                }
            }

            last = current;
        }

        if (!sorted) std::cerr << "error:  remainder not sorted" << std::endl;

        return consistent;
    }

    qf_histograms get_histograms()
    {
        qf_histograms histograms;

        size_t run           = 0;
        size_t cluster       = 0;
        size_t super_cluster = 0;
        cell   current;

        for (size_t pos = 0; pos < table_capacity; pos++)
        {
            current = table[pos];

            if (current.is_empty())
            {
                histograms.add_super_cluster(super_cluster);
                histograms.add_cluster(cluster);
                histograms.add_run(run);

                continue;
            }
            super_cluster++;

            if (current.is_cluster_start())
            {
                histograms.add_cluster(cluster);
                histograms.add_run(run);

                cluster = 1;
                run     = 1;
                continue;
            }
            cluster++;

            if (!current.is_continuation())
            {
                histograms.add_run(run);
                run = 1;
                continue;
            }
            run++;
        }

        return histograms;
    }
};

} // namespace qf
