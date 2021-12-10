#pragma once
/*******************************************************************************
 * implementation/base_filter/standard/standard_qfilter_conc.hpp
 *
 * concurrent quotient filter with grouped slots
 * (multiple slots per atomic non-templated)
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/



#include <iostream>
#include <memory>
#include <utility>

#include "utils/default_hash.hpp"
namespace htm = utils_tm::hash_tm;
#include "implementation/utilities.hpp"
#include "standard_qfilter_cell.hpp"


namespace qf
{

template <class Key, class Hash> class standard_qfilter_conc;

template <class Key, class Hash = htm::default_hash>
class standard_qfilter_conc_base
{
  protected:
    alignas(128) std::atomic<qf::TableStatus> table_status;

    standard_qfilter_conc_base(const standard_qfilter_conc_base& other)
        : table_status(other.table_status.load())
    {
    }

    standard_qfilter_conc_base&
    operator=(const standard_qfilter_conc_base& other)
    {
        table_status = other.table_status.load();
        return *this;
    }

    standard_qfilter_conc_base(standard_qfilter_conc_base&&) = default;
    standard_qfilter_conc_base&
    operator=(standard_qfilter_conc_base&&) = default;

    standard_qfilter_conc_base(
        qf::TableStatus table_status = qf::TableStatus::active)
        : table_status(table_status)
    {
    }

  public:
    using key_type           = Key;
    using hash_function_type = Hash;
    using hashed_type = decltype(std::declval<Hash>()(std::declval<Key>()));

    virtual ~standard_qfilter_conc_base() = default;

    virtual void               prefetch(const hashed_type& hashed) const    = 0;
    virtual qf::InsertResult   insert(const key_type& key)                  = 0;
    virtual qf::InsertResult   insert_hash(const hashed_type& hashed)       = 0;
    virtual bool               quick_insert(const key_type& key)            = 0;
    virtual bool               quick_insert_hash(const hashed_type& hashed) = 0;
    virtual qf::ContainsResult contains(const key_type& key)                = 0;
    virtual qf::ContainsResult contains_hash(const hashed_type& hashed)     = 0;
    virtual qf::ContainsResult unsafe_contains(const key_type& key) const   = 0;
    virtual qf::ContainsResult
                 unsafe_contains_hash(const hashed_type& hashed) const = 0;
    virtual bool check_consistency() const                             = 0;
    virtual standard_qfilter_conc_base* create_bigger_QF(
        const hash_function_type& hf = hash_function_type()) const = 0;
    virtual void   grow(standard_qfilter_conc_base<Key, Hash>& target_base,
                        size_t                                 block_size)                         = 0;
    virtual size_t capacity() const                                = 0;

    void set_table_status(qf::TableStatus status) { table_status = status; }
};


template <class Key, class Hash = htm::default_hash>
class standard_qfilter_conc : public standard_qfilter_conc_base<Key, Hash>
{
  public:
    static constexpr bool is_growing_compatible = true;
    static constexpr bool is_templated          = false;
    static constexpr bool is_sequential         = false;
    static constexpr bool is_growing            = false;
    static constexpr bool is_dynamic            = false;
    static constexpr bool uses_handle           = false;

    using this_type          = standard_qfilter_conc<Key, Hash>;
    using key_type           = Key;
    using hash_function_type = Hash;
    using hashed_type     = decltype(std::declval<Hash>()(std::declval<Key>()));
    using next_smaller_qf = standard_qfilter_conc<Key, Hash>;
    using next_bigger_qf  = standard_qfilter_conc<Key, Hash>;
    using base_qf         = standard_qfilter_conc_base<Key, Hash>;

    friend base_qf;
    // friend next_smaller_qf;

  private:
    using cell_base_type = qf::grouped_cell_base_type;
    using cell_type      = qf::standard_cell_non_atomic<cell_base_type>;
    using cell_atomic    = typename cell_type::atomic;

    using entry_type = typename cell_type::entry_type;
    using cell_stuff = typename cell_type::cstuff;
    using entry_pos  = std::pair<size_t, size_t>;

    using read_lock_guard  = qf::compact_read_lock_guard<cell_atomic>;
    using write_lock_guard = qf::compact_write_lock_guard<cell_atomic>;

    const cell_stuff cc;
    const size_t     quotient_bits;
    const size_t     table_capacity;
    const size_t     fingerprint_mask;
    alignas(128) std::atomic_size_t growing_block_counter  = 0;
    alignas(128) std::atomic_size_t finished_block_counter = 0;
    alignas(128) std::unique_ptr<cell_atomic[]> table;
    alignas(128) const hash_function_type hf;

  public:
    standard_qfilter_conc(size_t capacity       = qf::DEFAULT_MIN_CAPACITY,
                          size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS,
                          const hash_function_type& hf = hash_function_type());

    standard_qfilter_conc(const standard_qfilter_conc& other);
    standard_qfilter_conc& operator=(const standard_qfilter_conc& other);

    standard_qfilter_conc(standard_qfilter_conc&& other) = default;
    standard_qfilter_conc& operator=(standard_qfilter_conc&& other) = default;

    ~standard_qfilter_conc() = default;

    this_type* create_bigger_QF(
        const hash_function_type& hash = hash_function_type()) const override;

    void             prefetch(const hashed_type& hashed) const override;
    qf::InsertResult insert(const key_type& key) override;
    qf::InsertResult insert_hash(const hashed_type& hashed) override;
    bool             quick_insert(const key_type& key) override;
    bool             quick_insert_hash(const hashed_type& hashed) override;

    qf::ContainsResult contains(const key_type& key) override;
    qf::ContainsResult contains_hash(const hashed_type& hashed) override;
    qf::ContainsResult unsafe_contains(const key_type& key) const override;
    qf::ContainsResult
    unsafe_contains_hash(const hashed_type& hashed) const override;

    std::pair<size_t, entry_type>
    get_quotient_and_remainder(const hashed_type& hashed) const;

    size_t capacity() const;

    virtual void grow(base_qf& target_base, size_t block_size) override;

  private:
    entry_pos quotient_position(size_t quotient) const;
    void      increment(entry_pos& pos) const;
    void      decrement(entry_pos& pos) const;
    size_t    offset(const entry_pos& pos) const;
    size_t    distance(const entry_pos& start, const entry_pos& end) const;

    ContainsResultEnum
    try_contains_lock_elision(const cell_type cell, const entry_pos pos,
                              const entry_type remainder) const;
    qf::InsertResult
    try_insert_lock_elision(cell_type cell, const entry_pos pos,
                            const entry_type remainder) const;

    void increment_and_update(entry_pos& pos, cell_type& cell,
                              entry_type& entry) const;

    void decrement_and_update(entry_pos& pos, cell_type& cell,
                              entry_type& entry) const;

    void force_update(const entry_pos& pos, cell_type& cell,
                      entry_type& entry) const;

    bool
    try_empty_slot(cell_type local_cell, entry_pos pos, entry_type remainder);

    bool unsafe_try_empty_slot(entry_pos pos, entry_type remainder);

    // assumes that pos is not empty
    bool has_space_and_lock(entry_pos pos, write_lock_guard& write_lock);

    // assumes that pos is not empty
    void get_growing_write_lock(entry_pos pos, write_lock_guard& write_lock);

    entry_pos
    find_and_lock_cluster_start(entry_pos pos, read_lock_guard& read_lock);

    entry_pos find_cluster_start(entry_pos pos) const;

    entry_pos find_run_start(entry_pos cluster_start, entry_pos pos,
                             bool set_occupied = false) const;

    // does not use or check for locks
    entry_pos
    unsafe_find_run_start(entry_pos cluster_start, entry_pos pos) const;

    bool acquire_read_lock(entry_pos pos, read_lock_guard& read_lock);

    // insert
    void
    insert_and_shift(entry_type entry, entry_pos& insert_pos,
                     write_lock_guard& write_lock, read_lock_guard& read_lock);
    // inserts entry at pos and shifts the following entries cell at pos
    // must have enough space for the shift i.e. an empty entry (write lock)
    void insert_entry_not_full(entry_type entry, entry_pos insert_pos,
                               size_t write_lock_pos);

    // inserts entry at pos and shifts the following entries the last entry
    // which is lost due to the shift is written back to the entry reference
    void insert_entry_full(entry_type& entry, entry_pos insert_pos);

    entry_type wait_on_insert_write_lock(entry_pos pos) const;

    void grow(next_bigger_qf& target, size_t block_start, size_t block_size);
    entry_pos grow_cluster(next_bigger_qf&        target,
                           const read_lock_guard& read_lock) const;

  public:
    void print() const;
    void print(entry_pos start, entry_pos end) const;
    void print(size_t start, size_t count) const;
    void print_cell(size_t table_pos) const;
    void print_cell(cell_type cell, size_t table_pos) const;
    void print_entry(entry_pos pos) const;
    void print_entry(entry_type entry, entry_pos pos) const;
    void print_entry(entry_type entry, size_t index) const;
    void print_entry(entry_type entry) const;
    bool check_consistency() const override;
};



// *****************************************************************************
// *** CONSTRUCTORS ************************************************************
// *****************************************************************************

template <class K, class H>
standard_qfilter_conc<K, H>::standard_qfilter_conc(
    size_t capacity, size_t remainder_bits, const hash_function_type& _hf)
    : base_qf(), cc(remainder_bits),
      quotient_bits(qf::capacity_to_quotient_bits(capacity)),
      table_capacity((ONE << quotient_bits) / cc.capacity + qf::shift_buffer),
      fingerprint_mask(
          (ONE << static_cast<size_t>(remainder_bits + this->quotient_bits)) -
          1),
      hf(_hf)
{
    table = std::make_unique<cell_atomic[]>(table_capacity);
}

template <class K, class H>
standard_qfilter_conc<K, H>::standard_qfilter_conc(
    const standard_qfilter_conc& other)
    : base_qf(other), cc(other.cc), quotient_bits(other.quotient_bits),
      table_capacity(other.table_capacity),
      fingerprint_mask(other.fingerprint_mask),
      growing_block_counter(0),  // other.growing_block_counter),
      finished_block_counter(0), // other.finished_block_counter),
      table(new cell_atomic[table_capacity]), hf(other.hf)
{
    for (size_t i = 0; i < table_capacity; i++) { table[i] = other.table[i]; }
}

template <class K, class H>
standard_qfilter_conc<K, H>&
standard_qfilter_conc<K, H>::operator=(const standard_qfilter_conc& other)
{
    if (this == &other) return *this;

    this->~this_type();
    new (this) this_type(other);
    return *this;
}

template <class K, class H>
standard_qfilter_conc<K, H>* standard_qfilter_conc<K, H>::create_bigger_QF(
    const hash_function_type& hash) const
{
    if (cc.remainder_bits > 1)
    {
        return new next_bigger_qf(ONE << (this->quotient_bits + 1),
                                  cc.remainder_bits - 1, hash);
    }
    else
    {
        return nullptr;
    }
}



// *****************************************************************************
// *** MAIN FUNCTIONALITY ******************************************************
// *****************************************************************************

// *** INSERT ******************************************************************

template <class K, class H>
void standard_qfilter_conc<K, H>::prefetch(const hashed_type& hashed) const
{
    const auto [q, r]     = this->get_quotient_and_remainder(hashed);
    const entry_pos q_pos = quotient_position(q);
    __builtin_prefetch(&table[q_pos.first]);
}


template <class K, class H>
qf::InsertResult standard_qfilter_conc<K, H>::insert(const key_type& key)
{
    return insert_hash(hf(key));
}

template <class K, class H>
qf::InsertResult
standard_qfilter_conc<K, H>::insert_hash(const hashed_type& hashed)
{
    if (this->table_status != qf::TableStatus::active)
        return qf::FAILED_OPERATION;

    const auto [q, r]      = this->get_quotient_and_remainder(hashed);
    const entry_pos q_pos  = quotient_position(q);
    const cell_type q_cell = table[q_pos.first];

    if (try_empty_slot(q_cell, q_pos, r)) return 1;

    if (cc.capacity > 1)
    {
        auto lock_elision_result = try_insert_lock_elision(q_cell, q_pos, r);
        if (lock_elision_result.successful()) return lock_elision_result;
    }

    write_lock_guard write_lock(cc);
    if (!has_space_and_lock(q_pos, write_lock)) return qf::FAILED_OPERATION;

    // works even if q_pos is read_locked since it was occupied before locking
    // and read_lock mask set the occupied bit
    const bool run_already_existing =
        table[q_pos.first].is_occupied(cc, q_pos.second);

    read_lock_guard read_lock(cc);
    const auto cluster_start = find_and_lock_cluster_start(q_pos, read_lock);
    const auto run =
        find_run_start(cluster_start, q_pos, !run_already_existing);
    auto insert_pos = run;

    entry_type new_entry = 0;
    cell_type::entry_set_remainder(cc, new_entry, r);

    if (insert_pos == write_lock.pos)
    {
        cell_type::entry_set_shifted(cc, new_entry, true);
        table[write_lock.pos.first].set_entry(cc, new_entry,
                                              write_lock.pos.second);
        write_lock.clear();
        return distance(cluster_start, insert_pos);
    }
    else if (run.first == table_capacity - 1)
    {
        if (!run_already_existing)
            table[q_pos.first].set_occupied(cc, run_already_existing,
                                            q_pos.second);

        return qf::FAILED_OPERATION;
    }

    cell_type  current_cell  = table[insert_pos.first];
    entry_type current_entry = current_cell.entry(cc, insert_pos.second);

    if (run_already_existing)
    {
        // find position in existing run
        do {
            if (cell_type::entry_remainder(cc, current_entry) == r)
                // element already present
                return distance(cluster_start, insert_pos);
            else if (cell_type::entry_remainder(cc, current_entry) > r)
                break;

            increment_and_update(insert_pos, current_cell, current_entry);

            while (cell_type::entry_is_read_locked(cc, current_entry))
            {
                force_update(insert_pos, current_cell, current_entry);
            }
        } while (!cell_type::entry_is_locked(cc, current_entry) &&
                 cell_type::entry_is_continuation(cc, current_entry));

        if (insert_pos == run)
        {
            // note: possibly writing to read locked entry but only modifying
            //       the remainder
            table[run.first].set_remainder(cc, r, run.second);
            cell_type::entry_set_remainder(
                cc, new_entry, cell_type::entry_remainder(cc, current_entry));

            increment_and_update(insert_pos, current_cell, current_entry);

            while (cell_type::entry_is_read_locked(cc, current_entry))
            {
                force_update(insert_pos, current_cell, current_entry);
            }
        }

        cell_type::entry_set_continuation(cc, new_entry, true);
    }

    // TODO: remove this?
    if (cell_type::entry_is_read_locked(cc, current_entry) ||
        cell_type::entry_is_cluster_start(cc, current_entry))
        acquire_read_lock(insert_pos, read_lock);

    cell_type::entry_set_shifted(cc, new_entry, insert_pos != q_pos);

    insert_and_shift(new_entry, insert_pos, write_lock, read_lock);
    return distance(cluster_start, insert_pos);
}


template <class K, class H>
bool standard_qfilter_conc<K, H>::quick_insert(const key_type& key)
{
    return quick_insert_hash(hf(key));
}

// try insert without shifting (does not use or check for locks)
template <class K, class H>
bool standard_qfilter_conc<K, H>::quick_insert_hash(const hashed_type& hashed)
{
    const auto [q, r]     = this->get_quotient_and_remainder(hashed);
    const entry_pos q_pos = quotient_position(q);

    if (unsafe_try_empty_slot(q_pos, r)) return 1;

    return unsafe_contains_hash(hashed);
}


// *** CONTAINS ****************************************************************

template <class K, class H>
qf::ContainsResult standard_qfilter_conc<K, H>::contains(const key_type& key)
{
    return contains_hash(hf(key));
}

template <class K, class H>
qf::ContainsResult
standard_qfilter_conc<K, H>::contains_hash(const hashed_type& hashed)
{
    if (this->table_status != qf::TableStatus::active)
        return qf::FAILED_OPERATION;

    const auto [q, r]        = this->get_quotient_and_remainder(hashed);
    const entry_pos  q_pos   = quotient_position(q);
    const cell_type  q_cell  = table[q_pos.first];
    const entry_type q_entry = q_cell.entry(cc, q_pos.second);

    if (cell_type::entry_is_empty(cc, q_entry))
        return ContainsResultEnum::CANONICAL_SLOT_EMPTY;
    if (!cell_type::entry_is_occupied(cc, q_entry))
        return ContainsResultEnum::ENTRY_NOT_PRESENT;

    if (cc.capacity > 1)
    {
        auto lock_elision_result = try_contains_lock_elision(q_cell, q_pos, r);
        if (lock_elision_result != ContainsResultEnum::FAILED_OPERATION)
            return lock_elision_result;
    }

    read_lock_guard read_lock(cc);
    auto            iter =
        find_run_start(find_and_lock_cluster_start(q_pos, read_lock), q_pos);

    cell_type  current_cell  = table[iter.first];
    entry_type current_entry = current_cell.entry(cc, iter.second);

    do {
        if (cell_type::entry_remainder(cc, current_entry) == r)
            return ContainsResultEnum::ENTRY_PRESENT;

        // early abort due to sort is acutally slower in + and - contains
        // if (cell_type::entry_remainder(current_entry) > r)
        //	return ContainsResultEnum::ENTRY_NOT_PRESENT;

        increment_and_update(iter, current_cell, current_entry);
    } while (!cell_type::entry_is_locked(cc, current_entry) &&
             cell_type::entry_is_continuation(cc, current_entry));

    return ContainsResultEnum::ENTRY_NOT_PRESENT;
}

template <class K, class H>
qf::ContainsResult
standard_qfilter_conc<K, H>::unsafe_contains(const key_type& key) const
{
    return unsafe_contains_hash(hf(key));
}

// does not use or check for locks
template <class K, class H>
qf::ContainsResult standard_qfilter_conc<K, H>::unsafe_contains_hash(
    const hashed_type& hashed) const
{
    const auto [q, r]        = this->get_quotient_and_remainder(hashed);
    const entry_pos  q_pos   = quotient_position(q);
    const entry_type q_entry = table[q_pos.first].entry(cc, q_pos.second);

    if (cell_type::entry_is_empty(cc, q_entry))
        return ContainsResultEnum::CANONICAL_SLOT_EMPTY;
    if (!cell_type::entry_is_occupied(cc, q_entry))
        return ContainsResultEnum::ENTRY_NOT_PRESENT;

    auto iter = unsafe_find_run_start(find_cluster_start(q_pos), q_pos);

    entry_type current_entry = table[iter.first].entry(cc, iter.second);

    do {
        if (cell_type::entry_remainder(cc, current_entry) == r)
            return ContainsResultEnum::ENTRY_PRESENT;

        increment(iter);
        current_entry = table[iter.first].entry(cc, iter.second);
    } while (cell_type::entry_is_continuation(cc, current_entry));

    return ContainsResultEnum::ENTRY_NOT_PRESENT;
}


// *** BASIC FUNCTIONS *********************************************************

template <class K, class H>
std::pair<size_t, typename standard_qfilter_conc<K, H>::entry_type>
standard_qfilter_conc<K, H>::get_quotient_and_remainder(
    const hashed_type& hashed) const
{
    return qf::get_quotient_and_remainder<entry_type>(
        hashed, this->quotient_bits, cc.remainder_bits);
}


template <class K, class H> size_t standard_qfilter_conc<K, H>::capacity() const
{
    return ONE << quotient_bits;
}


template <class K, class H>
void standard_qfilter_conc<K, H>::grow(base_qf& target_base, size_t block_size)
{
    if (cc.remainder_bits <= 1) { return; }
    else
    {
        if (this->table_status != qf::TableStatus::growing) return;

        auto& target = dynamic_cast<next_bigger_qf&>(target_base);

        const size_t total_blocks =
            (table_capacity + block_size - 1) / block_size;
        size_t block = growing_block_counter++;

        while (block < total_blocks)
        {
            grow(target, block * block_size, block_size);
            finished_block_counter++;
            block = growing_block_counter++;
        }

        // wait for all threads to finish their blocks
        while (finished_block_counter != total_blocks) {}
    }
}




// *****************************************************************************
// *** PRIVATE HELPER FUNCTIONS ************************************************
// *****************************************************************************

// *** SOME POSITION STUFF *****************************************************

template <class K, class H>
typename standard_qfilter_conc<K, H>::entry_pos
standard_qfilter_conc<K, H>::quotient_position(size_t quotient) const
{
    return {quotient / cc.capacity, quotient % cc.capacity};
}

template <class K, class H>
void standard_qfilter_conc<K, H>::increment(entry_pos& pos) const
{
    pos.second++;
    if (pos.second == cc.capacity)
    {
        pos.second = 0;
        pos.first++;
    }
}

template <class K, class H>
void standard_qfilter_conc<K, H>::decrement(entry_pos& pos) const
{
    if (pos.second == 0)
    {
        pos.second = cc.capacity;
        ;
        pos.first--;
    }
    pos.second--;
}

template <class K, class H>
size_t standard_qfilter_conc<K, H>::offset(const entry_pos& pos) const
{
    return pos.first * cc.capacity + pos.second;
}

template <class K, class H>
size_t standard_qfilter_conc<K, H>::distance(const entry_pos& start,
                                             const entry_pos& end) const
{
    return static_cast<int>(offset(end)) - static_cast<int>(offset(start));
}





template <class K, class H>
ContainsResultEnum standard_qfilter_conc<K, H>::try_contains_lock_elision(
    const cell_type cell, const entry_pos pos, const entry_type remainder) const
{
    // ensure no locks exist to the right
    size_t iter = pos.second;
    do {
        if (cell.is_locked(cc, iter))
            return ContainsResultEnum::FAILED_OPERATION;

        iter++;
    } while (iter < cc.capacity && !cell.is_empty(cc, iter));

    // find cluster start && ensure no locks exist to the left
    size_t cluster_start = pos.second;

    while (!cell.is_cluster_start(cc, cluster_start))
    {
        if (cluster_start == 0 || cell.is_read_locked(cc, cluster_start))
            return ContainsResultEnum::FAILED_OPERATION;
        --cluster_start;
    }

    // find run start
    size_t bucket = cluster_start;
    size_t run    = bucket;

    while (bucket != pos.second)
    {
        do {
            ++run;
            if (run == cc.capacity) return ContainsResultEnum::FAILED_OPERATION;

        } while (cell.is_continuation(cc, run));

        do {
            ++bucket;
            if (bucket == cc.capacity)
                return ContainsResultEnum::FAILED_OPERATION;

        } while (!cell.is_occupied(cc, bucket));
    }

    // search run
    do {
        if (cell.remainder(cc, run) == remainder)
            return ContainsResultEnum::ENTRY_PRESENT;

        if (cell.remainder(cc, run) > remainder)
            return ContainsResultEnum::ENTRY_NOT_PRESENT;

        ++run;
        if (run == cc.capacity) return ContainsResultEnum::FAILED_OPERATION;

    } while (cell.is_continuation(cc, run));

    return ContainsResultEnum::ENTRY_NOT_PRESENT;
}


template <class K, class H>
qf::InsertResult standard_qfilter_conc<K, H>::try_insert_lock_elision(
    cell_type cell, const entry_pos pos, const entry_type remainder) const
{
    auto original_cell = cell;
    if (cell.is_locked(cc, pos.second)) return qf::FAILED_OPERATION;

    // find empty entry && ensure no locks exist to the right
    size_t empty_pos = pos.second;

    do {
        ++empty_pos;
        if (empty_pos == cc.capacity || cell.is_locked(cc, empty_pos))
            return qf::FAILED_OPERATION;
    } while (!cell.is_empty(cc, empty_pos));

    // find cluster start && ensure no locks exist to the left
    size_t cluster_start = pos.second;

    while (!cell.is_cluster_start(cc, cluster_start))
    {
        if (cluster_start == 0 || cell.is_read_locked(cc, cluster_start))
            return qf::FAILED_OPERATION;
        --cluster_start;
    }

    const bool run_already_existing = cell.is_occupied(cc, pos.second);
    if (!run_already_existing) cell.set_occupied(cc, true, pos.second);

    // find run start
    size_t bucket = cluster_start;
    size_t run    = bucket;

    while (bucket != pos.second)
    {
        do {
            ++run;
        } while (cell.is_continuation(cc, run));

        do {
            ++bucket;
        } while (!cell.is_occupied(cc, bucket));
    }

    // find insert position in run
    size_t insert_pos = run;

    entry_type new_entry = 0;
    cell_type::entry_set_remainder(cc, new_entry, remainder);

    if (run_already_existing)
    {
        do {
            if (cell.remainder(cc, insert_pos) == remainder)
                return insert_pos - cluster_start; // element already present

            if (cell.remainder(cc, insert_pos) > remainder) break;

            ++insert_pos;
            if (insert_pos == cc.capacity) return qf::FAILED_OPERATION;

        } while (cell.is_continuation(cc, insert_pos));

        if (insert_pos == run)
        {
            auto tmp_remainder = cell.remainder(cc, run);
            cell.set_remainder(cc, remainder, run);
            cell_type::entry_set_remainder(cc, new_entry, tmp_remainder);

            ++insert_pos;
            if (insert_pos == cc.capacity) return qf::FAILED_OPERATION;
        }

        cell_type::entry_set_continuation(cc, new_entry, true);
    }

    cell_type::entry_set_shifted(cc, new_entry, insert_pos != pos.second);

    // insert and shift
    cell_type shifted_cell;

    // write unshifted entries at beginning of cell
    for (size_t i = 0; i < insert_pos; i++)
    {
        shifted_cell.set_entry(cc, cell.entry(cc, i), i);
    }

    // write new entry
    shifted_cell.set_entry(cc, new_entry, insert_pos);

    // write shifted entries
    for (size_t i = insert_pos; i < empty_pos; i++)
    {
        shifted_cell.set_entry(cc, cell.entry(cc, i), i + 1);
        shifted_cell.set_shifted(cc, true, i + 1);
    }

    // write unshifted entries at end of cell
    for (size_t i = empty_pos + 1; i < cc.capacity; i++)
    {
        shifted_cell.set_entry(cc, cell.entry(cc, i), i);
    }

    shifted_cell.set_occupied_full(cc, cell.occupied_mask_full(cc));

    if (!table[pos.first].CAS(original_cell, shifted_cell))
        return qf::FAILED_OPERATION;

    return insert_pos - cluster_start;
}



template <class K, class H>
void standard_qfilter_conc<K, H>::increment_and_update(entry_pos&  pos,
                                                       cell_type&  cell,
                                                       entry_type& entry) const
{
    pos.second++;
    if (pos.second == cc.capacity)
    {
        pos.second = 0;
        pos.first++;
        cell = table[pos.first];
    }
    entry = cell.entry(cc, pos.second);
}

template <class K, class H>
void standard_qfilter_conc<K, H>::decrement_and_update(entry_pos&  pos,
                                                       cell_type&  cell,
                                                       entry_type& entry) const
{
    if (pos.second == 0)
    {
        pos.second = cc.capacity;
        ;
        pos.first--;
        cell = table[pos.first];
    }
    pos.second--;
    entry = cell.entry(cc, pos.second);
}

template <class K, class H>
void standard_qfilter_conc<K, H>::force_update(const entry_pos& pos,
                                               cell_type&       cell,
                                               entry_type&      entry) const
{
    cell  = table[pos.first];
    entry = cell.entry(cc, pos.second);
}

template <class K, class H>
bool standard_qfilter_conc<K, H>::try_empty_slot(cell_type  local_cell,
                                                 entry_pos  pos,
                                                 entry_type remainder)
{
    cell_type filled_cell;

    while (true)
    {
        while (local_cell.is_write_locked(cc, pos.second))
        {
            if (local_cell.remainder(cc, pos.second) == 1)
                return false; // write lock set by growing is final

            local_cell = table[pos.first];
        }

        if (!local_cell.is_empty(cc, pos.second)) return false;

        filled_cell = local_cell;
        filled_cell.set_remainder(cc, remainder, pos.second);
        filled_cell.set_occupied(cc, true, pos.second);

        if (table[pos.first].CAS(local_cell, filled_cell)) return true;

        local_cell = table[pos.first];
    }
}

// does not use or check for locks
template <class K, class H>
bool standard_qfilter_conc<K, H>::unsafe_try_empty_slot(entry_pos  pos,
                                                        entry_type remainder)
{
    cell_type local_cell = table[pos.first];
    cell_type filled_cell;

    while (local_cell.is_empty(cc, pos.second))
    {
        filled_cell = local_cell;
        filled_cell.set_remainder(cc, remainder, pos.second);
        filled_cell.set_occupied(cc, true, pos.second);

        if (table[pos.first].CAS(local_cell, filled_cell)) return true;

        local_cell = table[pos.first];
    }

    return false;
}

// assumes that pos is not empty
template <class K, class H>
bool standard_qfilter_conc<K, H>::has_space_and_lock(
    entry_pos pos, write_lock_guard& write_lock)
{
    increment(pos);

    cell_type  local_cell  = table[pos.first];
    entry_type local_entry = local_cell.entry(cc, pos.second);
    cell_type  locked_cell;

    while (pos.first != table_capacity - 1)
    {
        while (cell_type::entry_is_write_locked(cc, local_entry))
        {
            if (local_cell.remainder(cc, pos.second) == 1)
                return false; // write lock set by growing is final

            force_update(pos, local_cell, local_entry);
        }

        if (local_cell.is_empty(cc, pos.second))
        {
            // entry was empty and not locked -> acquire write lock
            locked_cell = local_cell;
            locked_cell.set_write_lock(cc, pos.second);

            if (table[pos.first].CAS(local_cell, locked_cell))
            {
                write_lock.set_data(&table[pos.first], pos,
                                    local_cell.status(cc, pos.second));

                if (this->table_status != qf::TableStatus::active)
                {
                    write_lock.release();
                    return false;
                }

                return true;
            }

            force_update(pos, local_cell, local_entry);
        }
        else
        {
            increment_and_update(pos, local_cell, local_entry);
        }
    }

    return false;
}

// assumes that pos is not empty
template <class K, class H>
void standard_qfilter_conc<K, H>::get_growing_write_lock(
    entry_pos pos, write_lock_guard& write_lock)
{
    cell_type  local_cell;
    cell_type  locked_cell;
    entry_type locked_entry = 0;
    cell_type::entry_set_status(cc, locked_entry, cell_type::write_lock_status);
    cell_type::entry_set_remainder(cc, locked_entry, 1u);

    increment(pos);

    while (true)
    {
        do {
            local_cell = table[pos.first];
        } while (local_cell.is_write_locked(cc, pos.second));

        if (local_cell.is_empty(cc, pos.second))
        {
            locked_cell = local_cell;
            locked_cell.set_entry(cc, locked_entry, pos.second);

            if (table[pos.first].CAS(local_cell, locked_cell))
            {
                write_lock.set_data(&table[pos.first], pos,
                                    cell_type::write_lock_status);
                return;
            }
        }
        else
        {
            increment(pos);
        }
    }
}

template <class K, class H>
typename standard_qfilter_conc<K, H>::entry_pos
standard_qfilter_conc<K, H>::find_and_lock_cluster_start(
    entry_pos pos, read_lock_guard& read_lock)
{
    cell_type  local_cell  = table[pos.first];
    entry_type local_entry = local_cell.entry(cc, pos.second);
    cell_type  locked_cell;

    while (true)
    {
        while (cell_type::entry_is_read_locked(cc, local_entry))
        {
            force_update(pos, local_cell, local_entry);
        }

        if (cell_type::entry_is_cluster_start(cc, local_entry))
        {
            // entry was cluster start and not locked -> acquire read lock
            locked_cell = local_cell;
            locked_cell.set_read_lock(cc, pos.second);

            if (table[pos.first].CAS(local_cell, locked_cell))
            {
                read_lock.set_data(&table[pos.first], pos,
                                   local_cell.status(cc, pos.second));
                return pos;
            }

            force_update(pos, local_cell, local_entry);
        }
        else
        {
            decrement_and_update(pos, local_cell, local_entry);
        }
    }
}


template <class K, class H>
typename standard_qfilter_conc<K, H>::entry_pos
standard_qfilter_conc<K, H>::find_cluster_start(entry_pos pos) const
{
    cell_type local_cell = table[pos.first];

    while (!local_cell.is_cluster_start(cc, pos.second))
    {
        decrement(pos);
        local_cell = table[pos.first];
    }

    return pos;
}


template <class K, class H>
typename standard_qfilter_conc<K, H>::entry_pos
standard_qfilter_conc<K, H>::find_run_start(entry_pos cluster_start,
                                            entry_pos pos,
                                            bool      set_occupied) const
{
    entry_pos bucket = cluster_start;
    entry_pos run    = bucket;

    if (set_occupied) table[pos.first].set_occupied(cc, true, pos.second);

    entry_type local_entry;

    while (bucket != pos)
    {
        do {
            increment(run);

            do {
                local_entry = table[run.first].entry(cc, run.second);
            } while (cell_type::entry_is_read_locked(cc, local_entry));
        } while (!cell_type::entry_is_write_locked(cc, local_entry) &&
                 cell_type::entry_is_continuation(cc, local_entry));

        do {
            increment(bucket);
            local_entry = table[bucket.first].entry(cc, bucket.second);
        } while (!cell_type::entry_is_occupied(cc, local_entry));
    }

    return run;
}

// does not use or check for locks
template <class K, class H>
typename standard_qfilter_conc<K, H>::entry_pos
standard_qfilter_conc<K, H>::unsafe_find_run_start(entry_pos cluster_start,
                                                   entry_pos pos) const
{
    entry_pos bucket = cluster_start;
    entry_pos run    = bucket;

    entry_type local_entry;

    while (bucket != pos)
    {
        do {
            increment(run);
            local_entry = table[run.first].entry(cc, run.second);
        } while (cell_type::entry_is_continuation(cc, local_entry));

        do {
            increment(bucket);
            local_entry = table[bucket.first].entry(cc, bucket.second);
        } while (!cell_type::entry_is_occupied(cc, local_entry));
    }

    return run;
}

template <class K, class H>
bool standard_qfilter_conc<K, H>::acquire_read_lock(entry_pos        pos,
                                                    read_lock_guard& read_lock)
{
    if (read_lock.locked() && read_lock.pos == pos) return true;

    cell_type local_cell;
    cell_type locked_cell;

    while (true)
    {
        do {
            local_cell = table[pos.first];
        } while (local_cell.is_read_locked(cc, pos.second));

        if (!local_cell.is_cluster_start(cc, pos.second)) { return false; }

        locked_cell = local_cell;
        locked_cell.set_read_lock(cc, pos.second);

        if (table[pos.first].CAS(local_cell, locked_cell))
        {
            read_lock.set_data(&table[pos.first], pos,
                               local_cell.status(cc, pos.second));
            return true;
        }
    }
}

// INSERT HELPER FUNCTIONS
template <class K, class H>
void standard_qfilter_conc<K, H>::insert_and_shift(entry_type        entry,
                                                   entry_pos&        insert_pos,
                                                   write_lock_guard& write_lock,
                                                   read_lock_guard&  read_lock)
{
    if (insert_pos == read_lock.pos)
    {
        // note: read lock is not on the run start of entry's canonical run

        const entry_type remainder =
            table[insert_pos.first].remainder(cc, insert_pos.second);
        table[insert_pos.first].set_remainder(
            cc, cell_type::entry_remainder(cc, entry), insert_pos.second);

        const auto status = read_lock.status;
        read_lock.status  = cell_type::entry_status(cc, entry);
        cell_type::status_set_occupied(
            cc, read_lock.status, cell_type::entry_is_occupied(cc, status));

        cell_type::entry_set_remainder(cc, entry, remainder);
        cell_type::entry_set_status(cc, entry,
                                    cell_type::entry_status(cc, status));
        cell_type::entry_set_shifted(cc, entry, true);
        cell_type::entry_set_occupied(cc, entry, false);

        increment(insert_pos);
    }

    while (insert_pos.first != write_lock.pos.first)
    {
        insert_entry_full(entry, insert_pos);
        insert_pos.first++;
        insert_pos.second = 0;
    };

    insert_entry_not_full(entry, insert_pos, write_lock.pos.second);
    write_lock.clear();
}

// inserts entry at pos and shifts the following entries cell at pos
// must have enough space for the shift i.e. an empty entry (write lock)
template <class K, class H>
void standard_qfilter_conc<K, H>::insert_entry_not_full(entry_type entry,
                                                        entry_pos  insert_pos,
                                                        size_t write_lock_pos)
{
    cell_type cell;
    cell_type shifted_cell;

    const size_t read_lock_check_start = insert_pos.second;

    do {
        shifted_cell.clear();

        do {
            cell = table[insert_pos.first];
        } while (
            cell.is_read_locked(cc, read_lock_check_start, write_lock_pos));

        // write unshifted entries at beginning of cell
        for (size_t i = 0; i < insert_pos.second; i++)
        {
            shifted_cell.set_entry(cc, cell.entry(cc, i), i);
        }

        // write new entry
        shifted_cell.set_entry(cc, entry, insert_pos.second);

        // write shifted entries
        for (size_t i = insert_pos.second; i < write_lock_pos; i++)
        {
            shifted_cell.set_entry(cc, cell.entry(cc, i), i + 1);
            shifted_cell.set_shifted(cc, true, i + 1);
        }

        // write unshifted entries at end of cell
        for (size_t i = write_lock_pos + 1; i < cc.capacity; i++)
        {
            shifted_cell.set_entry(cc, cell.entry(cc, i), i);
        }

        shifted_cell.set_occupied_full(cc, cell.occupied_mask_full(cc));
    } while (!table[insert_pos.first].CAS(cell, shifted_cell));
}

// inserts entry at pos and shifts the following entries the last
// entry which is lost due to the shift is written back to the entry reference
template <class K, class H>
void standard_qfilter_conc<K, H>::insert_entry_full(entry_type& entry,
                                                    entry_pos   insert_pos)
{
    cell_type cell;
    cell_type shifted_cell;

    const size_t read_lock_check_start = insert_pos.second;

    do {
        shifted_cell.clear();

        do {
            cell = table[insert_pos.first];
        } while (cell.is_read_locked(cc, read_lock_check_start, cc.capacity));

        // write unshifted entries at beginning of cell
        for (size_t i = 0; i < insert_pos.second; i++)
        {
            shifted_cell.set_entry(cc, cell.entry(cc, i), i);
        }

        // write new entry
        shifted_cell.set_entry(cc, entry, insert_pos.second);

        // write shifted entries
        for (size_t i = insert_pos.second + 1; i < cc.capacity; i++)
        {
            shifted_cell.set_entry(cc, cell.entry(cc, i - 1), i);
            shifted_cell.set_shifted(cc, true, i);
        }

        shifted_cell.set_occupied_full(cc, cell.occupied_mask_full(cc));
    } while (!table[insert_pos.first].CAS(cell, shifted_cell));

    entry = cell.entry(cc, cc.capacity - 1);
    cell_type::entry_set_shifted(cc, entry, true);
}

template <class K, class H>
typename standard_qfilter_conc<K, H>::entry_type
standard_qfilter_conc<K, H>::wait_on_insert_write_lock(entry_pos pos) const
{
    cell_type local_cell;

    do {
        local_cell = table[pos.first];
    } while (local_cell.is_write_locked(cc, pos.second) &&
             local_cell.remainder(cc, pos.second) == 0);

    return local_cell.entry(cc, pos.second);
}


// *****************************************************************************
// *** GROWING FUNCTION ********************************************************
// *****************************************************************************
template <class K, class H>
void standard_qfilter_conc<K, H>::grow(next_bigger_qf& target,
                                       size_t block_start, size_t block_size)
{
    read_lock_guard  read_lock(cc);
    write_lock_guard write_lock(cc);
    entry_pos        cluster_start{block_start, 0};
    const entry_pos  block_end{
        std::min(block_start + block_size, table_capacity), 0};
    entry_type local_entry =
        table[cluster_start.first].entry(cc, cluster_start.second);

    // skip write lock
    if (cell_type::entry_is_write_locked(cc, local_entry))
    {
        wait_on_insert_write_lock(cluster_start);
        increment(cluster_start);
        local_entry =
            table[cluster_start.first].entry(cc, cluster_start.second);
    }

    // skip partial cluster
    if (!cell_type::entry_is_empty(cc, local_entry) &&
        cluster_start != entry_pos{0, 0})
    {
        auto prev_pos = cluster_start;
        decrement(prev_pos);
        const entry_type prev_entry = wait_on_insert_write_lock(prev_pos);

        if (!cell_type::entry_is_empty(cc, prev_entry) &&
            !cell_type::entry_is_write_locked(cc, prev_entry))
        {
            // cluster start is in "middle" of cluster -> skip cluster
            do {
                increment(cluster_start);
                local_entry = wait_on_insert_write_lock(cluster_start);

                if (cell_type::entry_is_write_locked(cc, local_entry))
                {
                    // can only be a write lock from growing
                    increment(cluster_start);
                    break;
                }
            } while (cluster_start < block_end &&
                     !cell_type::entry_is_empty(cc, local_entry));
        }
    }

    // skip empty entries
    while (cluster_start < block_end &&
           table[cluster_start.first].is_empty(cc, cluster_start.second))
    {
        increment(cluster_start);
    }

    if (cluster_start >= block_end) return;

    do {
        get_growing_write_lock(cluster_start, write_lock);
        acquire_read_lock(cluster_start, read_lock);
        cluster_start = grow_cluster(target, read_lock);

    } while (cluster_start < block_end);
}

template <class K, class H>
typename standard_qfilter_conc<K, H>::entry_pos
standard_qfilter_conc<K, H>::grow_cluster(
    next_bigger_qf& target, const read_lock_guard& read_lock) const
{
    size_t remainder_bits_trgt = target.cc.remainder_bits;
    using cell_type_trgt       = typename next_bigger_qf::cell_type;
    using entry_pos_trgt       = typename next_bigger_qf::entry_pos;
    using entry_type_trgt      = typename next_bigger_qf::entry_type;

    entry_type_trgt target_remainder_mask = (1ull << remainder_bits_trgt) - 1;
    auto            run_start             = read_lock.pos;
    auto            run_iter              = read_lock.pos;

    entry_type      source_entry, tmp_entry;
    entry_type_trgt target_entry = 0;
    bool            leading_bit, last_leading_bit;
    entry_pos_trgt  target_pos, last_target_pos;

    auto read = [&]() {
        do {
            source_entry = table[run_iter.first].entry(cc, run_iter.second);
        } while (run_iter != read_lock.pos &&
                 cell_type::entry_is_read_locked(cc, source_entry));

        if (cell_type::entry_is_read_locked(cc, source_entry))
            cell_type::entry_set_status(cc, source_entry, read_lock.status);
        else if (cell_type::entry_is_write_locked(cc, source_entry))
            cell_type::entry_set_status(cc, source_entry, 0);

        leading_bit =
            cell_type::entry_remainder(cc, source_entry) >> remainder_bits_trgt;
        cell_type_trgt::entry_set_remainder(
            target.cc, target_entry,
            cell_type::entry_remainder(cc, source_entry) &
                target_remainder_mask);
    };

    auto advance = [&]() {
        // write
        target.table[target_pos.first].set_entry(target.cc, target_entry,
                                                 target_pos.second);

        // backup state
        last_leading_bit = leading_bit;
        last_target_pos  = target_pos;

        target.increment(target_pos);
        increment(run_iter);

        read();
    };

    read();
    target_pos = target.quotient_position(offset(run_iter) * 2 +
                                          static_cast<size_t>(leading_bit));

    bool second_pass = false;

    cell_type_trgt::entry_set_occupied(target.cc, target_entry, false);

    while (!cell_type::entry_is_empty(cc, source_entry))
    {
        // move run_iter start
        cell_type_trgt::entry_set_continuation(target.cc, target_entry, false);

        if (cell_type::entry_is_cluster_start(cc, source_entry))
        {
            cell_type_trgt::entry_set_shifted(target.cc, target_entry, false);
        }
        else if (target_pos <= last_target_pos)
        {
            target_pos = last_target_pos;
            target.increment(target_pos);
            cell_type_trgt::entry_set_shifted(target.cc, target_entry, true);
        }
        else
        {
            cell_type_trgt::entry_set_shifted(target.cc, target_entry, false);
        }

        advance();

        // move continuation of run_iter
        if (cell_type::entry_is_continuation(cc, source_entry))
        {
            cell_type_trgt::entry_set_continuation(target.cc, target_entry,
                                                   true);
            cell_type_trgt::entry_set_shifted(target.cc, target_entry, true);

            while (cell_type::entry_is_continuation(cc, source_entry) &&
                   leading_bit == last_leading_bit)
            {
                advance();
            }

            if (!second_pass &&
                cell_type::entry_is_continuation(cc, source_entry))
            {
                second_pass = true;
                target_pos  = target.quotient_position(
                    offset(run_start) * 2 + static_cast<size_t>(leading_bit));
                continue;
            }
        }

        if (second_pass)
        {
            target_pos = target.quotient_position(offset(run_start) * 2);

            target.table[target_pos.first].set_occupied(target.cc, true,
                                                        target_pos.second);
            target.increment(target_pos);
            target.table[target_pos.first].set_occupied(target.cc, true,
                                                        target_pos.second);
        }
        else
        {
            target_pos = target.quotient_position(
                offset(run_start) * 2 + static_cast<size_t>(last_leading_bit));

            target.table[target_pos.first].set_occupied(target.cc, true,
                                                        target_pos.second);
        }

        second_pass = false;

        increment(run_start);
        for (; run_start.first < table_capacity; increment(run_start))
        {
            do {
                tmp_entry = table[run_start.first].entry(cc, run_start.second);
            } while (run_start != read_lock.pos &&
                     cell_type::entry_is_read_locked(cc, tmp_entry));

            if (cell_type::entry_is_read_locked(cc, tmp_entry))
                cell_type::entry_set_status(cc, tmp_entry, read_lock.status);

            if (!cell_type::entry_is_write_locked(cc, tmp_entry) &&
                cell_type::entry_is_occupied(cc, tmp_entry))
                break;
        }

        // target_pos = entry_pos_trgt(offset(run_start) * 2
        //                             + static_cast<size_t>(leading_bit));

        target_pos = target.quotient_position(offset(run_start) * 2 +
                                              static_cast<size_t>(leading_bit));
    }

    return run_start;
}


// *****************************************************************************
// *** PRINT AND CHECKER FUNCTIONS *********************************************
// *****************************************************************************
template <class K, class H> void standard_qfilter_conc<K, H>::print() const
{
    print(0, table_capacity);
}

template <class K, class H>
void standard_qfilter_conc<K, H>::print(entry_pos start, entry_pos end) const
{
    print(start.first, end.first - start.first + 1);
}

template <class K, class H>
void standard_qfilter_conc<K, H>::print(size_t start, size_t count) const
{
    std::cout << "occu | cont | shif || qoutient | remainder\n";
    const size_t end = start + count;

    for (size_t table_pos = start; table_pos < end; table_pos++)
    {
        print_cell(table_pos);
        std::cout << "\n";
    }

    std::cout << std::endl;
}

template <class K, class H>
void standard_qfilter_conc<K, H>::print_cell(size_t table_pos) const
{
    print_cell(table[table_pos], table_pos);
}

template <class K, class H>
void standard_qfilter_conc<K, H>::print_cell(cell_type cell,
                                             size_t    table_pos) const
{
    for (size_t cell_pos = 0; cell_pos < cc.capacity; cell_pos++)
    {
        print_entry(cell.entry(cc, cell_pos), {table_pos, cell_pos});
        std::cout << "   ";
    }
}

template <class K, class H>
void standard_qfilter_conc<K, H>::print_entry(entry_pos pos) const
{
    print_entry(table[pos.first].entry(cc, pos.second), offset(pos));
}

template <class K, class H>
void standard_qfilter_conc<K, H>::print_entry(entry_type entry,
                                              entry_pos  pos) const
{
    print_entry(entry, offset(pos));
}

template <class K, class H>
void standard_qfilter_conc<K, H>::print_entry(entry_type entry,
                                              size_t     index) const
{
    std::cout << cell_type::entry_is_occupied(cc, entry) << "|"
              << cell_type::entry_is_continuation(cc, entry) << "|"
              << cell_type::entry_is_shifted(cc, entry) << " || "
              << qf::bitstring(index, this->quotient_bits) << "|"
              << qf::bitstring(cell_type::entry_remainder(cc, entry),
                               cc.remainder_bits);
}

template <class K, class H>
void standard_qfilter_conc<K, H>::print_entry(entry_type entry) const
{
    std::cout << cell_type::entry_is_occupied(cc, entry) << "|"
              << cell_type::entry_is_continuation(cc, entry) << "|"
              << cell_type::entry_is_shifted(cc, entry) << " ||     "
              << qf::bitstring(cell_type::entry_remainder(cc, entry),
                               cc.remainder_bits);
}

template <class K, class H>
bool standard_qfilter_conc<K, H>::check_consistency() const
{
    bool       consistent = true;
    bool       sorted     = true;
    entry_type last       = 0;
    entry_type current;

    for (entry_pos pos{0, 0}; pos.first < table_capacity; increment(pos))
    {
        current = table[pos.first].entry(cc, pos.second);

        /*if (cell_type::entry_is_read_locked(current))
          {
          print_entry(current);
          std::cerr << "-> error: " << pos << " read locked" << std::endl;
          consistent = false;
          }
          else if (cell_type::entry_is_write_locked(current))
          {
          print_entry(current);
          std::cerr << "-> error: " << pos << " write locked" << std::endl;
          consistent = false;
          }
          else */
        if (cell_type::entry_is_continuation(cc, current) &&
            !cell_type::entry_is_shifted(cc, current))
        {
            print_entry(current);
            std::cerr << "-> error: [" << pos.first << "," << pos.second
                      << "] cont =/=> shift" << std::endl;
            consistent = false;
        }
        else if (cell_type::entry_is_continuation(cc, current) &&
                 !cell_type::entry_is_shifted(cc, last) &&
                 !cell_type::entry_is_cluster_start(cc, last))
        {
            print_entry(current);
            std::cerr << "-> error2: [" << pos.first << "," << pos.second
                      << "] cont =/=> shift" << std::endl;
            consistent = false;
        }
        else if (cell_type::entry_is_continuation(cc, current))
        {
            if (cell_type::entry_remainder(cc, current) <
                cell_type::entry_remainder(cc, last))
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

} // namespace qf
