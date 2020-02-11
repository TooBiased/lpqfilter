#pragma once
#include <iostream>
#include <utility>
#include <memory>
#include "utils/utilities.h"
#include "standard_qfilter_cell.hpp"


namespace qf {

template <class Key, class Hash>
	class standard_qfilter_seq;

template <class Key, class Hash = default_hash>
class standard_qfilter_seq_base
{
public:
	using key_type = Key;
	using hash_function_type = Hash;
	using hashed_type = decltype(std::declval<Hash>()(std::declval<Key>()));

	virtual ~standard_qfilter_seq_base() = default;

	virtual qf::InsertResult insert(const key_type& key) = 0;
	virtual qf::InsertResult insert_hash(const hashed_type& hashed) = 0;
	virtual bool quick_insert(const key_type& key) = 0;
	virtual bool quick_insert_hash(const hashed_type& hashed) = 0;
	virtual qf::ContainsResult contains(const key_type& key) const = 0;
	virtual qf::ContainsResult contains_hash(const hashed_type& hashed) const=0;
	virtual size_t capacity() const = 0;
	virtual size_t memory_usage_bytes() const = 0;
	virtual size_t unused_memory_bits() const = 0;
	virtual double fill_level() const = 0;
	virtual bool check_consistency() = 0;
	virtual standard_qfilter_seq_base* create_bigger_QF(
        const hash_function_type& hf = hash_function_type()) const = 0;
	virtual void grow(standard_qfilter_seq_base<Key, Hash>& target_base) = 0;
};


template <class Key, class Hash>
class standard_qfilter_locking;

template <class Key,
          //size_t remainder_bits,
	class Hash = default_hash>
	class standard_qfilter_seq : public standard_qfilter_seq_base<Key, Hash>
{
public:
    static constexpr bool is_growing_compatible = true;
    static constexpr bool is_templated  = false;
    static constexpr bool is_sequential = true;
    static constexpr bool is_growing    = false;
    static constexpr bool is_dynamic    = false;
    static constexpr bool uses_handle   = false;

    using this_type = standard_qfilter_seq<Key, Hash>;
	using key_type = Key;
	using hash_function_type = Hash;
	using hashed_type = decltype(std::declval<Hash>()(std::declval<Key>()));
	using next_smaller_qf = standard_qfilter_seq<Key, Hash>;
	using next_bigger_qf  = standard_qfilter_seq<Key, Hash>;
	using base_qf = standard_qfilter_seq_base<Key, Hash>;

	friend base_qf;
    friend standard_qfilter_locking<Key, Hash>;
	// friend next_smaller_qf;

private:
	using cell_base_type = qf::grouped_cell_base_type;
	using cell_type = qf::standard_cell_non_atomic<cell_base_type>;
	using entry_type = typename cell_type::entry_type;
    using cell_stuff = typename cell_type::cstuff;
	using entry_pos = std::pair<size_t, size_t>;

    const cell_stuff             cc;
	const size_t                 quotient_bits;
	const size_t                 table_capacity;
	const size_t                 fingerprint_mask;
	std::unique_ptr<cell_type[]> table;
	const hash_function_type     hf;

public:
    // CONSTRUCTORS ************************************************************
	standard_qfilter_seq(size_t capacity = qf::DEFAULT_MIN_CAPACITY,
                         size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS,
                         const hash_function_type& hf = hash_function_type());

	standard_qfilter_seq(const standard_qfilter_seq& other);
	standard_qfilter_seq& operator=(const standard_qfilter_seq& other);

	standard_qfilter_seq(standard_qfilter_seq&& other) = default;
	standard_qfilter_seq& operator=(standard_qfilter_seq&& other) = default;

	~standard_qfilter_seq() = default;

    // BASIC TABLE FUNCTIONALITY ***********************************************
	this_type* create_bigger_QF(
        const hash_function_type& hash = hash_function_type()) const override;
	qf::InsertResult insert(const key_type& key) override;
	qf::InsertResult insert_hash(const hashed_type& hashed) override;
	bool quick_insert(const key_type& key) override;
    bool quick_insert_hash(const hashed_type& hashed) override;
    qf::ContainsResult contains(const key_type& key) const override;
    qf::ContainsResult contains_hash(const hashed_type& hashed) const override;

	std::pair<size_t, entry_type> get_quotient_and_remainder(
        const hashed_type& hashed) const;

    // STATS *******************************************************************
	size_t capacity()           const override;
	size_t memory_usage_bytes() const override;
	size_t unused_memory_bits() const override;
	double fill_level() const override;

	virtual void grow(base_qf& target_base) override;

private:
    // PRIVATE HELPER FUNCTIONS ************************************************
    entry_pos quotient_position(size_t quotient) const;
    void increment(entry_pos& pos) const;
    void decrement(entry_pos& pos) const;
    size_t offset(const entry_pos& pos) const;
    size_t distance(const entry_pos& start, const entry_pos& end) const;

    entry_pos find_cluster_start(entry_pos pos) const;
	entry_pos find_run_start(entry_pos cluster_start, entry_pos q) const;

	// assumes that pos is not empty
	bool has_space_for_shift(entry_pos pos) const;

	// inserts entry at pos and shifts the following entries cell at pos must
    // have enough space for the shift i.e. an empty entry
	void insert_entry_not_full(entry_type entry, entry_pos insert_pos,
                               size_t empty_pos);

	// inserts entry at pos and shifts the following entries the last entry
    // which is lost due to the shift is written back to the entry reference
	void insert_entry_full(entry_type& entry, entry_pos insert_pos);

	void insert_and_shift(entry_type entry, entry_pos& insert_pos);

	entry_pos grow_cluster(next_bigger_qf& target, entry_pos cluster_start);

public:
    // UNNECESSARY FUNCTIONS - CHECKER AND PRINT STUFF *************************
	void print();
	void print_entry(entry_pos pos);
	void print_entry(entry_type entry, size_t index);
	void print_entry(entry_type entry);
	bool check_consistency() override;
};







// *****************************************************************************
// *** CONSTRUCTORS ************************************************************
// *****************************************************************************
template<class K, class H>
standard_qfilter_seq<K,H>::standard_qfilter_seq(size_t capacity,
                                                              size_t remainder_bits,
                                                              const hash_function_type& hf)
    : cc(remainder_bits),
      quotient_bits(qf::capacity_to_quotient_bits(capacity)),
      table_capacity((ONE << quotient_bits) / cc.capacity + qf::shift_buffer),
      fingerprint_mask((ONE << (cc.remainder_bits + quotient_bits)) - 1),
      hf(hf)
{
    table = std::make_unique<cell_type[]>(table_capacity);
}

template<class K, class H>
standard_qfilter_seq<K,H>::standard_qfilter_seq(const standard_qfilter_seq& other)
    : cc(other.cc),
      quotient_bits(other.quotient_bits),
      table_capacity(other.table_capacity),
      fingerprint_mask(other.fingerprint_mask),
      table(new cell_type[table_capacity]),
      hf(other.hf)
{
    for (size_t i = 0; i < table_capacity; i++)
    {
        table[i] = other.table[i];
    }
}

template<class K, class H>
standard_qfilter_seq<K,H>&
standard_qfilter_seq<K,H>::operator=(const standard_qfilter_seq& other)
{
    if (this == &other)
        return *this;

    this->~this_type();
    new (this) this_type(other);
    return *this;
}

template<class K, class H>
standard_qfilter_seq<K,H>*
standard_qfilter_seq<K,H>::create_bigger_QF(const hash_function_type& hash) const
{
    if (cc.remainder_bits > 1)
    {
        return new next_bigger_qf(ONE << (quotient_bits + 1),
                                  cc.remainder_bits-1,
                                  hash);
    }
    else
    {
        return nullptr;
    }
}




// *****************************************************************************
// *** MAIN FUNCTIONALITY ******************************************************
// *****************************************************************************
template<class K, class H>
qf::InsertResult
standard_qfilter_seq<K,H>::insert(const key_type& key)
{
    return insert_hash(hf(key));
}

template<class K, class H>
qf::InsertResult
standard_qfilter_seq<K,H>::insert_hash(const hashed_type& hashed)
{
    const auto[q, r] = this->get_quotient_and_remainder(hashed);
    const auto q_pos = quotient_position(q);

    if (table[q_pos.first].is_empty(cc, q_pos.second))
    {
        table[q_pos.first].set_remainder(cc, r, q_pos.second);
        table[q_pos.first].set_occupied (cc, true, q_pos.second);
        return 1;
    }

    const bool run_already_existing = table[q_pos.first].is_occupied(cc, q_pos.second);
    table[q_pos.first].set_occupied(cc, true, q_pos.second);

    const auto cluster_start = find_cluster_start(q_pos);
    const auto run = find_run_start(cluster_start, q_pos);

    if (run.first == table_capacity - 1)
    {
        return qf::FAILED_OPERATION;
    }

    auto insert_pos = run;
    entry_type new_entry = 0;
    cell_type::entry_set_remainder(cc, new_entry, r);

    if (run_already_existing)
    {
        // find position in existing run
        do
        {
            if (table[insert_pos.first].remainder(cc, insert_pos.second) == r)
                // element already present
                return distance(cluster_start, insert_pos);
            else if (table[insert_pos.first].remainder(cc,
                                                       insert_pos.second) > r)
                break;

            increment(insert_pos);
        } while (table[insert_pos.first].is_continuation(cc,
                                                         insert_pos.second));

        // element is not present so we have to shift
        if (!has_space_for_shift(insert_pos))
        {
            return qf::FAILED_OPERATION;
        }

        if (insert_pos == run)
        {
            // old run start will be shifted
            table[run.first].set_continuation(cc, true, run.second);
        }
        else
        {
            // new entry is not run start
            cell_type::entry_set_continuation(cc, new_entry, true);
        }
    }
    else if (!has_space_for_shift(insert_pos))
    {
        table[q_pos.first].set_occupied(cc, false, q_pos.second);
        return qf::FAILED_OPERATION;
    }

    if (insert_pos != q_pos)
    {
        cell_type::entry_set_shifted(cc, new_entry, true);
    }

    insert_and_shift(new_entry, insert_pos);
    return distance(cluster_start, insert_pos);
}

template<class K, class H>
bool
standard_qfilter_seq<K,H>::quick_insert(const key_type& key)
{
    return quick_insert_hash(hf(key));
}

template<class K, class H>
bool
standard_qfilter_seq<K,H>::quick_insert_hash(const hashed_type& hashed)
{
    const auto[q, r] = this->get_quotient_and_remainder(hashed);
    const auto q_pos = quotient_position(q);

    if (table[q_pos.first].is_empty(cc, q_pos.second))
    {
        table[q_pos.first].set_remainder(cc, r, q_pos.second);
        table[q_pos.first].set_occupied(cc, true, q_pos.second);
        return true;
    }

    return contains_hash(hashed);
}

template<class K, class H>
qf::ContainsResult
standard_qfilter_seq<K,H>::contains(const key_type& key) const
{
    return contains_hash(hf(key));
}

template<class K, class H>
qf::ContainsResult
standard_qfilter_seq<K,H>::contains_hash(const hashed_type& hashed) const
{
    const auto[q, r] = this->get_quotient_and_remainder(hashed);
    const auto q_pos = quotient_position(q);
    const entry_type q_entry = table[q_pos.first].entry(cc, q_pos.second);

    if (cell_type::entry_is_empty(cc, q_entry))
        return ContainsResultEnum::CANONICAL_SLOT_EMPTY;
    if (!cell_type::entry_is_occupied(cc, q_entry))
        return ContainsResultEnum::ENTRY_NOT_PRESENT;

    auto iter = find_run_start(find_cluster_start(q_pos), q_pos);

    do
    {
        if (table[iter.first].remainder(cc, iter.second) == r)
            return ContainsResultEnum::ENTRY_PRESENT;

        increment(iter);
    } while (table[iter.first].is_continuation(cc, iter.second));

    return ContainsResultEnum::ENTRY_NOT_PRESENT;
}

template<class K, class H>
std::pair<size_t, typename standard_qfilter_seq<K,H>::entry_type>
standard_qfilter_seq<K,H>::get_quotient_and_remainder(
    const hashed_type& hashed) const
{
    return qf::get_quotient_and_remainder<entry_type>(hashed, quotient_bits,
                                                      cc.remainder_bits);
}




// *****************************************************************************
// *** STATS *******************************************************************
// *****************************************************************************
template<class K, class H>
size_t
standard_qfilter_seq<K,H>::capacity() const
{
    return ONE << quotient_bits;
}

template<class K, class H>
size_t
standard_qfilter_seq<K,H>::memory_usage_bytes() const
{
    return table_capacity * sizeof(cell_type);
}

template<class K, class H>
size_t
standard_qfilter_seq<K,H>::unused_memory_bits() const
{
    return table_capacity
           * (sizeof(cell_type) * 8
              - cc.capacity * (cc.remainder_bits + qf::status_bits));
}

template<class K, class H>
double
standard_qfilter_seq<K,H>::fill_level() const
{
    size_t used_slots = 0;

    for (entry_pos pos{ 0,0 }; pos.first < table_capacity; increment(pos))
    {
        if (!table[pos.first].is_empty(cc, pos.second))
            used_slots++;
    }

    return static_cast<double>(used_slots) / table_capacity / cc.capacity;
}

template<class K, class H>
void
standard_qfilter_seq<K,H>::grow(base_qf& target_base)
{
    if (cc.remainder_bits <= 1)
    {
        return;
    }
    else
    {
        auto& target = dynamic_cast<next_bigger_qf&>(target_base);

        entry_pos cluster_start{ 0,0 };

        while (table[cluster_start.first].is_empty(cc, cluster_start.second))
        {
            increment(cluster_start);
        }

        while (offset(cluster_start) < (ONE << quotient_bits))
        {
            cluster_start = this->grow_cluster(target, cluster_start);
        }
    }
}




// *****************************************************************************
// *** PRIVATE HELPER FUNCTIONS ************************************************
// *****************************************************************************
template<class K, class H>
typename standard_qfilter_seq<K,H>::entry_pos
standard_qfilter_seq<K,H>::quotient_position(size_t quotient) const
{
    return { static_cast<std::uint32_t>(quotient / cc.capacity),
             static_cast<std::uint8_t> (quotient % cc.capacity) };
}

template<class K, class H>
void
standard_qfilter_seq<K,H>::increment(entry_pos& pos) const
{
    pos.second++;
    if (pos.second == cc.capacity)
    {
        pos.second = 0;
        pos.first++;
    }
}

template<class K, class H>
void
standard_qfilter_seq<K,H>::decrement(entry_pos& pos) const
{
    if (pos.second == 0)
    {
        pos.second = cc.capacity;;
        pos.first--;
    }
    pos.second--;
}

template<class K, class H>
size_t
standard_qfilter_seq<K,H>::offset(const entry_pos& pos) const
{
    return pos.first*cc.capacity + pos.second;
}

template<class K, class H>
size_t
standard_qfilter_seq<K,H>::distance(const entry_pos& start,
                                           const entry_pos& end) const
{
    return static_cast<int>(offset(end)) - static_cast<int>(offset(start));
}

template<class K, class H>
typename standard_qfilter_seq<K,H>::entry_pos
standard_qfilter_seq<K,H>::find_cluster_start(entry_pos pos) const
{
    while (table[pos.first].is_shifted(cc, pos.second))
    {
        decrement(pos);
    }

    return pos;
}

template<class K, class H>
typename standard_qfilter_seq<K,H>::entry_pos
standard_qfilter_seq<K,H>::find_run_start(entry_pos cluster_start,
                                                 entry_pos q) const
{
    auto bucket = cluster_start;
    auto run = bucket;

    while (bucket != q)
    {
        do
        {
            if (run.first == table_capacity - 1)
            {
                return run; //invalid value
            }

            increment(run);
        } while (table[run.first].is_continuation(cc, run.second));

        do
        {
            increment(bucket);
        } while (!table[bucket.first].is_occupied(cc, bucket.second));
    }

    return run;
}

// assumes that pos is not empty
template<class K, class H>
bool
standard_qfilter_seq<K,H>::has_space_for_shift(entry_pos pos) const
{
    while (!table[pos.first].is_empty(cc, pos.second))
    {
        increment(pos);
    }

    return pos.first != table_capacity - 1;
}

// inserts entry at pos and shifts the following entries
// cell at pos must have enough space for the shift i.e. an empty entry
template<class K, class H>
void
standard_qfilter_seq<K,H>::insert_entry_not_full(entry_type entry,
                                                        entry_pos insert_pos,
                                                        size_t empty_pos)
{
    const cell_type cell = table[insert_pos.first];
    cell_type shifted_cell;

    // write unshifted entries at beginning of cell
    for (size_t i = 0; i < insert_pos.second; i++)
    {
        shifted_cell.set_entry(cc, cell.entry(cc, i), i);
    }

    // write new entry
    shifted_cell.set_entry(cc, entry, insert_pos.second);

    // write shifted entries
    for (size_t i = insert_pos.second + 1; i <= empty_pos; i++)
    {
        shifted_cell.set_entry(cc, cell.entry(cc, i - 1), i);
        shifted_cell.set_shifted(cc, true, i);
    }

    // write unshifted entries at end of cell
    for (size_t i = empty_pos + 1; i < cc.capacity; i++)
    {
        shifted_cell.set_entry(cc, cell.entry(cc, i), i);
    }

    shifted_cell.set_occupied_full(cc, cell.occupied_mask_full(cc));
    table[insert_pos.first] = shifted_cell;
}

// inserts entry at pos and shifts the following entries the last entry
// which is lost due to the shift is written back to the entry reference
template<class K, class H>
void
standard_qfilter_seq<K,H>::insert_entry_full(entry_type& entry,
                                                    entry_pos insert_pos)
{
    const cell_type cell = table[insert_pos.first];
    cell_type shifted_cell;

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
    table[insert_pos.first] = shifted_cell;

    entry = cell.entry(cc, cc.capacity - 1);
    cell_type::entry_set_shifted(cc, entry, true);
}

template<class K, class H>
void
standard_qfilter_seq<K,H>::insert_and_shift(entry_type entry,
                                                   entry_pos& insert_pos)
{
    entry_pos end_pos = insert_pos;
    while (!table[end_pos.first].is_empty(cc, end_pos.second))
    {
        increment(end_pos);
    }

    while (insert_pos.first != end_pos.first)
    {
        insert_entry_full(entry, insert_pos);
        insert_pos.first++;
        insert_pos.second = 0;
    };

    insert_entry_not_full(entry, insert_pos, end_pos.second);
    insert_pos = end_pos;
}

template<class K, class H>
typename standard_qfilter_seq<K,H>::entry_pos
standard_qfilter_seq<K,H>::grow_cluster(next_bigger_qf& target,
                                               entry_pos cluster_start)
{
    size_t remainder_bits_trgt = target.cc.remainder_bits;
    using cell_type_trgt       = typename next_bigger_qf::cell_type;
    using entry_pos_trgt       = typename next_bigger_qf::entry_pos;
    using entry_type_trgt      = typename next_bigger_qf::entry_type;

    entry_type_trgt target_remainder_mask = (1ull << remainder_bits_trgt) -1;
    auto run_start = cluster_start;
    auto run_iter  = cluster_start;

    entry_type source_entry;
    entry_type_trgt target_entry = 0;
    bool leading_bit, last_leading_bit;
    entry_pos_trgt target_pos, last_target_pos = {0,0};

    auto read = [&]()
		{
			source_entry = table[run_iter.first].entry(cc, run_iter.second);
			leading_bit = cell_type::entry_remainder(cc, source_entry)
                          >> remainder_bits_trgt;
			cell_type_trgt::entry_set_remainder(target.cc, target_entry,
                                    cell_type::entry_remainder(cc, source_entry)
                                     & target_remainder_mask);
		};

    auto advance = [&]()
		{
			// write
			target.table[target_pos.first].set_entry(target.cc,
                                                     target_entry,
                                                     target_pos.second);

			// backup state
			last_leading_bit = leading_bit;
			last_target_pos = target_pos;

			target.increment(target_pos);
			increment(run_iter);

			read();
		};

    read();
    target_pos = target.quotient_position(offset(run_iter) * 2
                                          + static_cast<size_t>(leading_bit));

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
            cell_type_trgt::entry_set_continuation(target.cc,
                                                   target_entry,
                                                   true);
            cell_type_trgt::entry_set_shifted(target.cc, target_entry, true);

            while (cell_type::entry_is_continuation(cc, source_entry)
                   && leading_bit == last_leading_bit)
            {
                advance();
            }

            if (!second_pass
                && cell_type::entry_is_continuation(cc, source_entry))
            {
                // leading bits differ -> new run_iter start
                second_pass = true;
                target_pos = target.quotient_position(offset(run_start) * 2
                                            + static_cast<size_t>(leading_bit));
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
            target_pos = target.quotient_position(offset(run_start) * 2
                                      + static_cast<size_t>(last_leading_bit));
            target.table[target_pos.first].set_occupied(target.cc, true,
                                                        target_pos.second);
        }

        second_pass = false;

        do
        {
            increment(run_start);
        } while (offset(run_start) < capacity()
                 && !table[run_start.first].is_occupied(cc, run_start.second));

        target_pos = target.quotient_position(offset(run_start) * 2
                                    + static_cast<size_t>(leading_bit));
    }

    return run_start;
}




// *****************************************************************************
// *** PRIVATE HELPER FUNCTIONS ************************************************
// *****************************************************************************
template<class K, class H>
void
standard_qfilter_seq<K,H>::print()
{
    std::cout << "occu | cont | shif || qoutient | remainder\n";

    for (size_t table_pos = 0; table_pos < table_capacity; table_pos++)
    {
        for (size_t cell_pos = 0; cell_pos < cc.capacity; cell_pos++)
        {
            print_entry({ table_pos, cell_pos });
            std::cout << "   ";
        }
        std::cout << "\n";
    }

    std::cout << std::endl;
}

template<class K, class H>
void
standard_qfilter_seq<K,H>::print_entry(entry_pos pos)
{
    print_entry(table[pos.first].entry(cc, pos.second), offset(pos));
}

template<class K, class H>
void
standard_qfilter_seq<K,H>::print_entry(entry_type entry, size_t index)
{
    std::cout << cell_type::entry_is_occupied(cc, entry)     << "|"
              << cell_type::entry_is_continuation(cc, entry) << "|"
              << cell_type::entry_is_shifted(cc, entry)      << " || "
              << qf::bitstring(index, quotient_bits)         << "|"
              << qf::bitstring(cell_type::entry_remainder(cc, entry),
                               cc.remainder_bits);
}

template<class K, class H>
void
standard_qfilter_seq<K,H>::print_entry(entry_type entry)
{
    std::cout << cell_type::entry_is_occupied(cc, entry)     << "|"
              << cell_type::entry_is_continuation(cc, entry) << "|"
              << cell_type::entry_is_shifted(cc, entry)      << " || "
              << "    "
              << qf::bitstring(cell_type::entry_remainder(cc, entry),
                               cc.remainder_bits);
}

template<class K, class H>
bool
standard_qfilter_seq<K,H>::check_consistency()
{
    bool consistent = true;
    bool sorted = true;
    entry_type last = 0;
    cell_type current;

    for (entry_pos pos{ 0,0 }; pos.first < table_capacity; increment(pos))
    {
        current = table[pos.first];

        if (current.is_continuation(cc, pos.second)
            && !current.is_shifted(cc, pos.second))
        {
            std::cerr << "error1: [" << pos.first << "," << pos.second
                      << "] cont =/=> shift" << std::endl;
            consistent = false;
        }
        else if (current.is_continuation(cc, pos.second)
                 && !cell_type::entry_is_shifted(cc, last)
                 && !cell_type::entry_is_cluster_start(cc, last))
        {
            std::cerr << "error2: [" << pos.first << "," << pos.second
                      << "] cont =/=> shift" << std::endl;
            consistent = false;
        }
        else if (current.is_continuation(cc, pos.second))
        {
            if (current.remainder(cc, pos.second)
                < cell_type::entry_remainder(cc, last))
            {
                sorted = false;
                consistent = false;
            }
        }

        last = current.entry(cc, pos.second);
    }

    if (!sorted)
        std::cerr << "error:  remainder not sorted" << std::endl;

    return consistent;
}

} // namespace qf
