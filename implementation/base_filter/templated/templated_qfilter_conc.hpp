#pragma once
/*******************************************************************************
 * implementation/base_filter/templated/templated_qfilter_conc.hpp
 *
 * concurrent quotient filter with templated grouped slots
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
	class templated_qfilter_conc;

template <class Key,
          class Hash = htm::default_hash>
	class templated_qfilter_conc_base
{
public:
    static constexpr bool is_growing_compatible = true;
    static constexpr bool is_templated  = true;
    static constexpr bool is_sequential = false;
    static constexpr bool is_growing    = false;
    static constexpr bool is_dynamic    = false;
    static constexpr bool uses_handle   = false;

protected:
	alignas(128) size_t                             quotient_bits;
	alignas(128) std::atomic<qf::TableStatus>       table_status;

	templated_qfilter_conc_base(const templated_qfilter_conc_base& other)
		: quotient_bits(other.quotient_bits), table_status(other.table_status)
    {}

	templated_qfilter_conc_base& operator=(const templated_qfilter_conc_base& other)
	{
		quotient_bits = other.quotient_bits;
		table_status = other.table_status.load();
		return *this;
	}

	templated_qfilter_conc_base(templated_qfilter_conc_base&& other) = default;
	templated_qfilter_conc_base& operator=(
        templated_qfilter_conc_base&& other) = default;
	templated_qfilter_conc_base(size_t capacity = qf::DEFAULT_MIN_CAPACITY,
                                qf::TableStatus table_status
                                  = qf::TableStatus::active)
	 : quotient_bits(qf::capacity_to_quotient_bits(capacity)),
                     table_status(table_status)
    {};

public:
	using key_type           = Key;
	using hash_function_type = Hash;
	using hashed_type        = decltype(std::declval<Hash>()(std::declval<Key>()));

    template <size_t rem_bits>
    using instanciated = templated_qfilter_conc<key_type,
                                                rem_bits,
                                                hash_function_type>;

	virtual ~templated_qfilter_conc_base() = default;

    virtual void prefetch(const hashed_type& hashed) const = 0;
	virtual qf::InsertResult insert(const key_type& key) = 0;
	virtual qf::InsertResult insert_hash(const hashed_type& hashed) = 0;
	virtual bool quick_insert(const key_type& key) = 0;
	virtual bool quick_insert_hash(const hashed_type& hashed) = 0;
	virtual qf::ContainsResult contains(const key_type& key) = 0;
	virtual qf::ContainsResult contains_hash(const hashed_type& hashed) = 0;
	virtual qf::ContainsResult unsafe_contains(const key_type& key) const = 0;
	virtual qf::ContainsResult unsafe_contains_hash(
        const hashed_type& hashed) const = 0;
	virtual bool check_consistency() const = 0;
	virtual templated_qfilter_conc_base* create_bigger_QF(
        const hash_function_type& hf = hash_function_type()) const = 0;
	virtual void grow(templated_qfilter_conc_base<Key, Hash>& target_base,
                      size_t block_size) = 0;

	size_t capacity() const
	{ return ONE << quotient_bits; }

	void set_table_status(qf::TableStatus status)
	{ table_status = status; }
};


template <class Key, size_t remainder_bits, class Hash = htm::default_hash>
class templated_qfilter_conc : public templated_qfilter_conc_base<Key, Hash>
{
public:
    template <size_t new_remainder>
    using other_remainder = templated_qfilter_conc<Key, new_remainder, Hash>;
    static constexpr bool is_growing_compatible = true;
    static constexpr bool is_templated  = true;
    static constexpr bool is_sequential = false;
    static constexpr bool is_growing    = false;
    static constexpr bool is_dynamic    = false;
    static constexpr bool uses_handle   = false;

	using key_type = Key;
	using hash_function_type = Hash;
	using hashed_type        = decltype(std::declval<Hash>()(
                                            std::declval<Key>()));
	using next_smaller_qf    = templated_qfilter_conc<Key,
                                                      remainder_bits + 1,
                                                      Hash>;
	using next_bigger_qf     = templated_qfilter_conc<Key,
                                                      remainder_bits - 1,
                                                      Hash>;
	using base_qf = templated_qfilter_conc_base<Key, Hash>;

	friend base_qf;
	friend next_smaller_qf;

private:
	using cell_base_type = qf::grouped_cell_base_type;
	using cell_type      = qf::templated_cell_non_atomic<cell_base_type,
                                                         remainder_bits>;
	using cell_atomic    = typename cell_type::atomic;

	using entry_type = typename cell_type::entry_type;
	using entry_pos = qf::entry_position<cell_type::capacity>;

	using read_lock_guard = qf::templated_read_lock_guard<cell_atomic>;
	using write_lock_guard = qf::templated_write_lock_guard<cell_atomic>;

	ContainsResultEnum try_contains_lock_elision(const cell_type cell,
                                               const entry_pos pos,
                                               const entry_type remainder) const
	{
		// ensure no locks exist to the right
		size_t iter = pos.cell;
		do
		{
			if (cell.is_locked(iter))
				return ContainsResultEnum::FAILED_OPERATION;

			iter++;
		} while (iter < cell_type::capacity && !cell.is_empty(iter));

		// find cluster start && ensure no locks exist to the left
		size_t cluster_start = pos.cell;

		while (!cell.is_cluster_start(cluster_start))
		{
			if (cluster_start == 0 || cell.is_read_locked(cluster_start))
				return ContainsResultEnum::FAILED_OPERATION;
			--cluster_start;
		}

		// find run start
		size_t bucket = cluster_start;
		size_t run = bucket;

		while (bucket != pos.cell)
		{
			do
			{
				++run;
				if (run == cell_type::capacity)
					return ContainsResultEnum::FAILED_OPERATION;

			} while (cell.is_continuation(run));

			do
			{
				++bucket;
				if (bucket == cell_type::capacity)
					return ContainsResultEnum::FAILED_OPERATION;

			} while (!cell.is_occupied(bucket));
		}

		// search run
		do
		{
			if (cell.remainder(run) == remainder)
				return ContainsResultEnum::ENTRY_PRESENT;

			if (cell.remainder(run) > remainder)
				return ContainsResultEnum::ENTRY_NOT_PRESENT;

			++run;
			if (run == cell_type::capacity)
				return ContainsResultEnum::FAILED_OPERATION;

		} while (cell.is_continuation(run));

		return ContainsResultEnum::ENTRY_NOT_PRESENT;
	}

	qf::InsertResult try_insert_lock_elision(cell_type cell,
                                             const entry_pos pos,
                                             const entry_type remainder) const
	{
		auto original_cell = cell;
		if (cell.is_locked(pos.cell))
			return qf::FAILED_OPERATION;

		// find empty entry && ensure no locks exist to the right
		size_t empty_pos = pos.cell;

		do
		{
			++empty_pos;
			if (empty_pos == cell_type::capacity || cell.is_locked(empty_pos))
				return qf::FAILED_OPERATION;
		} while (!cell.is_empty(empty_pos));

		// find cluster start && ensure no locks exist to the left
		size_t cluster_start = pos.cell;

		while (!cell.is_cluster_start(cluster_start))
		{
			if (cluster_start == 0 || cell.is_read_locked(cluster_start))
				return qf::FAILED_OPERATION;
			--cluster_start;
		}

		const bool run_already_existing = cell.is_occupied(pos.cell);
		if (!run_already_existing)
			cell.set_occupied(true, pos.cell);

		// find run start
		size_t bucket = cluster_start;
		size_t run = bucket;

		while (bucket != pos.cell)
		{
			do
			{
				++run;
			} while (cell.is_continuation(run));

			do
			{
				++bucket;
			} while (!cell.is_occupied(bucket));
		}

		// find insert position in run
		size_t insert_pos = run;

		entry_type new_entry = 0;
		cell_type::entry_set_remainder(new_entry, remainder);

		if (run_already_existing)
		{
			do
			{
				if (cell.remainder(insert_pos) == remainder)
					return insert_pos - cluster_start; //element already present

				if (cell.remainder(insert_pos) > remainder)
					break;

				++insert_pos;
				if (insert_pos == cell_type::capacity)
					return qf::FAILED_OPERATION;

			} while (cell.is_continuation(insert_pos));

			if (insert_pos == run)
			{
				auto tmp_remainder = cell.remainder(run);
				cell.set_remainder(remainder, run);
				cell_type::entry_set_remainder(new_entry, tmp_remainder);

				++insert_pos;
				if (insert_pos == cell_type::capacity)
					return qf::FAILED_OPERATION;
			}

			cell_type::entry_set_continuation(new_entry, true);
		}

		cell_type::entry_set_shifted(new_entry, insert_pos != pos.cell);

		// insert and shift
		cell_type shifted_cell;

		// write unshifted entries at beginning of cell
		for (size_t i = 0; i < insert_pos; i++)
		{
			shifted_cell.set_entry(cell.entry(i), i);
		}

		// write new entry
		shifted_cell.set_entry(new_entry, insert_pos);

		// write shifted entries
		for (size_t i = insert_pos; i < empty_pos; i++)
		{
			shifted_cell.set_entry(cell.entry(i), i + 1);
			shifted_cell.set_shifted(true, i + 1);
		}

		// write unshifted entries at end of cell
		for (size_t i = empty_pos + 1; i < cell_type::capacity; i++)
		{
			shifted_cell.set_entry(cell.entry(i), i);
		}

		shifted_cell.set_occupied_full(cell.occupied_mask_full());

		if (!table[pos.table].CAS(original_cell, shifted_cell))
			return qf::FAILED_OPERATION;

		return insert_pos - cluster_start;
	}


	void increment_and_update(entry_pos& pos,
                              cell_type& cell,
                              entry_type& entry) const
	{
		++pos;
		if (pos.cell == 0)
		{
			cell = table[pos.table];
		}
		entry = cell.entry(pos.cell);
	}

	void decrement_and_update(entry_pos& pos,
                              cell_type& cell,
                              entry_type& entry) const
	{
		--pos;
		if (pos.cell == cell_type::capacity - 1)
		{
			cell = table[pos.table];
		}
		entry = cell.entry(pos.cell);
	}

	void force_update(const entry_pos& pos,
                      cell_type& cell,
                      entry_type& entry) const
	{
		cell = table[pos.table];
		entry = cell.entry(pos.cell);
	}

	bool try_empty_slot(cell_type local_cell,
                        entry_pos pos,
                        entry_type remainder)
	{
		cell_type filled_cell;

		while (true)
		{
			while (local_cell.is_write_locked(pos.cell))
			{
				if (local_cell.remainder(pos.cell) == 1)
					return false; // write lock set by growing is final

				local_cell = table[pos.table];
			}

			if (!local_cell.is_empty(pos.cell))
				return false;

			filled_cell = local_cell;
			filled_cell.set_remainder(remainder, pos.cell);
			filled_cell.set_occupied(true, pos.cell);

			if (table[pos.table].CAS(local_cell, filled_cell))
				return true;

			local_cell = table[pos.table];
		}
	}

	// does not use or check for locks
	bool unsafe_try_empty_slot(entry_pos pos, entry_type remainder)
	{
		cell_type local_cell = table[pos.table];
		cell_type filled_cell;

		while (local_cell.is_empty(pos.cell))
		{
			filled_cell = local_cell;
			filled_cell.set_remainder(remainder, pos.cell);
			filled_cell.set_occupied(true, pos.cell);

			if (table[pos.table].CAS(local_cell, filled_cell))
				return true;

			local_cell = table[pos.table];
		}

		return false;
	}

	// assumes that pos is not empty
	bool has_space_and_lock(entry_pos pos, write_lock_guard& write_lock)
	{
		++pos;

		cell_type local_cell = table[pos.table];
		entry_type local_entry = local_cell.entry(pos.cell);
		cell_type locked_cell;

		while (pos.table != table_capacity - 1)
		{
			while (cell_type::entry_is_write_locked(local_entry))
			{
				if (local_cell.remainder(pos.cell) == 1)
					return false; // write lock set by growing is final

				force_update(pos, local_cell, local_entry);
			}

			if (local_cell.is_empty(pos.cell))
			{
				// entry was empty and not locked -> acquire write lock
				locked_cell = local_cell;
				locked_cell.set_write_lock(pos.cell);

				if (table[pos.table].CAS(local_cell, locked_cell))
				{
					write_lock.set_data(&table[pos.table],
                                        pos,
                                        local_cell.status(pos.cell));

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
	void get_growing_write_lock(entry_pos pos, write_lock_guard& write_lock)
	{
		cell_type local_cell;
		cell_type locked_cell;
		entry_type locked_entry = 0;
		cell_type::entry_set_status(locked_entry, cell_type::write_lock_status);
		cell_type::entry_set_remainder(locked_entry, 1u);

		++pos;

		while (true)
		{
			do
			{
				local_cell = table[pos.table];
			} while (local_cell.is_write_locked(pos.cell));

			if (local_cell.is_empty(pos.cell))
			{
				locked_cell = local_cell;
				locked_cell.set_entry(locked_entry, pos.cell);

				if (table[pos.table].CAS(local_cell, locked_cell))
				{
					write_lock.set_data(&table[pos.table],
                                        pos,
                                        cell_type::write_lock_status);
					return;
				}
			}
			else
			{
				++pos;
			}
		}
	}

	entry_pos find_and_lock_cluster_start(entry_pos pos,
                                          read_lock_guard& read_lock)
	{
		cell_type local_cell = table[pos.table];
		entry_type local_entry = local_cell.entry(pos.cell);
		cell_type locked_cell;

		while (true)
		{
			while (cell_type::entry_is_read_locked(local_entry))
			{
				force_update(pos, local_cell, local_entry);
			}

			if (cell_type::entry_is_cluster_start(local_entry))
			{
				// entry was cluster start and not locked -> acquire read lock
				locked_cell = local_cell;
				locked_cell.set_read_lock(pos.cell);

				if (table[pos.table].CAS(local_cell, locked_cell))
				{
					read_lock.set_data(&table[pos.table],
                                       pos,
                                       local_cell.status(pos.cell));
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

	entry_pos find_cluster_start(entry_pos pos) const
	{
		cell_type local_cell = table[pos.table];

		while (!local_cell.is_cluster_start(pos.cell))
		{
			--pos;
			local_cell = table[pos.table];
		}

		return pos;
	}

	entry_pos find_run_start(entry_pos cluster_start,
                             entry_pos pos,
                             bool set_occupied = false) const
	{
		entry_pos bucket = cluster_start;
		entry_pos run = bucket;

		if (set_occupied)
			table[pos.table].set_occupied(true, pos.cell);

		entry_type local_entry;

		while (bucket != pos)
		{
			do
			{
				++run;

				do
				{
					local_entry = table[run.table].entry(run.cell);
				} while (cell_type::entry_is_read_locked(local_entry));
			} while (!cell_type::entry_is_write_locked(local_entry)
                     && cell_type::entry_is_continuation(local_entry));

			do
			{
				++bucket;
				local_entry = table[bucket.table].entry(bucket.cell);
			} while (!cell_type::entry_is_occupied(local_entry));
		}

		return run;
	}

	// does not use or check for locks
	entry_pos unsafe_find_run_start(entry_pos cluster_start,entry_pos pos) const
	{
		entry_pos bucket = cluster_start;
		entry_pos run = bucket;

		entry_type local_entry;

		while (bucket != pos)
		{
			do
			{
				++run;
				local_entry = table[run.table].entry(run.cell);
			} while (cell_type::entry_is_continuation(local_entry));

			do
			{
				++bucket;
				local_entry = table[bucket.table].entry(bucket.cell);
			} while (!cell_type::entry_is_occupied(local_entry));
		}

		return run;
	}

	bool acquire_read_lock(entry_pos pos, read_lock_guard& read_lock)
	{
		if (read_lock.locked() && read_lock.pos == pos)
			return true;

		cell_type local_cell;
		cell_type locked_cell;

		while (true)
		{
			do
			{
				local_cell = table[pos.table];
			} while (local_cell.is_read_locked(pos.cell));

			if (!local_cell.is_cluster_start(pos.cell))
			{
				return false;
			}

			locked_cell = local_cell;
			locked_cell.set_read_lock(pos.cell);

			if (table[pos.table].CAS(local_cell, locked_cell))
			{
				read_lock.set_data(&table[pos.table],
                                   pos,
                                   local_cell.status(pos.cell));
				return true;
			}
		}
	}

	// insert

	void insert_and_shift(entry_type entry,
                          entry_pos& insert_pos,
                          write_lock_guard& write_lock,
                          read_lock_guard&  read_lock)
	{
		if (insert_pos == read_lock.pos)
		{
			// note: read lock is not on the run start of entry's canonical run

			const entry_type remainder = table[insert_pos.table].remainder(
                insert_pos.cell);
			table[insert_pos.table].set_remainder(
                cell_type::entry_remainder(entry),
                insert_pos.cell);

			const auto status = read_lock.status;
			read_lock.status = cell_type::entry_status(entry);
			cell_type::status_set_occupied(read_lock.status,
                                          cell_type::entry_is_occupied(status));

			cell_type::entry_set_remainder(entry, remainder);
			cell_type::entry_set_status(entry, cell_type::entry_status(status));
			cell_type::entry_set_shifted(entry, true);
			cell_type::entry_set_occupied(entry, false);

			++insert_pos;
		}

		while (insert_pos.table != write_lock.pos.table)
		{
			insert_entry_full(entry, insert_pos);
			insert_pos.table++;
			insert_pos.cell = 0;
		};

		insert_entry_not_full(entry, insert_pos, write_lock.pos.cell);
		write_lock.clear();
	}

	// inserts entry at pos and shifts the following entries
	// cell at pos must have enough space for the shift i.e. an empty entry
    // (write lock)
	void insert_entry_not_full(entry_type entry,
                               entry_pos insert_pos,
                               size_t write_lock_pos)
	{
		cell_type cell;
		cell_type shifted_cell;

		const size_t read_lock_check_start = insert_pos.cell;

		do
		{
			shifted_cell.clear();

			do
			{
				cell = table[insert_pos.table];
			} while (cell.is_read_locked(read_lock_check_start,write_lock_pos));

			// write unshifted entries at beginning of cell
			for (size_t i = 0; i < insert_pos.cell; i++)
			{
				shifted_cell.set_entry(cell.entry(i), i);
			}

			// write new entry
			shifted_cell.set_entry(entry, insert_pos.cell);

			// write shifted entries
			for (size_t i = insert_pos.cell; i < write_lock_pos; i++)
			{
				shifted_cell.set_entry(cell.entry(i), i + 1);
				shifted_cell.set_shifted(true, i + 1);
			}

			// write unshifted entries at end of cell
			for (size_t i = write_lock_pos + 1; i < cell_type::capacity; i++)
			{
				shifted_cell.set_entry(cell.entry(i), i);
			}

			shifted_cell.set_occupied_full(cell.occupied_mask_full());
		} while (!table[insert_pos.table].CAS(cell, shifted_cell));
	}

	// inserts entry at pos and shifts the following entries
	// the last entry which is lost due to the shift is written back to
    // the entry reference
	void insert_entry_full(entry_type& entry, entry_pos insert_pos)
	{
		cell_type cell;
		cell_type shifted_cell;

		const size_t read_lock_check_start = insert_pos.cell;

		do
		{
			shifted_cell.clear();

			do
			{
				cell = table[insert_pos.table];
			} while (cell.is_read_locked(read_lock_check_start,
                                         cell_type::capacity));

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
		} while (!table[insert_pos.table].CAS(cell, shifted_cell));

		entry = cell.entry(cell_type::capacity - 1);
		cell_type::entry_set_shifted(entry, true);
	}

	entry_type wait_on_insert_write_lock(entry_pos pos) const
	{
		cell_type local_cell;

		do
		{
			local_cell = table[pos.table];
		} while (local_cell.is_write_locked(pos.cell)
                 && local_cell.remainder(pos.cell) == 0);

		return local_cell.entry(pos.cell);
	}

	void grow(next_bigger_qf& target, size_t block_start, size_t block_size)
	{
		read_lock_guard read_lock;
		write_lock_guard write_lock;
		entry_pos cluster_start{ block_start, 0 };
		const entry_pos block_end{ std::min(block_start + block_size,
                                            table_capacity), 0 };
		entry_type local_entry = table[cluster_start.table].entry(
            cluster_start.cell);

		// skip write lock
		if (cell_type::entry_is_write_locked(local_entry))
		{
			wait_on_insert_write_lock(cluster_start);
			++cluster_start;
			local_entry = table[cluster_start.table].entry(cluster_start.cell);
		}

		// skip partial cluster
		if (!cell_type::entry_is_empty(local_entry)
            && cluster_start != entry_pos{0, 0})
		{
			auto prev_pos = cluster_start;
			--prev_pos;
			const entry_type prev_entry = wait_on_insert_write_lock(prev_pos);

			if (!cell_type::entry_is_empty(prev_entry)
                && !cell_type::entry_is_write_locked(prev_entry))
			{
				// cluster start is in "middle" of cluster -> skip cluster
				do
				{
					local_entry = wait_on_insert_write_lock(++cluster_start);

					if (cell_type::entry_is_write_locked(local_entry))
					{
						// can only be a write lock from growing
						++cluster_start;
						break;
					}
				} while (cluster_start < block_end
                         && !cell_type::entry_is_empty(local_entry));
			}
		}

		// skip empty entries
		while (cluster_start < block_end
               && table[cluster_start.table].is_empty(cluster_start.cell))
		{
			++cluster_start;
		}

		if (cluster_start >= block_end)
			return;

		do
		{
			get_growing_write_lock(cluster_start, write_lock);
			acquire_read_lock(cluster_start, read_lock);
			cluster_start = grow_cluster(target, read_lock);
		} while (cluster_start < block_end);
	}

	entry_pos grow_cluster(next_bigger_qf& target,
                           const read_lock_guard& read_lock) const
	{
		static constexpr size_t remainder_bits_trgt = remainder_bits - 1;
		using cell_type_trgt = typename next_bigger_qf::cell_type;
		using entry_pos_trgt = typename next_bigger_qf::entry_pos;
		using entry_type_trgt = typename next_bigger_qf::entry_type;

		static constexpr entry_type_trgt target_remainder_mask =
            ~((~0ull) << remainder_bits_trgt);
		auto run_start = read_lock.pos;
		auto run_iter = read_lock.pos;

		entry_type source_entry, tmp_entry;
		entry_type_trgt target_entry = 0;
		bool leading_bit, last_leading_bit;
		entry_pos_trgt target_pos, last_target_pos;

		auto read = [&]()
		{
			do
			{
				source_entry = table[run_iter.table].entry(run_iter.cell);
			} while (run_iter != read_lock.pos
                     && cell_type::entry_is_read_locked(source_entry));

			if (cell_type::entry_is_read_locked(source_entry))
				cell_type::entry_set_status(source_entry, read_lock.status);
			else if (cell_type::entry_is_write_locked(source_entry))
				cell_type::entry_set_status(source_entry, 0);

			leading_bit = cell_type::entry_remainder(source_entry)
                          >> remainder_bits_trgt;
			cell_type_trgt::entry_set_remainder(target_entry,
                                 cell_type::entry_remainder(source_entry)
                                  & target_remainder_mask);
		};

		auto advance = [&]()
		{
			// write
			target.table[target_pos.table].set_entry(target_entry,
                                                     target_pos.cell);

			// backup state
			last_leading_bit = leading_bit;
			last_target_pos = target_pos;

			target_pos++;
			run_iter++;

			read();
		};

		read();
		target_pos = entry_pos_trgt(run_iter.offset() * 2
                                    + static_cast<size_t>(leading_bit));

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

				while (cell_type::entry_is_continuation(source_entry)
                       && leading_bit == last_leading_bit)
				{
					advance();
				}

				if (!second_pass
                    && cell_type::entry_is_continuation(source_entry))
				{
					// leading bits differ -> new run_iter start
					second_pass = true;
					target_pos = entry_pos_trgt(run_start.offset() * 2
                                            + static_cast<size_t>(leading_bit));
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
				target_pos = entry_pos_trgt(run_start.offset() * 2
                                       + static_cast<size_t>(last_leading_bit));
				target.table[target_pos.table].set_occupied(true,
                                                            target_pos.cell);
			}

			second_pass = false;

			++run_start;
			for (; run_start.table < table_capacity; ++run_start)
			{
				do
				{
					tmp_entry = table[run_start.table].entry(run_start.cell);
				} while (run_start != read_lock.pos
                         && cell_type::entry_is_read_locked(tmp_entry));

				if (cell_type::entry_is_read_locked(tmp_entry))
					cell_type::entry_set_status(tmp_entry, read_lock.status);

				if (!cell_type::entry_is_write_locked(tmp_entry)
                    && cell_type::entry_is_occupied(tmp_entry))
					break;
			}

			target_pos = entry_pos_trgt(run_start.offset() * 2
                                        + static_cast<size_t>(leading_bit));
		}

		return run_start;
	}

	alignas(128) size_t				                table_capacity;
	alignas(128) std::unique_ptr<cell_atomic[]>     table;
	alignas(128) hash_function_type                             hf;
	alignas(128) size_t                             fingerprint_mask;
	alignas(128) std::atomic_size_t                 growing_block_counter  = 0;
	alignas(128) std::atomic_size_t                 finished_block_counter = 0;

public:

	templated_qfilter_conc(size_t capacity = qf::DEFAULT_MIN_CAPACITY,
                           const hash_function_type& hf = hash_function_type())
	 : base_qf(capacity), hf(hf)
	{
		table_capacity = (ONE << this->quotient_bits) / cell_type::capacity
                         + qf::shift_buffer;
		table = std::make_unique<cell_atomic[]>(table_capacity);
		fingerprint_mask = (ONE << static_cast<size_t>(remainder_bits
                                                       +this->quotient_bits))-1;
	}

	templated_qfilter_conc(const templated_qfilter_conc& other)
		: base_qf(other),
          table_capacity(other.table_capacity),
          table(new cell_atomic[table_capacity]),
          hf(other.hf),
          fingerprint_mask(other.fingerprint_mask),
          growing_block_counter(other.growing_block_counter),
          finished_block_counter(other.finished_block_counter)
	{
		for (size_t i = 0; i < table_capacity; i++)
		{
			table[i] = other.table[i];
		}
	}

	templated_qfilter_conc& operator=(const templated_qfilter_conc& other)
	{
		base_qf::operator=(other);
		table_capacity = other.table_capacity;
		table = std::make_unique<cell_atomic[]>(other.table_capacity);
		hf = other.hf;
		fingerprint_mask = other.fingerprint_mask;
		growing_block_counter = other.growing_block_counter.load();
		finished_block_counter = other.finished_block_counter.load();

		for (size_t i = 0; i < table_capacity; i++)
		{
			table[i] = other.table[i];
		}

		return *this;
	}

	templated_qfilter_conc(templated_qfilter_conc&& other) = default;
	templated_qfilter_conc& operator=(templated_qfilter_conc&& other) = default;

	~templated_qfilter_conc() = default;

	base_qf* create_bigger_QF(const hash_function_type& hash
                              = hash_function_type()) const override
	{
		if constexpr (remainder_bits > 1)
		{
			return new next_bigger_qf(ONE << (this->quotient_bits + 1), hash);
		}
		else
		{
			return nullptr;
		}
	}

    void prefetch(const hashed_type& hashed) const override
    {
        const auto[q, r] = this->get_quotient_and_remainder(hashed);
		const entry_pos q_pos(q);
        __builtin_prefetch(&table[q_pos.first]);
    }

	qf::InsertResult insert(const key_type& key) override
	{
		return insert_hash(hf(key));
	}

	qf::InsertResult insert_hash(const hashed_type& hashed) override
	{
		if (this->table_status != qf::TableStatus::active)
			return qf::FAILED_OPERATION;

		const auto[q, r] = this->get_quotient_and_remainder(hashed);
		const entry_pos q_pos(q);
		const cell_type q_cell = table[q_pos.table];

		if (try_empty_slot(q_cell, q_pos, r))
			return 1;

		if constexpr (cell_type::capacity > 1)
		{
			auto lock_elision_res = try_insert_lock_elision(q_cell, q_pos, r);
			if (lock_elision_res.successful())
				return lock_elision_res;
		}

		write_lock_guard write_lock;
		if (!has_space_and_lock(q_pos, write_lock))
			return qf::FAILED_OPERATION;

		// works even if q_pos is read_locked since it was occupied before
        // locking and read_lock mask set the occupied bit
		const bool run_already_existing = table[q_pos.table].is_occupied(q_pos.cell);

		read_lock_guard read_lock;
		const auto cluster_start = find_and_lock_cluster_start(q_pos,
                                                               read_lock);
		const auto run = find_run_start(cluster_start,
                                        q_pos,
                                        !run_already_existing);
		auto insert_pos = run;

		entry_type new_entry = 0;
		cell_type::entry_set_remainder(new_entry, r);

		if (insert_pos == write_lock.pos)
		{
			cell_type::entry_set_shifted(new_entry, true);
			table[write_lock.pos.table].set_entry(new_entry,
                                                  write_lock.pos.cell);
			write_lock.clear();
			return entry_pos::distance(cluster_start, insert_pos);
		}
		else if (run.table == table_capacity - 1)
		{
			if (!run_already_existing)
				table[q_pos.table].set_occupied(run_already_existing,
                                                q_pos.cell);

			return qf::FAILED_OPERATION;
		}

		cell_type current_cell = table[insert_pos.table];
		entry_type current_entry = current_cell.entry(insert_pos.cell);

		if (run_already_existing)
		{
			// find position in existing run
			do
			{
				if (cell_type::entry_remainder(current_entry) == r)
                    // element already present
					return entry_pos::distance(cluster_start, insert_pos);
				else if (cell_type::entry_remainder(current_entry) > r)
					break;

				increment_and_update(insert_pos, current_cell, current_entry);

				while (cell_type::entry_is_read_locked(current_entry))
				{
					force_update(insert_pos, current_cell, current_entry);
				}
			} while (!cell_type::entry_is_locked(current_entry)
                     && cell_type::entry_is_continuation(current_entry));

			if (insert_pos == run)
			{
				// note: possibly writing to read locked entry but only
                //       modifying the remainder
				table[run.table].set_remainder(r, run.cell);
				cell_type::entry_set_remainder(new_entry,
                                     cell_type::entry_remainder(current_entry));

				increment_and_update(insert_pos, current_cell, current_entry);

				while (cell_type::entry_is_read_locked(current_entry))
				{
					force_update(insert_pos, current_cell, current_entry);
				}
			}

			cell_type::entry_set_continuation(new_entry, true);
		}

		//TODO: remove this?
		if (cell_type::entry_is_read_locked(current_entry)
            || cell_type::entry_is_cluster_start(current_entry))
			acquire_read_lock(insert_pos, read_lock);

		cell_type::entry_set_shifted(new_entry, insert_pos != q_pos);

		insert_and_shift(new_entry, insert_pos, write_lock, read_lock);
		return entry_pos::distance(cluster_start, insert_pos);
	}

	bool quick_insert(const key_type& key) override
	{
		return quick_insert_hash(hf(key));
	}

	// try insert without shifting (does not use or check for locks)
	bool quick_insert_hash(const hashed_type& hashed) override
	{
		const auto[q, r] = this->get_quotient_and_remainder(hashed);
		const entry_pos q_pos(q);

		if (unsafe_try_empty_slot(q_pos, r))
			return 1;

		return unsafe_contains_hash(hashed);
	}

	qf::ContainsResult contains(const key_type& key) override
	{
		return contains_hash(hf(key));
	}

	qf::ContainsResult contains_hash(const hashed_type& hashed) override
	{
		if (this->table_status != qf::TableStatus::active)
			return qf::FAILED_OPERATION;

		const auto[q, r] = this->get_quotient_and_remainder(hashed);
		const entry_pos q_pos(q);
		const cell_type q_cell = table[q_pos.table];
		const entry_type q_entry = q_cell.entry(q_pos.cell);

		if (cell_type::entry_is_empty(q_entry))
			return ContainsResultEnum::CANONICAL_SLOT_EMPTY;
		if (!cell_type::entry_is_occupied(q_entry))
			return ContainsResultEnum::ENTRY_NOT_PRESENT;

		if constexpr (cell_type::capacity > 1)
		{
			auto lock_elision_result = try_contains_lock_elision(q_cell,
                                                                 q_pos,
                                                                 r);
			if (lock_elision_result != ContainsResultEnum::FAILED_OPERATION)
				return lock_elision_result;
		}

		read_lock_guard read_lock;
		auto iter = find_run_start(find_and_lock_cluster_start(q_pos,
                                                               read_lock),
                                   q_pos);

		cell_type current_cell = table[iter.table];
		entry_type current_entry = current_cell.entry(iter.cell);

		do
		{
			if (cell_type::entry_remainder(current_entry) == r)
				return ContainsResultEnum::ENTRY_PRESENT;

			// early abort due to sort is acutally slower in + and - contains
			//if (cell_type::entry_remainder(current_entry) > r)
			//	return ContainsResultEnum::ENTRY_NOT_PRESENT;

			increment_and_update(iter, current_cell, current_entry);
		} while (!cell_type::entry_is_locked(current_entry)
                 && cell_type::entry_is_continuation(current_entry));

		return ContainsResultEnum::ENTRY_NOT_PRESENT;
	}

	qf::ContainsResult unsafe_contains(const key_type& key) const override
	{
		return unsafe_contains_hash(hf(key));
	}
	// does not use or check for locks
	qf::ContainsResult unsafe_contains_hash(
        const hashed_type& hashed) const override
	{
		const auto[q, r] = this->get_quotient_and_remainder(hashed);
		const entry_pos q_pos(q);
		const entry_type q_entry = table[q_pos.table].entry(q_pos.cell);

		if (cell_type::entry_is_empty(q_entry))
			return ContainsResultEnum::CANONICAL_SLOT_EMPTY;
		if (!cell_type::entry_is_occupied(q_entry))
			return ContainsResultEnum::ENTRY_NOT_PRESENT;

		auto iter = unsafe_find_run_start(find_cluster_start(q_pos), q_pos);

		entry_type current_entry = table[iter.table].entry(iter.cell);

		do
		{
			if (cell_type::entry_remainder(current_entry) == r)
				return ContainsResultEnum::ENTRY_PRESENT;

			++iter;
			current_entry = table[iter.table].entry(iter.cell);
		} while (cell_type::entry_is_continuation(current_entry));

		return ContainsResultEnum::ENTRY_NOT_PRESENT;
	}

	std::pair<size_t, entry_type> get_quotient_and_remainder(
        const hashed_type& hashed) const
	{
		return qf::get_quotient_and_remainder<entry_type, remainder_bits>
            (hashed, this->quotient_bits);
	}

	virtual void grow(base_qf& target_base, size_t block_size) override
	{
		if constexpr (remainder_bits <= 1)
		{
			return;
		}
		else
		{
			if (this->table_status != qf::TableStatus::growing)
				return;

			auto& target = dynamic_cast<next_bigger_qf&>(target_base);

			const size_t total_blocks = (table_capacity + block_size - 1)
                                        / block_size;
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

	void print() const
	{
		print(0, table_capacity);
	}

	void print(entry_pos start, entry_pos end) const
	{
		print(start.table, end.table - start.table + 1);
	}

	void print(size_t start, size_t count) const
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

	void print_cell(size_t table_pos) const
	{
		print_cell(table[table_pos], table_pos);
	}

	void print_cell(cell_type cell, size_t table_pos) const
	{
		for (size_t cell_pos = 0; cell_pos < cell_type::capacity; cell_pos++)
		{
			print_entry(cell.entry(cell_pos), { table_pos, cell_pos });
			std::cout << "   ";
		}
	}

	void print_entry(entry_pos pos) const
	{
		print_entry(table[pos.table].entry(pos.cell), pos.offset());
	}

	void print_entry(entry_type entry, entry_pos pos) const
	{
		print_entry(entry, pos.offset());
	}

	void print_entry(entry_type entry, size_t index) const
	{
		std::cout << cell_type::entry_is_occupied(entry) << "|"
                  << cell_type::entry_is_continuation(entry) << "|"
                  << cell_type::entry_is_shifted(entry) << " || "
                  << qf::bitstring(index, this->quotient_bits) << "|"
                  << qf::bitstring(cell_type::entry_remainder(entry),
                                   remainder_bits);
	}

	void print_entry(entry_type entry) const
	{
		std::cout << cell_type::entry_is_occupied(entry) << "|"
                  << cell_type::entry_is_continuation(entry) << "|"
                  << cell_type::entry_is_shifted(entry) << " ||     "
                  << qf::bitstring(cell_type::entry_remainder(entry),
                                   remainder_bits);
	}

	bool check_consistency() const override
	{
		bool consistent = true;
		bool sorted = true;
		entry_type last = 0;
		entry_type current;

		for (entry_pos pos{ 0,0 }; pos.table < table_capacity; ++pos)
		{
			current = table[pos.table].entry(pos.cell);

			/*if (cell_type::entry_is_read_locked(current))
			{
				print_entry(current);
				std::cerr << "-> error: " << pos << " read locked" << std::endl;
				consistent = false;
			}
			else if (cell_type::entry_is_write_locked(current))
			{
				print_entry(current);
				std::cerr << "-> error: " << pos << " write locked"<< std::endl;
				consistent = false;
			}
			else */if (cell_type::entry_is_continuation(current)
                       && !cell_type::entry_is_shifted(current))
			{
				print_entry(current);
				std::cerr << "-> error: " << pos << " cont =/=> shift"
                          << std::endl;
				consistent = false;
			}
			else if (cell_type::entry_is_continuation(current)
                     && !cell_type::entry_is_shifted(last)
                     && !cell_type::entry_is_cluster_start(last))
			{
				print_entry(current);
				std::cerr << "-> error2: " << pos << " cont =/=> shift"
                          << std::endl;
				consistent = false;
			}
			else if (cell_type::entry_is_continuation(current))
			{
				if (cell_type::entry_remainder(current)
                    < cell_type::entry_remainder(last))
				{
					sorted = false;
					consistent = false;
				}
			}

			last = current;
		}

		if (!sorted)
			std::cerr << "error:  remainder not sorted" << std::endl;

		return consistent;
	}
};

} // namespace qf
