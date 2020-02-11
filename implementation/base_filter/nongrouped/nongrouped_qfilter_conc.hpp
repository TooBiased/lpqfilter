#pragma once
#include <iostream>
#include <memory>
#include "utils/utilities.h"
#include "nongrouped_qfilter_cell.hpp"


namespace qf
{

template < class Key,
	class Hash = default_hash>
class nongrouped_qfilter_conc
{
public:
	using key_type = Key;
	using hasher = Hash;
	using HashVal = decltype(std::declval<Hash>()(std::declval<Key>()));

private:
	using entry_type = std::uint64_t;

	using cell_atomic = qf::nongrouped_cell_atomic<entry_type>;
	using cell = cell_atomic::non_atomic;
	using cell_base = cell_atomic::base;

	using read_lock_guard  = qf::nongrouped_read_lock_guard<entry_type>;
	using write_lock_guard = qf::nongrouped_write_lock_guard<entry_type>;
	using lock_guard       = qf::nongrouped_cell_lock_guard<entry_type>;


	bool try_empty_slot(size_t pos, entry_type remainder)
	{
		cell local_entry;
		cell filled_entry;
		filled_entry.set_remainder(remainder);
		filled_entry.set_occupied(true);

		while (true)
		{
			for (local_entry = table[pos];
                 local_entry.is_write_locked();
                 local_entry = table[pos])
			{
				if (local_entry.remainder() == 1)
					return false; // write lock set by growing is final
			}

			if (!local_entry.is_empty())
				return false;

			if (table[pos].CAS(local_entry, filled_entry))
				return true;
		}
	}

	// does not use or check for locks
	bool unsafe_try_empty_slot(size_t pos, entry_type remainder)
	{
		cell local_entry;
		cell filled_entry;
		filled_entry.set_remainder(remainder);
		filled_entry.set_occupied(true);

		while (local_entry.is_empty())
		{
			if (table[pos].CAS(local_entry, filled_entry))
				return true;
		}

		return false;
	}

	// assumes that pos is not empty
	bool has_space_and_lock(size_t pos, write_lock_guard& write_lock)
	{
		cell local_entry;
		cell locked_entry;
		pos++;

		while (pos != table_capacity - 1)
		{
			for (local_entry = table[pos];
                 local_entry.is_write_locked();
                 local_entry = table[pos])
			{
				if (local_entry.remainder() == 1)
					return false; // write lock set by growing is final
			}

			if (local_entry.is_empty())
			{
				// entry was empty and not locked -> acquire writer lock
				locked_entry = local_entry;
				locked_entry.set_write_lock();

				if (table[pos].CAS(local_entry, locked_entry))
				{
					set_lock_data(write_lock, pos, local_entry.status());

					if (table_status != qf::TableStatus::active)
					{
						write_lock.release();
						return false;
					}

					return true;
				}
			}
			else
			{
				pos++;
			}
		}

		return false;
	}

	// assumes that pos is not empty
	void get_growing_write_lock(size_t pos, write_lock_guard& write_lock)
	{
		cell local_entry;
		cell locked_entry;
		locked_entry.set_write_lock();
		locked_entry.set_remainder(1u);

		++pos;

		while (true)
		{
			do
			{
				local_entry = table[pos];
			} while (local_entry.is_write_locked());

			if (!local_entry.is_empty())
			{
				++pos;
			}
			else if (table[pos].CAS(local_entry, locked_entry))
			{
				set_lock_data(write_lock, pos, cell_base::write_lock_status);
				return;
			}
		}
	}

	size_t find_and_lock_cluster_start(size_t pos, read_lock_guard& read_lock)
	{
		cell local_entry;
		cell locked_entry;

		while (true)
		{
			do
			{
				local_entry = table[pos];
			} while (local_entry.is_read_locked());

			if (local_entry.is_cluster_start())
			{
				// slot was cluster start and not locked -> acquire read lock
				locked_entry = local_entry;
				locked_entry.set_read_lock();

				if (table[pos].CAS(local_entry, locked_entry))
				{
					set_lock_data(read_lock, pos, local_entry.status());
					return pos;
				}
			}
			else
			{
				pos--;
			}
		}
	}

	size_t find_cluster_start(size_t pos) const
	{
		cell local_entry = table[pos];

		while (!local_entry.is_cluster_start())
		{
			local_entry = table[--pos];
		}

		return pos;
	}

	size_t find_run_start(size_t cluster_start,
                          size_t pos,
                          read_lock_guard& read_lock,
                          bool set_occupied = false)
	{
		size_t bucket = cluster_start;
		size_t run = bucket;

		if (set_occupied)
			set_occupied_atomic(pos, read_lock);

		cell local_entry;

		while (bucket != pos)
		{
			do
			{
				run++;

				do
				{
					local_entry = table[run];
				} while (local_entry.is_read_locked());
			} while (!local_entry.is_write_locked()
                     && local_entry.is_continuation());

			do
			{
				bucket++;
				local_entry = table[bucket];
			} while (!local_entry.is_occupied());
		}

		return run;
	}

	// does not use or check for locks
	size_t unsafe_find_run_start(size_t cluster_start, size_t pos) const
	{
		size_t bucket = cluster_start;
		size_t run = bucket;

		cell local_entry;

		while (bucket != pos)
		{
			do
			{
				local_entry = table[++run];
			} while (local_entry.is_continuation());

			do
			{
				local_entry = table[++bucket];
			} while (!local_entry.is_occupied());
		}

		return run;
	}

	// insert a remainder in a slot that is not its canonical slot
	void insert_and_shift(cell entry,
                          size_t& insert_pos,
                          write_lock_guard& write_lock,
                          read_lock_guard& read_lock)
	{
		if (insert_pos == read_lock.pos)
		{
			// note: read lock is not on the run start of entry's canonical run

			const entry_type remainder = table[insert_pos].remainder();
			table[insert_pos].set_remainder(entry.remainder());

			const auto old_read_lock_status = read_lock.status;
			read_lock.status = entry.status();
			cell::status_set_occupied(read_lock.status,
                                      cell::status_is_occupied(
                                          old_read_lock_status));

			entry.set_remainder(remainder);
			entry.set_status(old_read_lock_status);
			entry.set_shifted(true);
			entry.set_occupied(false);

			++insert_pos;
		}

		cell expected;
		cell current = entry;

		while (insert_pos != write_lock.pos)
		{
			do
			{
				do
				{
					expected = table[insert_pos];
				} while (expected.is_read_locked());

				if (expected.is_occupied())
				{
					current.set_occupied(true);
				}
			} while (!table[insert_pos].CAS(expected, current));

			current = expected;
			current.set_shifted(true);
			current.set_occupied(false);
			insert_pos++;
		}

		table[insert_pos] = current;
		write_lock.clear();
	}

	void set_occupied_atomic(size_t pos, lock_guard& lock)
	{
		if (pos == lock.pos)
		{
			lock.status &= ~cell::occupied_mask;
			lock.status |= cell::occupied_mask;
			return;
		}

		set_occupied_atomic(pos);
	}

	void set_occupied_atomic(size_t pos)
	{
		cell local_entry;
		cell new_entry;

		while (true)
		{
			do
			{
				local_entry = table[pos];
			} while (local_entry.is_locked());

			new_entry = local_entry;
			new_entry.set_occupied(true);

			if (table[pos].CAS(local_entry, new_entry))
				break;
		}
	}

	void set_lock_data(lock_guard& lock, size_t pos, qf::status_type status)
	{
		lock.set_data(pos, table[pos], status);
	}

	bool acquire_read_lock(size_t pos, read_lock_guard& read_lock)
	{
		if (read_lock.locked() && read_lock.pos == pos)
			return true;

		cell local_entry;
		cell locked_entry;

		while (true)
		{
			do
			{
				local_entry = table[pos];
			} while (local_entry.is_read_locked());

			if (!local_entry.is_cluster_start()) {
				return false;
			}

			locked_entry = local_entry;
			locked_entry.set_read_lock();

			if (table[pos].CAS(local_entry, locked_entry))
			{
				read_lock.set_data(pos, table[pos], local_entry.status());
				return true;
			}
		}
	}

	cell wait_on_insert_write_lock(size_t pos)
	{
		cell local_entry;

		do
		{
			local_entry = table[pos];
		} while (local_entry.is_write_locked() && local_entry.remainder() == 0);

		return local_entry;
	}

	void grow(nongrouped_qfilter_conc& target,
              size_t block_start,
              size_t block_size)
	{
		read_lock_guard read_lock;
		write_lock_guard write_lock;
		size_t cluster_start = block_start;
		const size_t block_end = std::min(block_start + block_size,
                                          table_capacity);
		cell local_entry = table[cluster_start];

		// skip write lock
		if (local_entry.is_write_locked())
		{
			wait_on_insert_write_lock(cluster_start);
			local_entry = table[++cluster_start];
		}

		// skip partial cluster
		if (!local_entry.is_empty() && cluster_start != 0)
		{
			cell prev_entry = wait_on_insert_write_lock(cluster_start - 1);

			if (!prev_entry.is_empty() && !prev_entry.is_write_locked())
			{
				// cluster start is in "middle" of cluster -> skip cluster
				do
				{
					local_entry = wait_on_insert_write_lock(++cluster_start);

					if (local_entry.is_write_locked())
					{
						// can only be a write lock from growing
						cluster_start++;
						break;
					}
				} while (cluster_start < block_end && !local_entry.is_empty());
			}
		}

		// skip empty entries
		while (cluster_start < block_end && table[cluster_start].is_empty())
		{
			cluster_start++;
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

	size_t grow_cluster(nongrouped_qfilter_conc& target,
                        read_lock_guard& read_lock)
	{
		const entry_type target_remainder_mask =
            ~((~static_cast<entry_type>(0)) << target.remainder_bits);
		size_t run_start = read_lock.pos;
		size_t run_iter = read_lock.pos;

		cell source_entry, target_entry, tmp_entry;
		bool leading_bit, last_leading_bit;
		size_t target_index, last_target_index = 0;

		auto read = [&]()
		{
			do
			{
				source_entry = table[run_iter];
			} while (run_iter != read_lock.pos
                     && source_entry.is_read_locked());

			if (source_entry.is_read_locked())
				source_entry.set_status(read_lock.status);
			else if (source_entry.is_write_locked())
				source_entry.set_status(0);

			leading_bit = source_entry.remainder() >> target.remainder_bits;
			target_entry.set_remainder(source_entry.remainder()
                                       & target_remainder_mask);
		};

		auto advance = [&]()
		{
			// write
			target.table[target_index] = target_entry;

			// backup state
			last_leading_bit = leading_bit;
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

				while (source_entry.is_continuation()
                       && leading_bit == last_leading_bit)
				{
					advance();
				}

				if (!second_pass && source_entry.is_continuation())
				{
					// leading bits differ -> new run_iter start
					second_pass = true;
					target_index = run_start * 2
                                   + static_cast<size_t>(leading_bit);
					continue;
				}
			}

			if (second_pass)
			{
				target.set_occupied_atomic(run_start * 2);
				target.set_occupied_atomic(run_start * 2 + 1);
			}
			else
			{
				target.set_occupied_atomic(run_start * 2
                                       + static_cast<size_t>(last_leading_bit));
			}

			second_pass = false;

			while (++run_start < table_capacity)
			{
				do
				{
					tmp_entry = table[run_start];
				} while (run_start != read_lock.pos
                         && tmp_entry.is_read_locked());

				if (tmp_entry.is_read_locked())
					tmp_entry.set_status(read_lock.status);

				if (!tmp_entry.is_write_locked() && tmp_entry.is_occupied())
					break;
			}

			target_index = run_start * 2 + static_cast<size_t>(leading_bit);
		}

		return run_start;
	}

	alignas(128) size_t                           table_capacity;
	alignas(128) std::unique_ptr<cell_atomic[]>   table;
	alignas(128) std::atomic<qf::TableStatus>     table_status{ qf::TableStatus::active };
	alignas(128) hasher                           hf;
	alignas(128) size_t                           fingerprint_mask;
	alignas(128) size_t                           remainder_bits;
	alignas(128) size_t                           quotient_bits;
	alignas(128) std::atomic_size_t               growing_block_counter  = 0;
	alignas(128) std::atomic_size_t               finished_block_counter = 0;

	static constexpr entry_type                one = static_cast<entry_type>(1);

public:

	nongrouped_qfilter_conc(size_t capacity = qf::DEFAULT_MIN_CAPACITY,
                            size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS,
                            const hasher& hf = hasher())
		: hf(hf),
          remainder_bits(remainder_bits),
          quotient_bits(qf::capacity_to_quotient_bits(capacity))
	{
		table_capacity = (one << quotient_bits) + qf::shift_buffer;
		table = std::make_unique<cell_atomic[]>(table_capacity);
		fingerprint_mask = (ONE << static_cast<size_t>(
                                       remainder_bits + quotient_bits)) - 1;
	}

	nongrouped_qfilter_conc(const nongrouped_qfilter_conc& other)
		: table_capacity(other.table_capacity),
          table(new cell_atomic[table_capacity]),
          table_status(other.table_status),
          hf(other.hf),
          fingerprint_mask(other.fingerprint_mask),
          remainder_bits(other.remainder_bits),
          quotient_bits(other.quotient_bits),
          growing_block_counter(other.growing_block_counter),
          finished_block_counter(other.finished_block_counter)
	{
		for (size_t i = 0; i < table_capacity; i++)
		{
			table[i].entry = other.table[i].entry;
		}
	}

	nongrouped_qfilter_conc& operator=(const nongrouped_qfilter_conc& other)
	{
		table_capacity = other.table_capacity;
		table = std::make_unique<cell_atomic[]>(other.table_capacity);
		table_status = other.table_status.load();
		hf = other.hf;
		fingerprint_mask = other.fingerprint_mask;
		remainder_bits = other.remainder_bits;
		quotient_bits = other.quotient_bits;
		growing_block_counter = other.growing_block_counter.load();
		finished_block_counter = other.finished_block_counter.load();

		for (size_t i = 0; i < table_capacity; i++)
		{
			table[i].entry = other.table[i].entry.load();
		}

		return *this;
	}

	nongrouped_qfilter_conc(nongrouped_qfilter_conc&&) = default;
	nongrouped_qfilter_conc& operator=(nongrouped_qfilter_conc&&) = default;

	bool insert(const Key& key)
	{
		return insert_hash(hf(key));
	}

	qf::InsertResult insert_hash(const HashVal& hashed)
	{
		if (table_status != qf::TableStatus::active)
			return qf::FAILED_OPERATION;

		auto[q, r] = get_quotient_and_remainder(hashed);

		if (try_empty_slot(q, r))
			return 1;

		write_lock_guard write_lock;
		if (!has_space_and_lock(q, write_lock))
			return qf::FAILED_OPERATION;

		bool run_already_existing = table[q].is_occupied();

		read_lock_guard read_lock;
		const size_t cluster_start = find_and_lock_cluster_start(q, read_lock);
		size_t run = find_run_start(cluster_start,
                                    q,
                                    read_lock,
                                    !run_already_existing);
		size_t insert_pos = run;

		cell new_entry;
		new_entry.set_remainder(r);

		if (insert_pos == write_lock.pos) {
			new_entry.set_shifted(true);
			table[insert_pos] = new_entry;
			write_lock.clear();
			return insert_pos - cluster_start;
		}
		else if (insert_pos == table_capacity - 1)
		{
			if (!run_already_existing)
				table[q].set_occupied(run_already_existing);

			return qf::FAILED_OPERATION;
		}

		cell current_entry = table[insert_pos];

		if (run_already_existing)
		{
			// find position in existing run
			do
			{
				if (current_entry.remainder() == r)
					return insert_pos - cluster_start; //element already present
				else if (current_entry.remainder() > r)
					break;

				++insert_pos;

				do
				{
					current_entry = table[insert_pos];
				} while (current_entry.is_read_locked());

			} while (!current_entry.is_locked()
                     && current_entry.is_continuation());

			if (insert_pos == run)
			{
				// note: possibly writing to read locked entry but only
                //       modifying the remainder
				table[run].set_remainder(r);
				new_entry.set_remainder(current_entry.remainder());

				++insert_pos;

				do
				{
					current_entry = table[insert_pos];
				} while (current_entry.is_read_locked());
			}

			new_entry.set_continuation(true);
		}

		//TODO: remove this?
		if (current_entry.is_read_locked() || current_entry.is_cluster_start())
			acquire_read_lock(insert_pos, read_lock);

		new_entry.set_shifted(insert_pos != q);

		insert_and_shift(new_entry, insert_pos, write_lock, read_lock);
		return insert_pos - cluster_start;
	}

	bool quick_insert(const Key& key)
	{
		return quick_insert_hash(hf(key));
	}

	// try insert without shifting (does not use or check for locks)
	bool quick_insert_hash(const HashVal& hashed)
	{
		auto[q, r] = get_quotient_and_remainder(hashed);

		if (unsafe_try_empty_slot(q, r))
			return 1;

		return unsafe_contains_hash(hashed);
	}

	qf::ContainsResult contains(const Key& key)
	{
		return contains_hash(hf(key));
	}

	qf::ContainsResult contains_hash(const HashVal& hashed)
	{
		if (table_status != qf::TableStatus::active)
			return qf::FAILED_OPERATION;

		auto[q, r] = get_quotient_and_remainder(hashed);
		cell current_entry = table[q];

		if (current_entry.is_empty())
			return ContainsResultEnum::CANONICAL_SLOT_EMPTY;
		if (!current_entry.is_occupied())
			return ContainsResultEnum::ENTRY_NOT_PRESENT;

		read_lock_guard read_lock;
		const size_t cluster_start = find_and_lock_cluster_start(q, read_lock);
		size_t iter = find_run_start(cluster_start, q, read_lock);

		current_entry = table[iter];

		do
		{
			// TODO: early abort because of sorting
			if (current_entry.remainder() == r)
				return ContainsResultEnum::ENTRY_PRESENT;

			iter++;

			do
			{
				current_entry = table[iter];
			} while (current_entry.is_read_locked());
		} while (!current_entry.is_write_locked()
                 && current_entry.is_continuation());

		return ContainsResultEnum::ENTRY_NOT_PRESENT;
	}

	qf::ContainsResult unsafe_contains(const Key& key)
	{
		return unsafe_contains_hash(hf(key));
	}

	// does not use or check for locks
	qf::ContainsResult unsafe_contains_hash(const HashVal& hashed) const
	{
		auto[q, r] = get_quotient_and_remainder(hashed);
		cell current_entry = table[q];

		if (current_entry.is_empty())
			return ContainsResultEnum::CANONICAL_SLOT_EMPTY;
		if (!current_entry.is_occupied())
			return ContainsResultEnum::ENTRY_NOT_PRESENT;

		const size_t cluster_start = find_cluster_start(q);
		size_t iter = unsafe_find_run_start(cluster_start, q);

		current_entry = table[iter];

		do
		{
			if (current_entry.remainder() == r)
				return ContainsResultEnum::ENTRY_PRESENT;

			current_entry = table[++iter];
		} while (current_entry.is_continuation());

		return ContainsResultEnum::ENTRY_NOT_PRESENT;
	}

	std::pair<size_t, entry_type> get_quotient_and_remainder(
        const HashVal& hashed) const
	{
		return qf::get_quotient_and_remainder<entry_type>(hashed, quotient_bits,
                                                          remainder_bits);
	}

	size_t capacity() const
	{
		return ONE << quotient_bits;
	}

	static size_t capacity(size_t min_capacity, size_t remainder_bits)
	{
		return qf::next_power_of_two(min_capacity);
	}

	void print()
	{
		print(0, table_capacity);
	}

	void print(size_t start, size_t count)
	{
		std::cout << "occu | cont | shif || qoutient | remainder\n";
		const size_t end = start + count;

		for (size_t i = start; i < end; i++)
		{
			print_entry(i);
			std::cout << "\n";
		}

		std::cout << std::endl;
	}

	void print_entry(size_t index)
	{
		print_entry(static_cast<cell>(table[index]), index);
	}

	void print_entry(cell entry, size_t index)
	{
		std::cout << entry.is_occupied() << "|"
                  << entry.is_continuation() << "|"
                  << entry.is_shifted() << " || "
                  << qf::bitstring(index, quotient_bits) << "|"
                  << qf::bitstring(entry.remainder(), remainder_bits);
	}

	void grow(nongrouped_qfilter_conc& target, size_t block_size)
	{
		if (table_status != qf::TableStatus::growing)
			return;

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

	bool check_consistency()
	{
		bool consistent = true;
		bool sorted = true;
		cell last;
		cell current;

		for (size_t pos = 0; pos < table_capacity; pos++)
		{
			current = table[pos];

			if (current.is_read_locked())
			{
				std::cerr << "error: [" << pos << "] read locked" << std::endl;
				consistent = false;
			}
			else if (current.is_write_locked())
			{
				std::cerr << "error: [" << pos << "] write locked" << std::endl;
				consistent = false;
			}
			else if (current.is_continuation() && !current.is_shifted())
			{
				std::cerr << "error1: [" << pos << "] cont =/=> shift"
                          << std::endl;
				consistent = false;
			}
			else if (current.is_continuation()
                     && !last.is_shifted()
                     && !last.is_cluster_start())
			{
				std::cerr << "error2: [" << pos << "] cont =/=> shift"
                          << std::endl;
				consistent = false;
			}
			else if (current.is_continuation())
			{
				if (current.remainder() < last.remainder())
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

	void set_table_status(qf::TableStatus status)
	{
		table_status = status;
	}
};

} // namespace qf
