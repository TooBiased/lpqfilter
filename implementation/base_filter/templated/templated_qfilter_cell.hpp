#pragma once
/*******************************************************************************
 * implementation/base_filter/templated/templated_qfilter_cell.hpp
 *
 * templated implementation of the grouped slot concept for quotient filters
 * (multiple slots per data element templated)
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/



#include <cstddef>
#include <atomic>

#include "implementation/definitions.hpp"

//#include "compact_quotient_filter_entry.h"

namespace qf {



// chooses the minimum type necessary to hold a certain number of bits
template <size_t entry_bits>
struct entry_type_chooser
{ using type = typename entry_type_chooser<entry_bits - 1>::type; };

template <>
struct entry_type_chooser<0>  { using type = std::uint8_t; };

template <>
struct entry_type_chooser<9>  { using type = std::uint16_t; };

template <>
struct entry_type_chooser<17> { using type = std::uint32_t; };

template <>
struct entry_type_chooser<33> { using type = std::uint64_t; };

template <size_t remainder_bits>
using min_entry_type = typename entry_type_chooser<remainder_bits
                                                   + status_bits>::type;
template <size_t remainder_bits>
using min_remainder_type = typename entry_type_chooser<remainder_bits>::type;






// represents the table and cell index for the exact location of an entry
template<size_t cell_capacity>
struct entry_position
{
	size_t table = 0;
	size_t cell = 0;

	entry_position() = default;
	explicit entry_position(size_t offset)
        : table(offset / cell_capacity), cell(offset % cell_capacity)
    {}
	entry_position(size_t table, size_t cell)
        : table(table), cell(cell)
    {}

	template<size_t other_cell_capacity>
	entry_position(const entry_position<other_cell_capacity>& other)
        : entry_position(other.offset())
    {}

	bool operator==(const entry_position& other)
    { return table == other.table && cell == other.cell; }
	bool operator!=(const entry_position& other)
    { return !(*this == other); }
	bool operator< (const entry_position& other)
    {
        return this->table < other.table
                 || (this->table == other.table && this->cell < other.cell);
    }
	bool operator> (const entry_position& other)
    {
        return this->table > other.table
                 || (this->table == other.table && this->cell > other.cell);
    }
	bool operator<=(const entry_position& other)
    { return !(*this > other);	}
	bool operator>=(const entry_position& other)
    { return !(*this < other); }

	entry_position& operator++()
	{
		if (++cell == cell_capacity)
		{
			cell = 0;
			table++;
		}

		return *this;
	}

	entry_position& operator--()
	{
		if (cell-- == 0)
		{
			cell = cell_capacity - 1;
			table--;
		}

		return *this;
	}

	entry_position operator++(int)
	{
		entry_position previous = *this;

		if (++cell == cell_capacity)
		{
			cell = 0;
			table++;
		}

		return previous;
	}

	entry_position operator--(int)
	{
		entry_position previous = *this;

		if (cell-- == 0)
		{
			cell = cell_capacity - 1;
			table--;
		}

		return previous;
	}

	friend std::ostream& operator<<(std::ostream& stream,
                                    const entry_position& pos)
    {
		return stream << "[" << pos.table << "," << pos.cell << "]";
	}

	size_t offset() const
	{
		return table * cell_capacity + cell;
	}

	static constexpr int distance(const entry_position& start,
                                  const entry_position& end)
	{
		return static_cast<int>(end.offset()) -static_cast<int>(start.offset());
	}
};






template<class cell_type, size_t remainder_bits,
         class intgeral_base_type = cell_type>
struct templated_cell_base
{
public:
	using cell_base_type = intgeral_base_type;
    using entry_type = min_entry_type<remainder_bits>;

//private:
	static constexpr cell_base_type generate_full_occupied_mask()
    {
		cell_base_type mask = occupied_mask;

		for (size_t i = 1; i < capacity; i++)
		{
			mask <<= remainder_bits; // two shifts to prevent shift overflow / undefined behavior
			mask <<= status_bits;
			mask |= occupied_mask;
		}

		return mask;
	}

	static constexpr cell_base_type generate_full_shifted_mask()
    {
		cell_base_type mask = shifted_mask;

		for (size_t i = 1; i < capacity; i++)
		{
			mask <<= remainder_bits; // two shifts to prevent shift overflow / undefined behavior
			mask <<= status_bits;
			mask |= shifted_mask;
		}

		return mask;
	}

public:
	templated_cell_base() = default;
	templated_cell_base(const cell_type& cell) : cell(cell) {};

    // getter

    entry_type entry(size_t pos) const
	{ return (cell >> shift_amount(pos)) & entry_mask; }

	entry_type remainder(size_t pos) const
	{ return (cell >> (shift_amount(pos) + status_bits)) & remainder_mask; }

	static entry_type entry_remainder(entry_type entry)
	{ return (entry >> status_bits) & remainder_mask; }

	status_type status(size_t pos) const
	{
        return static_cast<status_type>((cell >> shift_amount(pos))
                                        & status_mask);
    }

	static status_type entry_status(entry_type entry)
	{ return static_cast<status_type>(entry & status_mask); }

	bool is_occupied(size_t pos) const
    { return status(pos) & occupied_mask; }

	static bool entry_is_occupied(entry_type entry)
	{ return entry_status(entry) & occupied_mask; }

	bool is_continuation(size_t pos) const
	{ return status(pos) & continuation_mask; }

	static bool entry_is_continuation(entry_type entry)
	{ return entry_status(entry) & continuation_mask; }

	bool is_shifted(size_t pos) const
	{ return status(pos) & shifted_mask; }

	static bool entry_is_shifted(entry_type entry)
	{ return entry_status(entry) & shifted_mask; }

	bool is_cluster_start(size_t pos) const
	{ return status(pos) == cluster_start_status; }

	static bool entry_is_cluster_start(entry_type entry)
	{ return entry_status(entry) == cluster_start_status; }

	bool is_empty(size_t pos) const
	{ return status(pos) == 0; }

	static bool entry_is_empty(entry_type entry)
	{ return entry_status(entry) == 0; }

	cell_base_type occupied_mask_full() const
	{ return cell & full_occupied_mask; }

	bool is_read_locked(size_t pos) const
	{ return status(pos) == read_lock_status; }

	bool is_read_locked(size_t start_pos, size_t end_pos) const
	{
        for (size_t pos = start_pos; pos < end_pos; pos++)
		{
			if (status(pos) == read_lock_status)
				return true;
		}

		return false;
	}

	static bool entry_is_read_locked(entry_type entry)
	{ return entry_status(entry) == read_lock_status; }

	bool is_write_locked(size_t pos) const
	{ return status(pos) == write_lock_status; }

	static bool entry_is_write_locked(entry_type entry)
	{ return entry_status(entry) == write_lock_status; }

	bool is_locked(size_t pos) const
	{
        const status_type local_status = status(pos);
        return local_status == read_lock_status
               || local_status == write_lock_status;
	}

	static bool entry_is_locked(entry_type entry)
	{
		const status_type local_status = entry_status(entry);
		return local_status == read_lock_status
               || local_status == write_lock_status;
	}

	// setter

	void clear() {
        cell = 0;
    }

	static void entry_set_remainder(entry_type& entry, entry_type remainder)
	{
		entry &= ~(remainder_mask << status_bits);
		entry |= remainder << status_bits;
	}

	static void entry_set_status(entry_type& entry, status_type status)
	{
		entry &= ~status_mask;
		entry |= status;
	}

	static void entry_set_occupied(entry_type& entry, bool val)
	{
		entry &= ~occupied_mask;
		entry |= occupied_mask * static_cast<int>(val);
	}

	static void entry_set_continuation(entry_type& entry, bool val)
	{
		entry &= ~continuation_mask;
		entry |= continuation_mask * static_cast<int>(val);
	}

	static void entry_set_shifted(entry_type& entry, bool val)
	{
		entry &= ~shifted_mask;
		entry |= shifted_mask * static_cast<int>(val);
	}

	static void status_set_occupied(status_type& status, bool val)
	{
		status &= ~occupied_mask;
		status |= occupied_mask * static_cast<int>(val);
	}

	static void status_set_continuation(status_type& status, bool val)
	{
		status &= ~continuation_mask;
		status |= continuation_mask * static_cast<int>(val);
	}

	static void status_set_shifted(status_type& status, bool val)
	{
		status &= ~shifted_mask;
		status |= shifted_mask * static_cast<int>(val);
	}

    // misc
    size_t shift_amount(size_t pos) const {
        return pos * (remainder_bits + status_bits);
    }

	// members

	cell_type cell{ 0 };

    static constexpr size_t          capacity             = sizeof(cell_base_type) * 8u / (remainder_bits + status_bits);
    static constexpr cell_base_type  one                  = static_cast<cell_base_type>(1u);

	static constexpr cell_base_type  entry_mask           = (~static_cast<cell_base_type>(0)) >> (sizeof(cell_base_type) * 8u - (status_bits + remainder_bits));
	static constexpr cell_base_type  remainder_mask       = (one << remainder_bits) - 1;
	static constexpr cell_base_type  status_mask          = (one << status_bits) - 1;
	static constexpr cell_base_type  occupied_mask        = one << 0;
	static constexpr cell_base_type  continuation_mask    = one << 1;
	static constexpr cell_base_type  shifted_mask         = one << 2;
	static constexpr cell_base_type  full_occupied_mask   = generate_full_occupied_mask();
	static constexpr cell_base_type  full_shifted_mask    = generate_full_shifted_mask();

	static constexpr cell_base_type  cluster_start_status = occupied_mask;
	static constexpr cell_base_type  write_lock_status    = continuation_mask;
	static constexpr cell_base_type  read_lock_status     = occupied_mask | continuation_mask;
};






template<class cell_type, size_t remainder_bits>
struct templated_cell_atomic;

template<class cell_type, size_t remainder_bits>
struct templated_cell_non_atomic
    : public templated_cell_base<cell_type, remainder_bits>
{
	using atomic = templated_cell_atomic<cell_type, remainder_bits>;
	using base = templated_cell_base<cell_type, remainder_bits>;

    using typename base::entry_type;
	using base::base;

	templated_cell_non_atomic() = default;
	templated_cell_non_atomic(const templated_cell_non_atomic& other) = default;

	templated_cell_non_atomic(const atomic& other) : base(other.cell.load()) {}

	templated_cell_non_atomic& operator=(const templated_cell_non_atomic&) = default;

	templated_cell_non_atomic& operator=(const atomic& other)
	{
		this->cell = other.cell.load();
		return *this;
	}

	// setter

    void set_entry(entry_type entry, size_t pos)
	{
		this->cell &= ~(base::entry_mask << base::shift_amount(pos));
		this->cell |= static_cast<cell_type>(entry) << base::shift_amount(pos);
	}

	void set_remainder(entry_type remainder, size_t pos)
	{
		this->cell &= ~(base::remainder_mask
                          << (base::shift_amount(pos) + status_bits));
		this->cell |= static_cast<cell_type>(remainder)
                          << (base::shift_amount(pos) + status_bits);
	}

	void set_status(status_type status, size_t pos)
	{
		this->cell &= ~(base::status_mask << base::shift_amount(pos));
		this->cell |= static_cast<cell_type>(status) << base::shift_amount(pos);
	}

	void set_continuation(bool val, size_t pos)
	{
		const cell_type mask = base::continuation_mask
                                 << base::shift_amount(pos);
		this->cell = (this->cell & ~mask) | (mask * static_cast<int>(val));
	}

	void set_shifted(bool val, size_t pos)
	{
		const cell_type mask = base::shifted_mask << base::shift_amount(pos);
		this->cell = (this->cell & ~mask) | (mask * static_cast<int>(val));
	}

	void set_occupied(bool val, size_t pos)
	{
		const cell_type mask = base::occupied_mask << base::shift_amount(pos);
		this->cell = (this->cell & ~mask) | (mask * static_cast<int>(val));
	}

	void set_occupied_full(cell_type occupied_values)
	{
		this->cell = (this->cell & ~base::full_occupied_mask)
                       | (base::full_occupied_mask & occupied_values);
	}

	void set_read_lock(size_t pos)
	{
		this->cell &= ~(base::status_mask << base::shift_amount(pos));
		this->cell |= base::read_lock_status << base::shift_amount(pos);
	}

	void set_write_lock(size_t pos)
	{
		this->cell &= ~(base::status_mask << base::shift_amount(pos));
		this->cell |= base::write_lock_status << base::shift_amount(pos);
	}
};





template<class cell_type, size_t remainder_bits>
struct templated_cell_atomic
    : templated_cell_base<std::atomic<cell_type>, remainder_bits, cell_type>
{
	using non_atomic = templated_cell_non_atomic<cell_type, remainder_bits>;
	using base = templated_cell_base<std::atomic<cell_type>,
                                     remainder_bits,
                                     cell_type>;

    using typename base::entry_type;
	using base::base;

	templated_cell_atomic() = default;
	templated_cell_atomic(const templated_cell_atomic& other)
        : base(other.cell.load())
    {}
	templated_cell_atomic(const non_atomic& other)
        : base(other.cell)
    {}

	templated_cell_atomic& operator=(const templated_cell_atomic& other)
	{
		this->cell.store(other.cell.load());
		return *this;
	}

	templated_cell_atomic& operator=(const non_atomic& other)
	{
		this->cell.store(other.cell);
		return *this;
	}

	operator non_atomic()
	{
		return { this->cell.load() };
	}

	// setter

	void set_entry(entry_type entry, size_t pos)
	{
		update_cell(static_cast<cell_type>(entry) << base::shift_amount(pos),
                    base::entry_mask << base::shift_amount(pos));
	}

	void set_remainder(entry_type remainder, size_t pos)
	{
		update_cell(static_cast<cell_type>(remainder)
                      << (base::shift_amount(pos) + status_bits),
                    base::remainder_mask
                      << (base::shift_amount(pos) + status_bits));
	}

	void set_status(status_type status, size_t pos)
	{
		update_cell(static_cast<cell_type>(status) << base::shift_amount(pos),
                    base::status_mask << base::shift_amount(pos));
	}

	void set_continuation(bool val, size_t pos)
	{
		const cell_type mask = base::continuation_mask
                                 << base::shift_amount(pos);
		update_cell(mask * static_cast<int>(val), mask);
	}

	void set_shifted(bool val, size_t pos)
	{
		const cell_type mask = base::shifted_mask << base::shift_amount(pos);
		update_cell(mask * static_cast<int>(val), mask);
	}

	void set_occupied(bool val, size_t pos)
	{
		const cell_type mask = base::occupied_mask << base::shift_amount(pos);
		update_cell(mask * static_cast<int>(val), mask);
	}

	void set_read_lock(size_t pos)
	{
		set_status(base::read_lock_status, pos);
	}

	void set_write_lock(size_t pos)
	{
		set_status(base::write_lock_status, pos);
	}

	void update_cell(const cell_type& val, const cell_type& mask)
	{
		cell_type cell;
		cell_type new_cell;

		do
		{
			cell = this->cell.load();
			new_cell = (cell & ~mask) | val;
		} while (!CAS(cell, new_cell));
	}

    bool CAS(cell_type& expected, const cell_type& desired)
	{
		return this->cell.compare_exchange_strong(expected, desired);
	}

	bool CAS(non_atomic& expected, const non_atomic& desired)
	{
		return CAS(expected.cell, desired.cell);
	}

};






template<class cell_type>
struct templated_cell_lock_guard
{
	using entry_pos = entry_position<cell_type::capacity>;

	templated_cell_lock_guard() = default;
	templated_cell_lock_guard(entry_pos pos,
                              cell_type& locked_cell,
                              status_type status)
        : pos(pos), cell(&locked_cell), status(status)
    {};
	~templated_cell_lock_guard() { release(); }

	templated_cell_lock_guard(const templated_cell_lock_guard&) = delete;
	templated_cell_lock_guard& operator= (const templated_cell_lock_guard&) = delete;

	templated_cell_lock_guard(templated_cell_lock_guard&& other)
        : pos(std::move(other.pos)),
          cell(std::move(other.cell)),
          status(std::move(other.status))
	{
		other.clear();
	};

	templated_cell_lock_guard& operator= (templated_cell_lock_guard&& other)
	{
		if (pos != other.pos)
			release();

		pos = std::move(other.pos);
		cell = std::move(other.cell);
		status = std::move(other.status);

		other.clear();
		return *this;
	}

	void set_data(cell_type* cell, const entry_pos& pos, status_type status)
	{
		if (this->pos != pos)
		{
			release();
			this->pos = pos;
		}

		this->cell = cell;
		this->status = status;
	}

	bool locked() const { return cell != nullptr; }
	void clear() { cell = nullptr; }

	void release()
	{
		if (cell)
		{
			cell->set_status(status, pos.cell);
			cell = nullptr;
		}
	}

	entry_pos pos;
	cell_type* cell = nullptr;
	status_type status{};
};

template<class cell_type>
struct templated_read_lock_guard  : templated_cell_lock_guard<cell_type> {};

template<class cell_type>
struct templated_write_lock_guard : templated_cell_lock_guard<cell_type> {};



} // namespace qf
