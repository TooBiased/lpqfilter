#pragma once
/*******************************************************************************
 * implementation/base_filter/standard/standard_qfilter_cell.hpp
 *
 * slot implementation for the quotient filter variant with grouped slots
 * (multiple slots per atomic non-templated)
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/



#include <atomic>
#include <cstddef>
#include <utility>
// #include "nt_compact_quotient_filter_entry.h"

namespace qf
{

// represents the table and cell index for the exact location of an entry
template <size_t cell_capacity> struct entry_position
{
    size_t table = 0;
    size_t cell  = 0;

    entry_position() = default;
    explicit entry_position(size_t offset)
        : table(offset / cell_capacity), cell(offset % cell_capacity)
    {
    }
    entry_position(size_t table, size_t cell) : table(table), cell(cell) {}

    template <size_t other_cell_capacity>
    entry_position(const entry_position<other_cell_capacity>& other)
        : entry_position(other.offset())
    {
    }

    bool operator==(const entry_position& other)
    {
        return table == other.table && cell == other.cell;
    }
    bool operator!=(const entry_position& other) { return !(*this == other); }
    bool operator<(const entry_position& other)
    {
        return this->table < other.table ||
               (this->table == other.table && this->cell < other.cell);
    }
    bool operator>(const entry_position& other)
    {
        return this->table > other.table ||
               (this->table == other.table && this->cell > other.cell);
    }
    bool operator<=(const entry_position& other) { return !(*this > other); }
    bool operator>=(const entry_position& other) { return !(*this < other); }

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

    friend std::ostream&
    operator<<(std::ostream& stream, const entry_position& pos)
    {
        return stream << "[" << pos.table << "," << pos.cell << "]";
    }

    size_t offset() const { return table * cell_capacity + cell; }

    static constexpr int
    distance(const entry_position& start, const entry_position& end)
    {
        return static_cast<int>(end.offset()) -
               static_cast<int>(start.offset());
    }
};




template <class cell_base_type> class cell_stuff
{
  public:
    static constexpr size_t base_type_bits = sizeof(cell_base_type) << 3;
    static constexpr cell_base_type one    = static_cast<cell_base_type>(1u);
    static constexpr cell_base_type occupied_mask = one << 0;
    static constexpr cell_base_type shifted_mask  = one << 2;

    cell_stuff(size_t remainder)
        : capacity(base_type_bits / (remainder + status_bits)),
          remainder_bits(remainder),
          entry_mask((~static_cast<cell_base_type>(0)) >>
                     (base_type_bits - (status_bits + remainder))),
          remainder_mask((one << remainder) - 1),
          full_occupied_mask(generate_full_occupied_mask(remainder)),
          full_shifted_mask(generate_full_shifted_mask(remainder))
    {
    }

    const size_t         capacity;
    const size_t         remainder_bits;
    const cell_base_type entry_mask;
    const cell_base_type remainder_mask;
    const cell_base_type full_occupied_mask;
    const cell_base_type full_shifted_mask;

    inline cell_base_type generate_full_occupied_mask(size_t remainder)
    {
        cell_base_type mask = occupied_mask;

        for (size_t i = 1; i < capacity; i++)
        {
            // two shifts to prevent shift overflow / undefined behavior
            mask <<= remainder;
            mask <<= status_bits;
            mask |= occupied_mask;
        }

        return mask;
    }

    inline cell_base_type generate_full_shifted_mask(size_t remainder)
    {
        cell_base_type mask = shifted_mask;

        for (size_t i = 1; i < capacity; i++)
        {
            // two shifts to prevent shift overflow / undefined behavior
            mask <<= remainder;
            mask <<= status_bits;
            mask |= shifted_mask;
        }

        return mask;
    }
};

template <class cell_type, class integral_base_type = cell_type>
struct standard_cell_base
{
  public:
    using cell_base_type = integral_base_type;
    using entry_type     = size_t; // min_entry_type<remainder_bits>;


    // members
    cell_type cell{0};

    static constexpr cell_base_type one = static_cast<cell_base_type>(1u);
    static constexpr cell_base_type status_mask   = (one << status_bits) - 1;
    static constexpr cell_base_type occupied_mask = one << 0;
    static constexpr cell_base_type continuation_mask = one << 1;
    static constexpr cell_base_type shifted_mask      = one << 2;

    static constexpr cell_base_type cluster_start_status = occupied_mask;
    static constexpr cell_base_type write_lock_status    = continuation_mask;
    static constexpr cell_base_type read_lock_status =
        occupied_mask | continuation_mask;

    using cstuff = cell_stuff<cell_base_type>;

  public:
    standard_cell_base() = default;
    standard_cell_base(const cell_type& cell) : cell(cell){};

    void clear() { cell = 0; }

    // getter

    entry_type entry(const cstuff& cc, size_t pos) const
    {
        return (cell >> shift_amount(cc, pos)) & cc.entry_mask;
    }
    entry_type remainder(const cstuff& cc, size_t pos) const
    {
        return (cell >> (shift_amount(cc, pos) + status_bits)) &
               cc.remainder_mask;
    }
    status_type status(const cstuff& cc, size_t pos) const
    {
        return static_cast<status_type>((cell >> shift_amount(cc, pos)) &
                                        status_mask);
    }
    bool is_occupied(const cstuff& cc, size_t pos) const
    {
        return status(cc, pos) & occupied_mask;
    }
    bool is_continuation(const cstuff& cc, size_t pos) const
    {
        return status(cc, pos) & continuation_mask;
    }
    bool is_shifted(const cstuff& cc, size_t pos) const
    {
        return status(cc, pos) & shifted_mask;
    }
    bool is_cluster_start(const cstuff& cc, size_t pos) const
    {
        return status(cc, pos) == cluster_start_status;
    }
    bool is_empty(const cstuff& cc, size_t pos) const
    {
        return status(cc, pos) == 0;
    }
    cell_base_type occupied_mask_full(const cstuff& cc) const
    {
        return cell & cc.full_occupied_mask;
    }
    bool is_read_locked(const cstuff& cc, size_t pos) const
    {
        return status(cc, pos) == read_lock_status;
    }
    bool
    is_read_locked(const cstuff& cc, size_t start_pos, size_t end_pos) const
    {
        for (size_t pos = start_pos; pos < end_pos; pos++)
        {
            if (status(cc, pos) == read_lock_status) return true;
        }
        return false;
    }
    bool is_write_locked(const cstuff& cc, size_t pos) const
    {
        return status(cc, pos) == write_lock_status;
    }
    bool is_locked(const cstuff& cc, size_t pos) const
    {
        const status_type local_status = status(cc, pos);
        return (local_status == read_lock_status) ||
               (local_status == write_lock_status);
    }

    // static entry getter - with cc
    static entry_type entry_remainder(const cstuff& cc, entry_type entry)
    {
        return (entry >> status_bits) & cc.remainder_mask;
    }
    static status_type entry_status(const cstuff&, entry_type entry)
    {
        return static_cast<status_type>(entry & status_mask);
    }
    static bool entry_is_occupied(const cstuff& cc, entry_type entry)
    {
        return entry_status(cc, entry) & occupied_mask;
    }
    static bool entry_is_continuation(const cstuff& cc, entry_type entry)
    {
        return entry_status(cc, entry) & continuation_mask;
    }
    static bool entry_is_shifted(const cstuff& cc, entry_type entry)
    {
        return entry_status(cc, entry) & shifted_mask;
    }
    static bool entry_is_cluster_start(const cstuff& cc, entry_type entry)
    {
        return entry_status(cc, entry) == cluster_start_status;
    }
    static bool entry_is_empty(const cstuff& cc, entry_type entry)
    {
        return entry_status(cc, entry) == 0;
    }
    static bool entry_is_read_locked(const cstuff& cc, entry_type entry)
    {
        return entry_status(cc, entry) == read_lock_status;
    }
    static bool entry_is_write_locked(const cstuff& cc, entry_type entry)
    {
        return entry_status(cc, entry) == write_lock_status;
    }
    static bool entry_is_locked(const cstuff& cc, entry_type entry)
    {
        const status_type local_status = entry_status(cc, entry);
        return local_status == read_lock_status ||
               local_status == write_lock_status;
    }

    // setter with cstuff
    static void entry_set_remainder(const cstuff& cc, entry_type& entry,
                                    entry_type remainder)
    {
        entry &= ~(cc.remainder_mask << status_bits);
        entry |= remainder << status_bits;
    }

    static void
    entry_set_status(const cstuff&, entry_type& entry, status_type status)
    {
        entry &= ~status_mask;
        entry |= status;
    }

    static void entry_set_occupied(const cstuff&, entry_type& entry, bool val)
    {
        entry &= ~occupied_mask;
        entry |= occupied_mask * static_cast<int>(val);
    }

    static void
    entry_set_continuation(const cstuff&, entry_type& entry, bool val)
    {
        entry &= ~continuation_mask;
        entry |= continuation_mask * static_cast<int>(val);
    }

    static void entry_set_shifted(const cstuff&, entry_type& entry, bool val)
    {
        entry &= ~shifted_mask;
        entry |= shifted_mask * static_cast<int>(val);
    }

    static void
    status_set_occupied(const cstuff&, status_type& status, bool val)
    {
        status &= ~occupied_mask;
        status |= occupied_mask * static_cast<int>(val);
    }

    static void
    status_set_continuation(const cstuff&, status_type& status, bool val)
    {
        status &= ~continuation_mask;
        status |= continuation_mask * static_cast<int>(val);
    }

    static void status_set_shifted(const cstuff&, status_type& status, bool val)
    {
        status &= ~shifted_mask;
        status |= shifted_mask * static_cast<int>(val);
    }

    // misc
    size_t shift_amount(const cstuff& cc, size_t pos) const
    {
        return pos * (cc.remainder_bits + status_bits);
    }
};

template <class cell_type> struct standard_cell_atomic;

template <class cell_type>
struct standard_cell_non_atomic : public standard_cell_base<cell_type>
{
    using atomic = standard_cell_atomic<cell_type>;
    using base   = standard_cell_base<cell_type>;

    using base::base;
    using typename base::cstuff;
    using typename base::entry_type;


    standard_cell_non_atomic()                                      = default;
    standard_cell_non_atomic(const standard_cell_non_atomic& other) = default;

    standard_cell_non_atomic(const atomic& other) : base(other.cell.load()) {}

    standard_cell_non_atomic&
    operator=(const standard_cell_non_atomic&) = default;

    standard_cell_non_atomic& operator=(const atomic& other)
    {
        this->cell = other.cell.load();
        return *this;
    }

    // setter

    void set_entry(const cstuff& cc, entry_type entry, size_t pos)
    {
        this->cell &= ~(cc.entry_mask << base::shift_amount(cc, pos));
        this->cell |= static_cast<cell_type>(entry)
                      << base::shift_amount(cc, pos);
    }

    void set_remainder(const cstuff& cc, entry_type remainder, size_t pos)
    {
        this->cell &=
            ~(cc.remainder_mask << (base::shift_amount(cc, pos) + status_bits));
        this->cell |= static_cast<cell_type>(remainder)
                      << (base::shift_amount(cc, pos) + status_bits);
    }

    void set_status(const cstuff& cc, status_type status, size_t pos)
    {
        this->cell &= ~(base::status_mask << base::shift_amount(cc, pos));
        this->cell |= static_cast<cell_type>(status)
                      << base::shift_amount(cc, pos);
    }

    void set_continuation(const cstuff& cc, bool val, size_t pos)
    {
        const cell_type mask = base::continuation_mask
                               << base::shift_amount(cc, pos);
        this->cell = (this->cell & ~mask) | (mask * static_cast<int>(val));
    }

    void set_shifted(const cstuff& cc, bool val, size_t pos)
    {
        const cell_type mask = base::shifted_mask
                               << base::shift_amount(cc, pos);
        this->cell = (this->cell & ~mask) | (mask * static_cast<int>(val));
    }

    void set_occupied(const cstuff& cc, bool val, size_t pos)
    {
        const cell_type mask = base::occupied_mask
                               << base::shift_amount(cc, pos);
        this->cell = (this->cell & ~mask) | (mask * static_cast<int>(val));
    }

    void set_occupied_full(const cstuff& cc, cell_type occupied_values)
    {
        this->cell = (this->cell & ~cc.full_occupied_mask) |
                     (cc.full_occupied_mask & occupied_values);
    }

    void set_read_lock(const cstuff& cc, size_t pos)
    {
        this->cell &= ~(base::status_mask << base::shift_amount(cc, pos));
        this->cell |= base::read_lock_status << base::shift_amount(cc, pos);
    }

    void set_write_lock(const cstuff& cc, size_t pos)
    {
        this->cell &= ~(base::status_mask << base::shift_amount(cc, pos));
        this->cell |= base::write_lock_status << base::shift_amount(cc, pos);
    }
};

template <class cell_type>
struct standard_cell_atomic
    : standard_cell_base<std::atomic<cell_type>, cell_type>
{
    using non_atomic = standard_cell_non_atomic<cell_type>;
    using base       = standard_cell_base<std::atomic<cell_type>, cell_type>;

    using base::base;
    using typename base::cstuff;
    using typename base::entry_type;

    standard_cell_atomic() = default;
    standard_cell_atomic(const standard_cell_atomic& other)
        : base(other.cell.load())
    {
    }

    standard_cell_atomic(const non_atomic& other) : base(other.cell) {}

    standard_cell_atomic& operator=(const standard_cell_atomic& other)
    {
        this->cell.store(other.cell.load());
        return *this;
    }

    standard_cell_atomic& operator=(const non_atomic& other)
    {
        this->cell.store(other.cell);
        return *this;
    }

    operator non_atomic() { return {this->cell.load()}; }

    // setter

    void set_entry(const cstuff& cc, entry_type entry, size_t pos)
    {
        update_cell(static_cast<cell_type>(entry)
                        << base::shift_amount(cc, pos),
                    cc.entry_mask << base::shift_amount(cc, pos));
    }

    void set_remainder(const cstuff& cc, entry_type remainder, size_t pos)
    {
        update_cell(static_cast<cell_type>(remainder)
                        << (base::shift_amount(cc, pos) + status_bits),
                    cc.remainder_mask
                        << (base::shift_amount(cc, pos) + status_bits));
    }

    void set_status(const cstuff& cc, status_type status, size_t pos)
    {
        update_cell(static_cast<cell_type>(status)
                        << base::shift_amount(cc, pos),
                    base::status_mask << base::shift_amount(cc, pos));
    }

    void set_continuation(const cstuff& cc, bool val, size_t pos)
    {
        const cell_type mask = base::continuation_mask
                               << base::shift_amount(cc, pos);
        update_cell(mask * static_cast<int>(val), mask);
    }

    void set_shifted(const cstuff& cc, bool val, size_t pos)
    {
        const cell_type mask = base::shifted_mask
                               << base::shift_amount(cc, pos);
        update_cell(mask * static_cast<int>(val), mask);
    }

    void set_occupied(const cstuff& cc, bool val, size_t pos)
    {
        const cell_type mask = base::occupied_mask
                               << base::shift_amount(cc, pos);
        update_cell(mask * static_cast<int>(val), mask);
    }

    void set_read_lock(const cstuff& cc, size_t pos)
    {
        set_status(cc, base::read_lock_status, pos);
    }

    void set_write_lock(const cstuff& cc, size_t pos)
    {
        set_status(cc, base::write_lock_status, pos);
    }

    void update_cell(const cell_type& val, const cell_type& mask)
    {
        cell_type cell;
        cell_type new_cell;

        do {
            cell     = this->cell.load();
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


enum class lock_tag
{
    undefined_tag,
    read_lock_tag,
    write_lock_tag
};

template <class cell_type, lock_tag tag> class standard_cell_lock_guard
{
  public:
    using this_type = standard_cell_lock_guard<cell_type, tag>;
    using entry_pos = std::pair<size_t, size_t>;
    using cstuff    = typename cell_type::cstuff;

    standard_cell_lock_guard(const cstuff& cc) : cc(cc) {}
    standard_cell_lock_guard(const cstuff& cc, entry_pos pos,
                             cell_type& locked_cell, status_type status)
        : cc(cc), pos(pos), cell(&locked_cell), status(status){};
    ~standard_cell_lock_guard() { release(); }

    standard_cell_lock_guard(const standard_cell_lock_guard&) = delete;
    standard_cell_lock_guard&
    operator=(const standard_cell_lock_guard&) = delete;

    standard_cell_lock_guard(standard_cell_lock_guard&& other)
        : cc(other.cc), pos(std::move(other.pos)), cell(std::move(other.cell)),
          status(std::move(other.status))
    {
        other.clear();
    }

    standard_cell_lock_guard& operator=(standard_cell_lock_guard&& other)
    {
        if (this == &other) return *this;

        this->~this_type();
        new (this) this_type(other);
        return *this;
    }

    void set_data(cell_type* cell, const entry_pos& pos, status_type status)
    {
        if (this->pos != pos)
        {
            release();
            this->pos = pos;
        }

        this->cell   = cell;
        this->status = status;
    }

    bool locked() const { return cell != nullptr; }
    void clear() { cell = nullptr; }

    void release()
    {
        if (cell)
        {
            cell->set_status(cc, status, pos.second);
            cell = nullptr;
        }
    }

    const cstuff& cc;
    entry_pos     pos;
    cell_type*    cell = nullptr;
    status_type   status{};
};

template <class cell_type>
using compact_read_lock_guard =
    standard_cell_lock_guard<cell_type, lock_tag::read_lock_tag>;

template <class cell_type>
using compact_write_lock_guard =
    standard_cell_lock_guard<cell_type, lock_tag::write_lock_tag>;


} // namespace qf
