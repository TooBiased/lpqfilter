#pragma once
/*******************************************************************************
 * implementation/base_filter/templated/templated_lpfilter_cell.hpp
 *
 * templated grouped slot implemtation for the linear probing filter
 * (multiple slots per data element templated)
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/



#include <atomic>
#include "templated_qfilter_cell.hpp"

namespace qf {

template<class cell_type, size_t remainder_bits,
         class intgeral_base_type = cell_type>
struct templated_lpfilter_cell_base
{
public:
	using cell_base_type = intgeral_base_type;
    using remainder_type = min_remainder_type<remainder_bits>;

public:
	templated_lpfilter_cell_base() = default;
	templated_lpfilter_cell_base(const cell_type& cell) : cell(cell) {};

    // getter

	remainder_type remainder(size_t pos) const
	{
		return (cell >> shift_amount(pos)) & remainder_mask;
	}

	bool is_empty(size_t pos) const
	{
		return remainder(pos) == 0;
	}

	static bool remainder_is_empty(remainder_type remainder)
	{
		return remainder == 0;
	}

	// setter

	void clear() {
        cell = 0;
    }

    // misc
    size_t shift_amount(size_t pos) const {
        return pos * remainder_bits;
    }

	// members

	cell_type cell{ 0 };

    static constexpr size_t          capacity = sizeof(cell_base_type) * 8u
                                                               / remainder_bits;
    static constexpr cell_base_type  one      = static_cast<cell_base_type>(1u);
	static constexpr cell_base_type  remainder_mask = (one << remainder_bits)-1;
};

template<class cell_type, size_t remainder_bits>
struct templated_lpfilter_cell_atomic;

template<class cell_type, size_t remainder_bits>
struct templated_lpfilter_cell_non_atomic
    : public templated_lpfilter_cell_base<cell_type, remainder_bits>
{
	using atomic = templated_lpfilter_cell_atomic<cell_type, remainder_bits>;
	using base = templated_lpfilter_cell_base<cell_type, remainder_bits>;

    using typename base::remainder_type;
	using base::base;

	templated_lpfilter_cell_non_atomic() = default;
	templated_lpfilter_cell_non_atomic(
        const templated_lpfilter_cell_non_atomic&) = default;
	templated_lpfilter_cell_non_atomic(const atomic& other)
        : base(other.cell.load())
    {}

	templated_lpfilter_cell_non_atomic& operator=(
        const templated_lpfilter_cell_non_atomic&) = default;
	templated_lpfilter_cell_non_atomic& operator=(const atomic& other)
	{
		this->cell = other.cell.load();
		return *this;
	}

	// setter

	void set_remainder(remainder_type remainder, size_t pos)
	{
		this->cell &= ~(base::remainder_mask << base::shift_amount(pos));
		this->cell |= static_cast<cell_type>(remainder)
                      << base::shift_amount(pos);
	}
};

template<class cell_type, size_t remainder_bits>
struct templated_lpfilter_cell_atomic
    : templated_lpfilter_cell_base<std::atomic<cell_type>,
                                   remainder_bits,
                                   cell_type>
{
	using non_atomic = templated_lpfilter_cell_non_atomic<cell_type,
                                                          remainder_bits>;
	using base = templated_lpfilter_cell_base<std::atomic<cell_type>,
                                        remainder_bits,
                                        cell_type>;

    using typename base::remainder_type;
	using base::base;

	templated_lpfilter_cell_atomic() = default;
	templated_lpfilter_cell_atomic(const templated_lpfilter_cell_atomic& other)
        : base(other.cell.load())
    {}
	templated_lpfilter_cell_atomic(const non_atomic& other)
        : base(other.cell)
    {}

	templated_lpfilter_cell_atomic& operator=(
        const templated_lpfilter_cell_atomic& other)
	{
		this->cell.store(other.cell.load());
		return *this;
	}

	templated_lpfilter_cell_atomic& operator=(const non_atomic& other)
	{
		this->cell.store(other.cell);
		return *this;
	}

	operator non_atomic()
	{
		return { this->cell.load() };
	}

	// setter

	void set_remainder(remainder_type remainder, size_t pos)
	{
		update_cell(static_cast<cell_type>(remainder) <<base::shift_amount(pos),
                    base::remainder_mask <<base::shift_amount(pos));
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
} // namespace qf
