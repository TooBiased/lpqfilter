#pragma once
#include <atomic>
#include "standard_qfilter_cell.hpp"

namespace qf {

template<class cell_base_type>
class lp_cell_stuff
{
public:
    static constexpr size_t base_type_bits = sizeof(cell_base_type) << 3;
    static constexpr cell_base_type  one   = static_cast<cell_base_type>(1u);

    lp_cell_stuff(size_t remainder)
        : capacity(base_type_bits/remainder),
          remainder_bits(remainder),
          remainder_mask((one<<remainder)-1)
    { }

    const size_t         capacity;
    const size_t         remainder_bits;
    const cell_base_type remainder_mask;

};

template<class cell_type, class intgeral_base_type = cell_type>
struct standard_lpfilter_cell_base
{
public:
	using cell_base_type = intgeral_base_type;
    using entry_type     = size_t;

    using cstuff = lp_cell_stuff<cell_base_type>;
public:

	cell_type cell{ 0 };

    static constexpr cell_base_type  one = static_cast<cell_base_type>(1u);

	standard_lpfilter_cell_base() = default;
	standard_lpfilter_cell_base(const cell_type& cell) : cell(cell) {};

	inline void clear() {
        cell = 0;
    }

    // getter

	inline entry_type remainder(const cstuff& cc, size_t pos) const
	{
		return (cell >> shift_amount(cc, pos)) & cc.remainder_mask;
	}

	inline bool is_empty(const cstuff& cc, size_t pos) const
	{
		return remainder(cc, pos) == 0;
	}

	inline static bool remainder_is_empty(entry_type remainder)
	{
		return remainder == 0;
	}

    // misc
    inline size_t shift_amount(const cstuff& cc, size_t pos) const {
        return pos * cc.remainder_bits;
    }
};

template<class cell_type>
struct standard_lpfilter_cell_atomic;

template<class cell_type>
struct standard_lpfilter_cell_non_atomic
    : public standard_lpfilter_cell_base<cell_type>
{
	using atomic = standard_lpfilter_cell_atomic<cell_type>;
	using base = standard_lpfilter_cell_base<cell_type>;

    using typename base::cstuff;
    using typename base::entry_type;
	using base::base;

	standard_lpfilter_cell_non_atomic() = default;
	standard_lpfilter_cell_non_atomic(
        const standard_lpfilter_cell_non_atomic&) = default;
	standard_lpfilter_cell_non_atomic(const atomic& other)
        : base(other.cell.load())
    {}

	standard_lpfilter_cell_non_atomic& operator=(
        const standard_lpfilter_cell_non_atomic&) = default;

	standard_lpfilter_cell_non_atomic& operator=(const atomic& other)
	{
		this->cell = other.cell.load();
		return *this;
	}

	// setter

	inline void set_remainder(const cstuff& cc,
                              entry_type remainder,
                              size_t pos)
	{
		this->cell &= ~(cc.remainder_mask << base::shift_amount(cc, pos));
		this->cell |= static_cast<cell_type>(remainder)
                        << base::shift_amount(cc, pos);
	}
};

template<class cell_type>
struct standard_lpfilter_cell_atomic
    : standard_lpfilter_cell_base<std::atomic<cell_type>, cell_type>
{
	using non_atomic = standard_lpfilter_cell_non_atomic<cell_type>;
	using base = standard_lpfilter_cell_base<std::atomic<cell_type>, cell_type>;

    using typename base::cstuff;
    using typename base::entry_type;
	using base::base;

	standard_lpfilter_cell_atomic() = default;
	standard_lpfilter_cell_atomic(const standard_lpfilter_cell_atomic& other)
        : base(other.cell.load()) {}

	standard_lpfilter_cell_atomic(const non_atomic& other) : base(other.cell) {}

	inline standard_lpfilter_cell_atomic& operator=(
        const standard_lpfilter_cell_atomic& other)
	{
		this->cell.store(other.cell.load());
		return *this;
	}

	inline standard_lpfilter_cell_atomic& operator=(const non_atomic& other)
	{
		this->cell.store(other.cell);
		return *this;
	}

	inline operator non_atomic()
	{
		return { this->cell.load() };
	}

	// setter

	inline void set_remainder(const cstuff& cc,
                              entry_type remainder,
                              size_t pos)
	{
		update_cell(static_cast<cell_type>(remainder)
                      << base::shift_amount(cc, pos),
                    cc.remainder_mask << base::shift_amount(cc, pos));
	}

	inline void update_cell(const cell_type& val, const cell_type& mask)
	{
		cell_type cell;
		cell_type new_cell;

		do
		{
			cell = this->cell.load();
			new_cell = (cell & ~mask) | val;
		} while (!CAS(cell, new_cell));
	}

    inline bool CAS(cell_type& expected, const cell_type& desired)
	{
		return this->cell.compare_exchange_strong(expected, desired);
	}

	inline bool CAS(non_atomic& expected, const non_atomic& desired)
	{
		return CAS(expected.cell, desired.cell);
	}

};
} // namespace qf
