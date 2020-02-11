#pragma once
#include <atomic>
#include <limits>
#include "utils/definitions.h"

namespace qf {

template<class entry_type, class basic_entry_type = entry_type>
struct nongrouped_cell_base
{
	nongrouped_cell_base() = default;
	nongrouped_cell_base(const entry_type& entry) : entry(entry) {};


	// setter

	void clear()
	{
		entry = 0;
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

	// getter

	bool is_occupied() const
	{
		return entry & occupied_mask;
	}

	bool is_continuation() const
	{
		return entry & continuation_mask;
	}

	bool is_shifted() const
	{
		return entry & shifted_mask;
	}

	bool is_empty() const
	{
		return status() == 0;
	}

	bool is_cluster_start() const
	{
		return status() == cluster_start_status;
	}

	bool is_write_locked() const
	{
		return status() == write_lock_status;
	}

	bool is_read_locked() const
	{
		return status() == read_lock_status;
	}

	bool is_locked() const
	{
		return is_locked(status());
	}

	static bool is_locked(status_type status)
	{
		return (status == read_lock_status) || (status == write_lock_status);
	}

	basic_entry_type remainder() const
	{
		return entry >> status_bits;
	}

	status_type status() const
	{
		return static_cast<status_type>(entry & status_mask);
	}

	static bool status_is_occupied(status_type status)
	{
		return status & occupied_mask;
	}

	static bool status_is_continuation(status_type status)
	{
		return status & continuation_mask;
	}

	static bool status_is_shifted(status_type status)
	{
		return status & shifted_mask;
	}

	// members

	entry_type entry{ 0 };

	static constexpr basic_entry_type one = static_cast<basic_entry_type>(1);
	static constexpr basic_entry_type status_mask = (one << status_bits) - 1;
	static constexpr basic_entry_type occupied_mask = one << 0;
	static constexpr basic_entry_type continuation_mask = one << 1;
	static constexpr basic_entry_type shifted_mask = one << 2;

	static constexpr basic_entry_type cluster_start_status = occupied_mask;
	static constexpr basic_entry_type write_lock_status = continuation_mask;
	static constexpr basic_entry_type read_lock_status = occupied_mask
                                                            | continuation_mask;
};

template<class entry_type>
struct nongrouped_cell_atomic;

template<class entry_type>
struct nongrouped_cell_non_atomic : nongrouped_cell_base<entry_type, entry_type>
{
	using atomic = nongrouped_cell_atomic<entry_type>;
	using base = nongrouped_cell_base<entry_type, entry_type>;

	using base::base;

	nongrouped_cell_non_atomic() = default;
	nongrouped_cell_non_atomic(const nongrouped_cell_non_atomic&) = default;

	nongrouped_cell_non_atomic(const atomic& other)
        : base(other.entry.load())
    {}

	nongrouped_cell_non_atomic& operator=(
        const nongrouped_cell_non_atomic&) = default;
	nongrouped_cell_non_atomic& operator=(const atomic& other)
	{
		this->entry = other.entry.load();
		return *this;
	}

	// setter

	void set_remainder(entry_type remainder)
	{
		this->entry = (this->entry & base::status_mask)
                      | (remainder << status_bits);
	}

	void set_status(status_type status)
	{
		this->entry = (this->entry & ~base::status_mask) | status;
	}

	void set_occupied(bool val)
	{
		this->entry = (this->entry & ~base::occupied_mask)
                      | (base::occupied_mask * static_cast<int>(val));
	}

	void set_continuation(bool val)
	{
		this->entry = (this->entry & ~base::continuation_mask)
                      | (base::continuation_mask * static_cast<int>(val));
	}

	void set_shifted(bool val)
	{
		this->entry = (this->entry & ~base::shifted_mask)
                      | (base::shifted_mask * static_cast<int>(val));
	}

	void set_write_lock()
	{
		set_status(base::write_lock_status);
	}

	void set_read_lock()
	{
		set_status(base::read_lock_status);
	}

	void or_status(status_type status)
	{
		this->entry |= status;
	}
};

template<class entry_type>
struct nongrouped_cell_atomic
    : nongrouped_cell_base<std::atomic<entry_type>, entry_type>
{
	using non_atomic = nongrouped_cell_non_atomic<entry_type>;
	using base = nongrouped_cell_base<std::atomic<entry_type>, entry_type>;

	using base::base;

	nongrouped_cell_atomic() = default;
	nongrouped_cell_atomic(const nongrouped_cell_atomic& other)
        : base(other.entry.load())
    {}

	nongrouped_cell_atomic(const non_atomic& other) : base(other.entry) {}
	nongrouped_cell_atomic& operator=(const nongrouped_cell_atomic& other)
	{
		this->entry.store(other.entry.load());
		return *this;
	}

	nongrouped_cell_atomic& operator=(const non_atomic& other)
	{
		this->entry.store(other.entry);
		return *this;
	}

	operator non_atomic()
	{
		return { this->entry.load() };
	}

	void set_remainder(entry_type remainder)
	{
		update_entry(remainder << status_bits, ~base::status_mask);
	}

	void set_status(status_type status)
	{
		update_entry(static_cast<entry_type>(status), base::status_mask);
	}

	void set_continuation(bool val)
	{
		update_entry(base::continuation_mask * static_cast<int>(val),
                     base::continuation_mask);
	}

	void set_shifted(bool val)
	{
		update_entry(base::shifted_mask * static_cast<int>(val),
                     base::shifted_mask);
	}

	void set_occupied(bool val)
	{
		update_entry(base::occupied_mask * static_cast<int>(val),
                     base::occupied_mask);
	}

	void set_read_lock()
	{
		set_status(base::read_lock_status);
	}

	void set_write_lock()
	{
		set_status(base::write_lock_status);
	}

	void update_entry(const entry_type& val, const entry_type& mask)
	{
		entry_type local_entry;
		entry_type new_entry;

		do
		{
			local_entry = this->entry.load();
			new_entry = (local_entry & ~mask) | val;
		} while (!CAS(local_entry, new_entry));
	}

	bool CAS(entry_type& expected, const entry_type& desired)
	{
		return this->entry.compare_exchange_strong(expected, desired);
	}

	bool CAS(non_atomic& expected, const non_atomic& desired)
	{
		return CAS(expected.entry, desired.entry);
	}
};

template<class entry_type>
struct nongrouped_cell_lock_guard
{
	nongrouped_cell_lock_guard() = default;
	nongrouped_cell_lock_guard(size_t pos,
                               nongrouped_cell_atomic<entry_type>& locked_slot,
                               status_type status)
        : pos(pos), slot(&locked_slot), status(status)
    {};
	~nongrouped_cell_lock_guard() { release(); }

	nongrouped_cell_lock_guard(const nongrouped_cell_lock_guard&) = delete;
	nongrouped_cell_lock_guard& operator= (
        const nongrouped_cell_lock_guard&) = delete;

	nongrouped_cell_lock_guard(nongrouped_cell_lock_guard&& other)
        : pos(std::move(other.pos)),
          slot(std::move(other.slot)),
          status(std::move(other.status))
	{
		other.clear();
	};

	nongrouped_cell_lock_guard& operator= (nongrouped_cell_lock_guard&& other)
	{
		if (pos != other.pos)
			release();

		pos = std::move(other.pos);
		slot = std::move(other.slot);
		status = std::move(other.status);

		other.clear();
		return *this;
	}

	void set_data(size_t pos,
                  nongrouped_cell_atomic<entry_type>& locked_slot,
                  status_type status)
	{
		if (this->pos != pos)
		{
			release();
			this->pos = pos;
		}

		this->slot = &locked_slot;
		this->status = status;
	}

	bool locked() const { return slot != nullptr; }
	void clear() { slot = nullptr; }

	void release()
	{
		if (slot)
		{
			slot->set_status(status);
			slot = nullptr;
		}
	}

	size_t pos = std::numeric_limits<size_t>::max();
	nongrouped_cell_atomic<entry_type>* slot = nullptr;
	status_type status{};
};

template<class entry_type>
struct nongrouped_read_lock_guard : nongrouped_cell_lock_guard<entry_type> {};

template<class entry_type>
struct nongrouped_write_lock_guard : nongrouped_cell_lock_guard<entry_type> {};


} // namespace qf
