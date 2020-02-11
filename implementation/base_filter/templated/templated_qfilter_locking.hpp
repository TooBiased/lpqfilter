#pragma once
#include <mutex>

#include "templated_qfilter_seq.hpp"
#include "utils/utilities.h"
#include "utils/definitions.h"


namespace qf
{

template <class Key,
          size_t remainder_bits,
          class Hash = default_hash>
class templated_qfilter_locking
{
public:
    static constexpr bool is_growing_compatible = false;
    static constexpr bool is_templated  = true;
    static constexpr bool is_sequential = false;
    static constexpr bool is_growing    = false;
    static constexpr bool is_dynamic    = false;
    static constexpr bool uses_handle   = false;
    template <size_t rem>
    using instanciated = templated_qfilter_locking<Key, rem, Hash>;

	using key_type = Key;
	using hasher = Hash;
private:

	static constexpr size_t lock_range     = 4096;
	static constexpr size_t log_lock_range = 12;

    size_t                                              nlocks;
	Hash                                                hf;
	templated_qfilter_seq<Key, remainder_bits, Hash>    filter;
	mutable std::unique_ptr<std::mutex[]>               locks;


public:

	templated_qfilter_locking(size_t min_capacity = lock_range,
                              const hasher& hf = hasher())
		: hf(hf), filter(min_capacity, hf)
	{
		constexpr std::uint64_t ONE = 1;
		const std::uint64_t capacity = ONE
                                     << capacity_to_quotient_bits(min_capacity);
        nlocks = (capacity>>log_lock_range) + 2;
        locks  = std::make_unique<std::mutex[]>(nlocks);

	}

	templated_qfilter_locking(const templated_qfilter_locking& other)
		: nlocks(other.nlocks), hf(other.hf), filter(other.filter)
	{
        locks = std::make_unique<std::mutex[]>(nlocks);
    }

    templated_qfilter_locking& operator=(const templated_qfilter_locking& other)
	{
        nlocks = other.nlocks;
		hf     = other.hf;
		filter = other.filter;
        locks  = std::make_unique<std::mutex[]>(nlocks);

		return *this;
	}

	templated_qfilter_locking(templated_qfilter_locking&& other) = default;
	templated_qfilter_locking& operator=(templated_qfilter_locking&&) = default;

	bool insert(const key_type& key)
	{
		const auto hashed = hf(key);
		const auto qr = filter.get_quotient_and_remainder(hashed);

		std::lock_guard<std::mutex> lock1(locks[ qr.first >>log_lock_range]);
		std::lock_guard<std::mutex> lock2(locks[(qr.first >>log_lock_range)+1]);
		return filter.insert_hash(hashed).successful();
	}

	bool contains(const key_type& key) const
	{
		const auto hashed = hf(key);
		const auto qr = filter.get_quotient_and_remainder(hashed);

		std::lock_guard<std::mutex> lock1(locks[ qr.first >>log_lock_range]);
		std::lock_guard<std::mutex> lock2(locks[(qr.first >>log_lock_range)+1]);
		return filter.contains_hash(hashed);
	}

	size_t capacity() const
	{
		return filter.capacity();
	}


	size_t memory_usage_bytes() const
	{
		return filter.memory_usage_bytes() + sizeof(locks);
	}

	size_t unused_memory_bits() const
	{
		return filter.unused_memory_bits() + sizeof(locks) * 8;
	}
};

} // namespace qf
