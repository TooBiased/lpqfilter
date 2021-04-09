#pragma once
/*******************************************************************************
 * implementation/base_filter/nongrouped/nongrouped_qfilter_locking.hpp
 *
 * external lock based concurrent quotient filter with non-grouped slots
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/



#include <mutex>

#include "implementation/utilities.hpp"
#include "implementation/definitions.hpp"
#include "nongrouped_qfilter_seq.hpp"


namespace qf
{

template <class Key, class Hash>
	class nongrouped_qfilter_locking
{
public:
	using key_type = Key;
	using hash_function_type = Hash;

private:
	static constexpr size_t lock_count = 1025;

	Hash                                  hf;
	QuotientFilterSeq<Key, Hash>          filter;
	mutable std::unique_ptr<std::mutex[]> locks = std::make_unique<std::mutex[]>
        (lock_count);
	size_t                                lock_range;


public:

	nongrouped_qfilter_locking(size_t min_capacity = DEFAULT_MIN_CAPACITY,
                            size_t remainder_bits = DEFAULT_REMAINDER_BITS,
                            const hash_function_type& hf = hash_function_type())
		: hf(hf), filter(min_capacity, remainder_bits, hf)
	{
		constexpr std::uint64_t ONE = 1;
		const std::uint64_t capacity =
            ONE << capacity_to_quotient_bits(min_capacity);
		lock_range = std::max(ONE, capacity / (lock_count - 1));
	}

	nongrouped_qfilter_locking(const nongrouped_qfilter_locking& other)
		: hf(other.hf), filter(other.filter), lock_range(other.lock_range)
	{}

	nongrouped_qfilter_locking& operator=(
        const nongrouped_qfilter_locking& other)
	{
		hf = other.hf;
		filter = other.filter;
		lock_range = other.lock_range;

		return *this;
	}

	nongrouped_qfilter_locking(nongrouped_qfilter_locking&&) = default;
	nongrouped_qfilter_locking& operator=(nongrouped_qfilter_locking&&)=default;

    bool insert(const key_type& key)
	{
		const auto hashed = hf(key);
		const auto qr = filter.get_quotient_and_remainder(hashed);

		std::lock_guard<std::mutex> lock1(locks[qr.first / lock_range]);
		std::lock_guard<std::mutex> lock2(locks[qr.first / lock_range + 1]);
		return filter.insert_hash(hashed).successful();
	}

	bool contains(const key_type& key) const
	{
		const auto hashed = hf(key);
		const auto qr = filter.get_quotient_and_remainder(hashed);

		std::lock_guard<std::mutex> lock1(locks[qr.first / lock_range]);
		std::lock_guard<std::mutex> lock2(locks[qr.first / lock_range + 1]);
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

	template <class Key, class Hash>
	struct isLockingQuotientFilter<nongrouped_qfilter_locking<Key, Hash>>
	{
		static constexpr bool value = true;
	};

} // namespace qf
