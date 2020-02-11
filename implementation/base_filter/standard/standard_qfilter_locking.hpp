#pragma once
#include <mutex>

#include "standard_qfilter_seq.hpp"
#include "utils/utilities.h"
#include "utils/definitions.h"


namespace qf
{

template <class Key, class Hash = default_hash>
class standard_qfilter_locking
    : public standard_qfilter_seq<Key, Hash>
{
private:
    using this_type = standard_qfilter_locking<Key, Hash>;
    using base_type = standard_qfilter_seq<Key, Hash>;

public:
    static constexpr bool is_growing_compatible = false;
    static constexpr bool is_templated  = false;
    static constexpr bool is_sequential = false;
    static constexpr bool is_growing    = false;
    static constexpr bool is_dynamic    = false;
    static constexpr bool uses_handle   = false;

    using key_type           = typename base_type::key_type;
    using hash_function_type = typename base_type::hash_function_type;
    using hashed_type        = typename base_type::hashed_type;
    using next_smaller_qf    = this_type;
    using next_bigger_qf     = this_type;

private:

    using base_type::quotient_bits;
    using base_type::table_capacity;
    using base_type::hf;

    static constexpr size_t lock_range     = 4096;
    static constexpr size_t log_lock_range = 12;

    mutable std::unique_ptr<std::mutex[]>  locks;


public:

    standard_qfilter_locking(size_t min_capacity = qf::DEFAULT_MIN_CAPACITY,
                             size_t remainder = qf::DEFAULT_REMAINDER_BITS,
                             const hash_function_type& hf =hash_function_type())
        : base_type(min_capacity, remainder, hf)
    {
        size_t nlocks = (base_type::capacity() >> log_lock_range) + 2;
        //std::cout << "nlocks " <<nlocks << std::endl;
        locks  = std::make_unique<std::mutex[]>(nlocks);
    }

    standard_qfilter_locking(base_type&& base)
        : base_type(base)
    {
        size_t nlocks = (base_type::capacity() >> log_lock_range) + 2;
        //std::cout << "nlocks " <<nlocks << std::endl;
        locks = std::make_unique<std::mutex[]>(nlocks);
    }

    standard_qfilter_locking(const standard_qfilter_locking& other)
        : base_type(other)
    {
        size_t nlocks = (base_type::capacity() >> log_lock_range) + 2;
        //std::cout << "nlocks " <<nlocks << std::endl;
        locks = std::make_unique<std::mutex[]>(nlocks);
    }

    standard_qfilter_locking& operator=(const standard_qfilter_locking& other)
    {
        if (&other == this) return *this;

        this->~this_type();
        new (this) this_type(other);
        return *this;
    }

    standard_qfilter_locking(standard_qfilter_locking&& other)
        : base_type(other)
    {
        size_t nlocks = (base_type::capacity() >> log_lock_range) + 2;
        //std::cout << "nlocks " <<nlocks << std::endl;
        locks = std::make_unique<std::mutex[]>(nlocks);
    }

    standard_qfilter_locking& operator=(standard_qfilter_locking&& other)
    {
        if (&other == this) return *this;
        this->~this_type();
        new (this) this_type(other);
        return *this;
    }



    using base_type::get_quotient_and_remainder;
    using base_type::capacity;



    inline qf::InsertResult insert(const key_type& key) override
    {
        return insert_hash(hf(key));
    }

    inline qf::InsertResult insert_hash(const hashed_type& hashed) override
    {
        auto [qr, re] = get_quotient_and_remainder(hashed);
        qr >>= log_lock_range;

        //std::cout << "qr " << qr << std::endl;

        std::lock_guard<std::mutex> lock0(locks[qr]);
        std::lock_guard<std::mutex> lock1(locks[qr+1]);
        return base_type::insert_hash(hashed);
    }



    inline bool quick_insert(const key_type& key) override
    {
        return quick_insert_hash(hf(key));
    }

    inline bool quick_insert_hash(const hashed_type& hashed) override
    {
        auto [qr, re] = get_quotient_and_remainder(hashed);
        qr >>= log_lock_range;

        // qi only looks at canonical slot >> no need for a second lock
        std::lock_guard<std::mutex> lock0(locks[qr]);
        std::lock_guard<std::mutex> lock1(locks[qr+1]);
        return base_type::quick_insert_hash(hashed);
    }



    inline qf::ContainsResult contains(const key_type& key) const override
    {
        return contains_hash(hf(key));
    }

    inline qf::ContainsResult contains_hash(
        const hashed_type& hashed) const override
    {
        auto [qr, re] = get_quotient_and_remainder(hashed);
        qr >>= log_lock_range;

        std::lock_guard<std::mutex> lock0(locks[qr]);
        std::lock_guard<std::mutex> lock1(locks[qr+1]);
        return base_type::contains_hash(hashed);
    }

    inline qf::ContainsResult unsafe_contains(const key_type& key) const
    {
        return base_type::contains(key);
    }

    inline qf::ContainsResult unsafe_contains_hash(
        const hashed_type& hashed) const
    {
        return base_type::contains_hash(hashed);
    }



    inline this_type* create_bigger_QF(
        const hash_function_type& hash = hash_function_type()) const override
    {
        auto temp = base_type::create_bigger_QF(hash);
        auto result = new this_type(std::move(*temp));
        delete temp;
        return result;
    }

    inline virtual void grow(this_type& target_base)
    {
        base_type::grow(target_base);
    }



    inline size_t memory_usage_bytes() const
    {
        return base_type::memory_usage_bytes()
            + (table_capacity >> log_lock_range) * sizeof(std::mutex);
    }

    inline size_t unused_memory_bits() const
    {
        return base_type::unused_memory_bits()
            + (table_capacity >> log_lock_range) * sizeof(std::mutex) * 8;
    }
};

} // namespace qf
