#pragma once
/*******************************************************************************
 * implementation/composite_filter/dynamic_qfilter.hpp
 *
 * a dynamic quotient filter can grow indefinitely while keeping a
 * false positive bound by allocating new filter levels once one level
 * is full each level consists of a growing filter.  The dynamic
 * quotient filter can be instantiated with many different quotient
 * filter implementations from the basic_filter folder.
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/


#include <array>
#include <atomic>
#include "growing_qfilter.hpp"


namespace qf
{
template<class Dynamic_QF>
class dynamic_quotient_filter_handle
{
private:
    friend Dynamic_QF;
    using handle_type = typename Dynamic_QF::level_handle_type;
    using key_type = typename Dynamic_QF::key_type;

    Dynamic_QF& qf;
    handle_type handle;
    size_t level;

public:
    dynamic_quotient_filter_handle(Dynamic_QF& qf, handle_type&& handle, size_t level)
        : qf(qf), handle(std::move(handle)), level(level) {};

    bool insert(const key_type& key)
    { return qf.insert(key, *this); }

    auto contains(const key_type& key)
    { return qf.contains(key, *this); }

    size_t capacity() const
    { return qf.capacity(); }

    size_t memory_usage_bytes() const
    { return qf.memory_usage_bytes(); }

    size_t unused_memory_bits() const
    {return qf.unused_memory_bits(); }

    double fill_level() const
    { return qf.fill_level(); }
};




#ifdef NO_QI
    static constexpr bool qi_default = false;
#else
    static constexpr bool qi_default = true;
#endif



template <class BaseFilter, bool use_quick_insert = qi_default>
class dynamic_quotient_filter_base
{
    static_assert(BaseFilter::is_growing_compatible,
                  "base table in dynamic_quotient_filter_base non compatible");
    using this_type = dynamic_quotient_filter_base<BaseFilter, use_quick_insert>;
    using base_filter_type  = BaseFilter;
    using growing_filter_type = qf::growing_quotient_filter<BaseFilter>;

    using level_handle_type = typename growing_filter_type::handle_type;

    static constexpr bool quick_insert_enabled = use_quick_insert;

public:
    static constexpr bool is_growing_compatible = false;
    static constexpr bool is_templated  = base_filter_type::is_templated;
    static constexpr bool is_sequential = base_filter_type::is_sequential;
    static constexpr bool is_growing    = true;
    static constexpr bool is_dynamic    = true;
    static constexpr bool uses_handle   = true;

    using key_type           = typename base_filter_type::key_type;
    using hash_function_type = typename base_filter_type::hash_function_type;
    using handle_type = dynamic_quotient_filter_handle<this_type>;
    friend handle_type;


    dynamic_quotient_filter_base(double false_positive_rate = 0.01,
                                 size_t min_capacity = 2ull << 8,
                                 const hash_function_type& hf = hash_function_type());
    dynamic_quotient_filter_base(dynamic_quotient_filter_base&& other);
    dynamic_quotient_filter_base(const dynamic_quotient_filter_base&) = delete;

    dynamic_quotient_filter_base& operator =(dynamic_quotient_filter_base&& other);
    dynamic_quotient_filter_base& operator =(const dynamic_quotient_filter_base&) = delete;

    bool insert(const key_type& key, handle_type& dyn_handle);
    bool contains(const key_type& key, handle_type& dyn_handle);
    // only used for grow_test (should never be called on concurrent table)
    bool insert_grow_test(const key_type& key);
    // only used for grow_test (should never be called on concurrent table)
    void grow();

    size_t capacity() const;
    size_t memory_usage_bytes() const;
    size_t unused_memory_bits() const;
    double fill_level() const;

    handle_type get_handle();

private:
    bool add_level(size_t next_level);
    void update_handle(handle_type& dyn_handle);

    hash_function_type hf;

    size_t lvl0_initial_q_bits;
    size_t lvl0_initial_r_bits;

    using level_ptr_type = std::unique_ptr<growing_filter_type>;
    std::unique_ptr<level_ptr_type[]> levels;
    size_t level_count;

    std::atomic_size_t completed_index = 0;
    std::atomic_size_t in_progress_index = 0;

    static constexpr size_t MIN_INITIAL_QUOTIENT_BITS = 8;
    static constexpr size_t QUOTIENT_BIT_STEP_SIZE = 3;
};




template <class QF, bool qi>
dynamic_quotient_filter_base<QF,qi>::dynamic_quotient_filter_base(
    double false_positive_rate,
    size_t min_capacity,
    const hash_function_type& hf)
    : hf(hf)
{
    const size_t lvl0_min_r_bits = 1 + false_positive_rate_to_remainder_bits(false_positive_rate);
    lvl0_initial_r_bits = lvl0_min_r_bits + QUOTIENT_BIT_STEP_SIZE;
    //std::cout << "starting rbits" << lvl0_min_r_bits << std::endl;
    lvl0_initial_q_bits = std::max(MIN_INITIAL_QUOTIENT_BITS, capacity_to_quotient_bits(min_capacity));

    level_count = (BASE_INTEGER_BITS - lvl0_initial_r_bits - lvl0_initial_q_bits) / (QUOTIENT_BIT_STEP_SIZE + 1);
    levels = std::make_unique<level_ptr_type[]>(level_count);

    levels[0] = std::make_unique<growing_filter_type>(lvl0_initial_q_bits, lvl0_initial_r_bits, lvl0_min_r_bits);
}



template <class QF, bool qi>
dynamic_quotient_filter_base<QF,qi>::dynamic_quotient_filter_base(
    dynamic_quotient_filter_base&& other)
    : hf(other.hf),
      lvl0_initial_q_bits(other.lvl0_initial_q_bits),
      lvl0_initial_r_bits(other.lvl0_initial_r_bits),
      level_count(other.level_count),
      levels(std::move(other.levels)),
      completed_index(other.completed_index.load()),
      in_progress_index(other.in_progress_index.load())
{}



template <class QF, bool qi>
dynamic_quotient_filter_base<QF,qi>&
dynamic_quotient_filter_base<QF,qi>::operator=(dynamic_quotient_filter_base&& other)
{
    hf = other.hf;
    lvl0_initial_q_bits = other.lvl0_initial_q_bits;
    lvl0_initial_r_bits = other.lvl0_initial_r_bits;

    levels = std::move(other.levels);
    level_count = other.level_count;

    completed_index = other.completed_index.load();
    in_progress_index = other.in_progress_index.load();

    return *this;
}



template <class QF, bool qi>
bool
dynamic_quotient_filter_base<QF,qi>::insert(const key_type& key, handle_type& dyn_handle)
{
    const auto hashed = hf(key);
    const size_t local_completed_index = completed_index;

    if constexpr (quick_insert_enabled)
    {
        for (size_t current_level = 0; current_level < local_completed_index; current_level++)
        {
            if (levels[current_level]->quick_insert_hash(hashed))
                return true;
        }
    }

    do
    {
        update_handle(dyn_handle);
        if (dyn_handle.handle.insert_hash(hashed))
            return true;
    } while (add_level(dyn_handle.level + 1));

    return false;
}



template <class QF, bool qi>
bool
dynamic_quotient_filter_base<QF,qi>::insert_grow_test(const key_type& key)
{
    if constexpr (!is_sequential) return true;

    if constexpr (quick_insert_enabled)
                 {
                     for (size_t current_level = 0; current_level < completed_index; current_level++)
                     {
                         if (levels[current_level]->quick_insert(key))
                             return true;
                     }
                 }

    return levels[completed_index]->insert_grow_test(key);
}



template <class QF, bool qi>
void
dynamic_quotient_filter_base<QF,qi>::grow()
{
    if constexpr (!is_sequential) return;

    if (!levels[completed_index]->grow())
        add_level();
}



template <class QF, bool qi>
bool
dynamic_quotient_filter_base<QF,qi>::contains(const key_type& key, handle_type& dyn_handle)
{
    const auto hashed = hf(key);
    const size_t local_completed_index = completed_index;

    for (size_t current_level = 0; current_level < local_completed_index; current_level++)
    {
        const auto contains_result = levels[current_level]->unsafe_contains_hash(hashed).result();

        if (contains_result == qf::ContainsResultEnum::ENTRY_PRESENT)
            return true;

        if constexpr (quick_insert_enabled)
                     {
                         if (contains_result == qf::ContainsResultEnum::CANONICAL_SLOT_EMPTY)
                             return false;
                     }
    }

    update_handle(dyn_handle);
    return dyn_handle.handle.contains_hash(hashed);
}



template <class QF, bool qi>
size_t
dynamic_quotient_filter_base<QF,qi>::capacity() const
{
    return levels[completed_index]->capacity();
}



template <class QF, bool qi>
size_t
dynamic_quotient_filter_base<QF,qi>::memory_usage_bytes() const
{
    size_t mem_usage = 0;

    for (size_t current_level = 0; current_level <= completed_index; current_level++)
    {
        mem_usage += levels[current_level]->memory_usage_bytes();
    }

    return mem_usage;
}



template <class QF, bool qi>
size_t
dynamic_quotient_filter_base<QF,qi>::unused_memory_bits() const
{
    size_t unused_bits = 0;

    for (size_t current_level = 0; current_level <= completed_index; current_level++)
    {
        unused_bits += levels[current_level]->unused_memory_bits();
    }

    return unused_bits;
}



template <class QF, bool qi>
double
dynamic_quotient_filter_base<QF,qi>::fill_level() const
{
    return levels[completed_index]->fill_level();
}



template <class QF, bool qi>
typename dynamic_quotient_filter_base<QF,qi>::handle_type
dynamic_quotient_filter_base<QF,qi>::get_handle()
{
    const size_t local_completed_index = completed_index;
    return handle_type(*this, create_handle(*levels[local_completed_index].get()),
                       local_completed_index);
    return handle_type(*this, create_handle(*levels[local_completed_index].get()),
                       local_completed_index);
}



template <class QF, bool qi>
void
dynamic_quotient_filter_base<QF,qi>::update_handle(handle_type& dyn_handle)
{
    const size_t local_completed_index = completed_index;

    if (dyn_handle.level != local_completed_index)
    {
        dyn_handle.level = local_completed_index;
        dyn_handle.handle = create_handle(*levels[local_completed_index].get());
    }
}



template <class QF, bool qi>
bool
dynamic_quotient_filter_base<QF,qi>::add_level(size_t next_level)
{
    if (next_level <= completed_index)
        return true;

    if (next_level >= level_count)
        return false;

    size_t current_max_level = next_level - 1;

    if (!in_progress_index.compare_exchange_strong(current_max_level, next_level))
    {
        // next level is already being constructed, wait for it to finish
        while (completed_index.load() < next_level) {};
        return true;
    }

    // construct next level

    const size_t initial_q_bits = lvl0_initial_q_bits + next_level * QUOTIENT_BIT_STEP_SIZE;
    const size_t initial_r_bits = lvl0_initial_r_bits + next_level;
    const size_t min_r_bits = initial_r_bits - QUOTIENT_BIT_STEP_SIZE;

    levels[next_level] = std::make_unique<growing_filter_type>(initial_q_bits, initial_r_bits, min_r_bits);
    completed_index = next_level;
    return true;
}



template <class base_filter>
using dynamic_quotient_filter = dynamic_quotient_filter_base<base_filter>;

} // namespace qf
