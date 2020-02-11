#pragma once

#include <array>
#include <atomic>
//#include "benchmarks/filter_select.h"
#include "utils/utilities.h"
#include "utils/quotient_filter_type_traits.h"
//#include "nt_quotient_filter.h"
#include "implementation/reclamation_strategies/seq_reclamation.hpp"

#ifdef HAZARD
#include "implementation/reclamation_strategies/hazard_reclamation.hpp"
using default_reclamation_strategy = qf::hazard::hazard_reclamation_strategy;

#else
#include "implementation/reclamation_strategies/delayed_reclamation.hpp"
using default_reclamation_strategy = qf::delayed::delayed_reclamation_strategy;
#endif



namespace qf
{

class concurrent_growing_handler
{
public:
	enum class growing_task { master, worker, old_data };

	class growing_status_data
	{
	private:
		size_t data;
	public:
		growing_status_data(size_t data = 0) noexcept : data(data) {};

		size_t index() const { return data & index_mask; }
		bool is_growing() const {return static_cast<bool>(data & growing_mask);}
		bool can_grow() const { return index() < max_index; }

		void set_growing(bool value)
		{
			data &= index_mask;
			data |= growing_mask * static_cast<size_t>(value);
		}
	};

private:
	std::atomic<growing_status_data> atomic_data;

	static constexpr size_t max_index = REMAINDER_COUNT - 1;
	static constexpr size_t growing_mask = ONE << ((sizeof(size_t) * 8) - 1);
	static constexpr size_t index_mask = ~growing_mask;

public:
	concurrent_growing_handler(size_t index)
        : atomic_data(index & index_mask) {}
	concurrent_growing_handler(const concurrent_growing_handler& other)
        : atomic_data(other.atomic_data.load()) {}

	concurrent_growing_handler& operator=(const concurrent_growing_handler& other)
	{
		atomic_data = other.atomic_data.load();
		return *this;
	}

	concurrent_growing_handler(concurrent_growing_handler&&) = default;
	concurrent_growing_handler& operator=(concurrent_growing_handler&&) = default;

	growing_status_data status() const { return atomic_data.load(); };

	growing_task start_growing(growing_status_data status_data)
	{
		if (status_data.is_growing())
		{
			if (status_data.index() != atomic_data.load().index())
				return growing_task::old_data;

			return growing_task::worker;
		}

		auto growing = status_data;
		growing.set_growing(true);

		if (atomic_data.compare_exchange_strong(status_data, growing))
		{
			return growing_task::master;
		}
		else if (status_data.index() == growing.index())
		{
			return growing_task::worker;
		}

		return growing_task::old_data;
	}

	// should only be called once by master
	void end_growing()
	{
		atomic_data = atomic_data.load().index() + 1;
	}
};





template <class growing_filter_base>
class non_templated_filter_creator
{
public:
    static void create_filter(growing_filter_base& ref, size_t current_index)
    {
        using quotient_filter = typename growing_filter_base::quotient_filter;

        ref.filters[current_index].store(
            new quotient_filter(ONE << ref.quotient_bits,
                                ref.index_to_reamainder_bits(current_index),
                                ref.hf) );
    }

    non_templated_filter_creator() = delete;
};



template <class growing_filter_base>
class templated_filter_creator
{
public:

    template <size_t rbits>
    using quotient_filter_template =
        typename growing_filter_base::quotient_filter::template instanciated<rbits>;

    template <size_t... Indices>
    static void
    create_filter(growing_filter_base& ref, size_t current_index,
                  std::index_sequence<Indices...>)
    {

        // only necessary for the templated
        ((current_index == Indices
          && (ref.filters[Indices] =
              new quotient_filter_template
              <growing_filter_base::index_to_reamainder_bits(Indices)>
              ( ONE << ref.quotient_bits, ref.hf),
              true)
            )
         || ...);

    }

    static void create_filter(growing_filter_base& ref, size_t current_index)
    {
        create_filter(ref, current_index, REMAINDER_SEQUENCE);
    }

    templated_filter_creator() = delete;
};




template <class base_filter,
          class reclamation_strategy = default_reclamation_strategy>
class growing_quotient_filter_base
{
private:
    static_assert(base_filter::is_growing_compatible,
                  "base table in growing_quotient_filter_base non compatible");
    using this_type       = growing_quotient_filter_base<base_filter>;
    using derived_type    = base_filter;
    using quotient_filter = base_filter;


    using rec_strat      = typename std::conditional
                                < quotient_filter::is_sequential,
                                  qf::seq_reclamation::seq_reclamation_strategy,
                                  reclamation_strategy
                                >::type;

    using rec_pointer = typename rec_strat::
                            template atomic_reclamation <quotient_filter>;
    using rec_handle  = typename rec_strat::
                            template reclamation_handle <quotient_filter>;
    using rec_manager = typename rec_strat::
                            template reclamation_manager<quotient_filter>;


    using filter_creator = typename std::conditional
                                < quotient_filter::is_templated,
                                  templated_filter_creator<this_type>,
                                  non_templated_filter_creator<this_type>
                                >::type;
    friend filter_creator;


    using growing_handler_type = concurrent_growing_handler;
    using growing_status_data  = typename growing_handler_type::growing_status_data;

public:
    static constexpr bool is_growing_compatible = false;
    static constexpr bool is_templated  = quotient_filter::is_templated;
    static constexpr bool is_sequential = quotient_filter::is_sequential;
    static constexpr bool is_growing    = true;
    static constexpr bool is_dynamic    = false;
    static constexpr bool uses_handle   = true;

    using key_type           = typename quotient_filter::key_type;
    using hash_function_type = typename quotient_filter::hash_function_type;
    using hashed_type        = typename quotient_filter::hashed_type;

    class handle_type
    {
    private:
        friend growing_quotient_filter_base;
        growing_quotient_filter_base& qf;
        rec_handle handle;

    public:
        handle_type(growing_quotient_filter_base& qf, rec_handle&& handle);
        handle_type(const handle_type& ) = delete;
        handle_type& operator =(const handle_type& ) = delete;
        handle_type(handle_type&& other);
        handle_type& operator =(handle_type&& other);

        inline bool insert(const key_type& key);
        inline bool insert_hash(const hashed_type& hashed);
        inline bool quick_insert(const key_type& key);
        inline bool quick_insert_hash(const hashed_type& hashed);
        inline qf::ContainsResult contains(const key_type& key);
        inline qf::ContainsResult contains_hash(const hashed_type& hashed);
        inline qf::ContainsResult unsafe_contains(const key_type& key);
        inline qf::ContainsResult unsafe_contains_hash(const hashed_type& hashed);
        inline size_t capacity() const;
    };
    friend handle_type;

private:
    std::atomic_bool workers_start = false;
    size_t quotient_bits;
    size_t max_index;
    size_t cluster_length_limit;
    size_t block_size;

    hash_function_type    hf;
    growing_handler_type  growing_handler;
    rec_manager           filter_manager;
    std::array<rec_pointer, REMAINDER_COUNT> filters{};

public:
    growing_quotient_filter_base(size_t min_capacity = DEFAULT_MIN_CAPACITY,
                    size_t remainder_bits = DEFAULT_REMAINDER_BITS,
                    const hash_function_type& hf = hash_function_type());

    growing_quotient_filter_base(size_t quotient_bits,
                    size_t remainder_bits,
                    size_t min_remainder_bits,
                    const hash_function_type& hf = hash_function_type());

    ~growing_quotient_filter_base();

    growing_quotient_filter_base(const growing_quotient_filter_base&) = delete;
    growing_quotient_filter_base& operator=(const growing_quotient_filter_base&) = delete;

    growing_quotient_filter_base(growing_quotient_filter_base&&);
    growing_quotient_filter_base& operator=(growing_quotient_filter_base&&);


    inline bool insert_hash      (const hashed_type& hashed,
                                  rec_handle& handle);
    inline qf::ContainsResult contains_hash(const hashed_type& hashed,
                                  rec_handle& handle);
    inline bool quick_insert_hash(const hashed_type& hashed);
    inline qf::ContainsResult unsafe_contains_hash(const hashed_type& hashed);


    inline bool insert                (const key_type& key, rec_handle& handle)
    { return insert_hash(hf(key), handle); }
    inline qf::ContainsResult contains(const key_type& key, rec_handle& handle)
    { return contains_hash(hf(key), handle); }
    inline bool quick_insert          (const key_type& key)
    { return quick_insert(hf(key)); }
    inline qf::ContainsResult unsafe_contains(const key_type& key)
    { return unsafe_contains_hash(hf(key)); }


    bool grow(growing_status_data& status, rec_handle& handle);


    inline size_t capacity() const
    {
        auto index = growing_handler.status().index();
        return filters[index].load()->capacity();
    }

    // size_t memory_usage_bytes() const;
    // size_t unused_memory_bits() const;
    // double fill_level() const;

    handle_type get_handle()
    { return handle_type(*this, get_handle_internal()); }

private:
    static constexpr size_t remainder_bits_to_index(size_t remainder_bits)
    { return REMAINDER_COUNT - remainder_bits + MIN_REMAINDER_BITS - 1; }

    static constexpr size_t index_to_reamainder_bits(size_t index)
    { return REMAINDER_COUNT - index + MIN_REMAINDER_BITS - 1; }

    rec_handle get_handle_internal()
    { return filter_manager.get_handle(); }

    void create_next_filter(size_t current_filter);
    void create_filter(size_t current_filter);
};





// GROWING_QUOTIENT_FILTER_BASE:   CONSTRUCTORS *********************************************

template <class base_filter, class rec_strat>
growing_quotient_filter_base<base_filter, rec_strat>::
growing_quotient_filter_base(size_t min_capacity,
                             size_t remainder_bits,
                             const hash_function_type& hf)
    : max_index(filters.size() -1),
      hf(hf),
      growing_handler(remainder_bits_to_index(remainder_bits))
{
    quotient_bits = static_cast<size_t>(std::ceil(std::log2(min_capacity)));
    auto index = remainder_bits_to_index(remainder_bits);

    create_filter(index);

    auto capacity        = filters[index].load()->capacity();
    cluster_length_limit = std::min(MAX_RUN_SIZE, capacity/2);
    block_size           = std::min(capacity+shift_buffer, GROWING_BLOCK_SIZE);

}

template <class base_filter, class rec_strat>
growing_quotient_filter_base<base_filter, rec_strat>::
growing_quotient_filter_base(size_t quotient_bits,
                             size_t remainder_bits,
                             size_t min_remainder_bits,
                             const hash_function_type& hf)
    : quotient_bits(quotient_bits),
      max_index(remainder_bits_to_index(min_remainder_bits)),
      hf(hf),
      growing_handler(remainder_bits_to_index(remainder_bits))
{
    auto index = remainder_bits_to_index(remainder_bits);

    create_filter(index);

    auto capacity = filters[index].load()->capacity();
    cluster_length_limit = std::min(MAX_RUN_SIZE, capacity/2);
    block_size = std::min(capacity + shift_buffer, GROWING_BLOCK_SIZE);
}


template <class base_filter, class rec_strat>
growing_quotient_filter_base<base_filter, rec_strat>::
growing_quotient_filter_base(growing_quotient_filter_base&& source)
    : workers_start(false),
      quotient_bits(source.quotient_bits),
      max_index(source.max_index),
      cluster_length_limit(source.cluster_length_limit),
      block_size(source.block_size),
      hf(source.hf),
      growing_handler(std::move(source.growing_handler)),
      filter_manager()
{
    for (size_t i = 0; i < REMAINDER_COUNT; ++i)
    {
        filters[i].store(source.filters[i].load());
        source.filters[i].store(nullptr);
    }
}


template <class base_filter, class rec_strat>
growing_quotient_filter_base<base_filter, rec_strat>&
growing_quotient_filter_base<base_filter, rec_strat>::
operator=(growing_quotient_filter_base&& source)
{
    if (this == &source) return *this;

    this->~this_type();
    new (this) this_type(std::move(source));
    return *this;
}


template <class base_filter, class rec_strat>
growing_quotient_filter_base<base_filter, rec_strat>::
~growing_quotient_filter_base()
{
    for (auto& filter : filters)
    {
        auto ptr = filter.load();
        if (ptr) { delete ptr; }
    }
}




// GROWING_QUOTIENT_FILTER_BASE:   MAIN FUNCTIONALITY ***************************************

template <class base_filter, class rec_strat>
bool
growing_quotient_filter_base<base_filter, rec_strat>::
insert_hash(const hashed_type& hashed, rec_handle& handle)
{
    // Special case for the sequential implementation
    if constexpr (is_sequential)
    {
        auto current_status = growing_handler.status();
        auto result = filters[current_status.index()].load()->insert_hash(hashed);

		if (!result.successful())
			return grow(current_status, handle) && insert_hash(hashed, handle);

		if (result.cluster_length() > cluster_length_limit)
			return grow(current_status, handle);
        //FIX: possibly returns false even though insert was successful!

		return true;
    }

    // Concurrent implementation
    while (true)
    {
        auto before_status = growing_handler.status();
        while (before_status.is_growing())
        {
            grow(before_status, handle);
            before_status = growing_handler.status();
        }

        auto qf_ptr = handle.protect(filters[before_status.index()]);
        if (qf_ptr == nullptr)
            continue;

        auto result = qf_ptr->insert_hash(hashed);
        handle.unprotect_last();

        auto after_status = growing_handler.status();

        while (after_status.is_growing())
        {
            grow(after_status, handle);
            after_status = growing_handler.status();
        }

        if (before_status.index() == after_status.index())
        {
            if (result.successful())
            {
                if (result.cluster_length() > cluster_length_limit)
                    return grow(after_status, handle);
                return true;
            }
            else if (!grow(after_status, handle))
            {
                return false;
            }
        }
    }

    return false;
}


template <class base_filter, class rec_strat>
qf::ContainsResult
growing_quotient_filter_base<base_filter, rec_strat>::
contains_hash(const hashed_type& hashed, rec_handle& handle)
{
    if constexpr (is_sequential)
    {
        auto current_status = growing_handler.status();
        // std::cout << "contains start  cindex:" << current_status.index() << std::flush;
        auto result = filters[current_status.index()].load()->contains_hash(hashed);
        // std::cout << "  result:" << result << std::endl;
        // exit(66);
        return result;
    }

    while (true)
    {
        auto current_status = growing_handler.status();
        while (current_status.is_growing())
        {
            grow(current_status, handle);
            current_status = growing_handler.status();
        }


        auto qf_ptr = handle.protect(filters[current_status.index()]);
        if (qf_ptr == nullptr)
            continue;

        auto result = qf_ptr->contains_hash(hashed);
        handle.unprotect_last();
        return result;
    }
}


template <class base_filter, class rec_strat>
bool
growing_quotient_filter_base<base_filter, rec_strat>::
quick_insert_hash(const hashed_type& hashed)
{
    return filters[max_index].load()->quick_insert_hash(hashed);
}


template <class base_filter, class rec_strat>
qf::ContainsResult
growing_quotient_filter_base<base_filter, rec_strat>::
unsafe_contains_hash(const hashed_type& hashed)
{
    if constexpr (is_sequential)
        return filters[max_index].load()->contains_hash(hashed);
    else
        return filters[max_index].load()->unsafe_contains_hash(hashed);
}



// GROWING AND MULTIPLE FILTER STUFF *******************************************


// bool grow(GrowingStatusData current_status, ReclamationHandle& handle)
// {
//     if (!current_status.can_grow() || current_status.index() == max_index)
//         return false;

//     auto task = growing_handler.start_growing(current_status);

//     while (task == GrowingTask::old_data)
//     {
//         // growing for the old index was already finished successfully
//         return true;
//     }

//     const auto current_index = current_status.index();

//     auto current_qf_ptr = handle.protect(filters[current_index]);
//     if (current_qf_ptr == nullptr)
//         return true;

//     if (task == GrowingTask::master)
//     {

//         static_cast<base_filter*>(this)->create_next_filter(current_index);
//         auto next_qf_ptr = handle.protect(filters[current_index + 1]);

//         current_qf_ptr->set_table_status(TableStatus::growing);
//         workers_start = true;
//         current_qf_ptr->grow(*next_qf_ptr, block_size);
//         workers_start = false;
//         current_qf_ptr->set_table_status(TableStatus::inactive);

//         const auto new_capacity = next_qf_ptr->capacity();
//         cluster_length_limit = std::min(MAX_RUN_SIZE, new_capacity / 2);
//         block_size = std::min(new_capacity + shift_buffer, GROWING_BLOCK_SIZE);
//         quotient_bits++;

//         handle.safe_delete(filters[current_index]);

//         growing_handler.end_growing();
//         handle.unprotect_last(2);
//     }
//     else // worker
//     {
//         while (!workers_start)
//         {
//             if (growing_handler.status().index() > current_index)
//             {
//                 handle.unprotect_last();
//                 return true;  // growing has already finished
//             }
//         };

//         auto next_qf_ptr = handle.protect(filters[current_index + 1]);
//         if (next_qf_ptr == nullptr)
//         {
//             handle.unprotect_last();
//             return true;
//         }

//         current_qf_ptr->grow(*next_qf_ptr, block_size);

//         handle.unprotect_last(2);
//     }

//     return true;
// }

// bool grow()
// {
//     if (current_index >= max_index)
//         return false;

//     static_cast<base_filter*>(this)->create_next_filter(current_index);
//     filters[current_index]->grow(*filters[current_index + 1]);

//     cluster_length_limit = std::min(MAX_RUN_SIZE, filters[current_index + 1]->capacity() / 2);
//     quotient_bits++;
//     current_index++;

//     return true;
// }


template <class base_filter, class rec_strat>
bool
growing_quotient_filter_base<base_filter, rec_strat>::
grow(growing_handler_type::growing_status_data& current_status,
     rec_handle& handle)
{
    // std::cout << "growing start " << std::flush;
    if constexpr (is_sequential)
    {
        size_t current_index = current_status.index();
        // std::cout << "cindex:" << current_index;

        if (current_index >= max_index)
			return false;

	    create_next_filter(current_index);
		filters[current_index].load()->grow(*filters[current_index + 1].load());

		cluster_length_limit = std::min(
                             MAX_RUN_SIZE,
                             filters[current_index + 1].load()->capacity() / 2);
		quotient_bits++;
		growing_handler.end_growing();
        // std::cout << " nindex:" << current_index+1 << std::endl;
        //if (!filters[current_index+1].load()->check_consistency()) exit(66);

		return true;
    }
    else
    {
        if (!current_status.can_grow() || current_status.index() == max_index)
            return false;

        auto task = growing_handler.start_growing(current_status);

        while (task == growing_handler_type::growing_task::old_data)
        {
            // growing for the old index was already finished successfully
            return true;
        }

        const auto current_index = current_status.index();

        auto current_qf_ptr = handle.protect(filters[current_index]);
        if (current_qf_ptr == nullptr)
            return true;

        if (task == growing_handler_type::growing_task::master)
        {

            create_next_filter(current_index);
            auto next_qf_ptr = handle.protect(filters[current_index + 1]);

            current_qf_ptr->set_table_status(TableStatus::growing);
            workers_start = true;
            current_qf_ptr->grow(*next_qf_ptr, block_size);
            workers_start = false;
            current_qf_ptr->set_table_status(TableStatus::inactive);

            const auto new_capacity = next_qf_ptr->capacity();
            cluster_length_limit = std::min(MAX_RUN_SIZE, new_capacity / 2);
            block_size = std::min(new_capacity + shift_buffer, GROWING_BLOCK_SIZE);
            quotient_bits++;

            handle.safe_delete(filters[current_index]);

            growing_handler.end_growing();

            // if (!next_qf_ptr->check_consistency()) exit(66);
            // else std::cout << "growK" << std::endl;

            handle.unprotect_last(2);
        }
        else // worker
        {
            while (!workers_start)
            {
                if (growing_handler.status().index() > current_index)
                {
                    handle.unprotect_last();
                    return true;  // growing has already finished
                }
            };

            auto next_qf_ptr = handle.protect(filters[current_index + 1]);
            if (next_qf_ptr == nullptr)
            {
                handle.unprotect_last();
                return true;
            }

            current_qf_ptr->grow(*next_qf_ptr, block_size);

            handle.unprotect_last(2);
        }

        return true;
    }
}


template <class base_filter, class rec_strat>
void
growing_quotient_filter_base<base_filter, rec_strat>::
create_next_filter(size_t current_index)
{
    // should work in templated or nontemplated scenarios
    auto ptr = filters[current_index].load();
    filters[current_index+1].store(ptr->create_bigger_QF(this->hf));
}

// template <class base_filter, class rec_strat> template<size_t... Indices>
// void
// growing_quotient_filter_base<base_filter, rec_strat>::
// create_filter(size_t current_index, std::index_sequence<Indices...>)
// {
//     // only necessary for the templated
//     if constexpr (is_templated)
//     {
//         ((current_index == Indices
// 		  && (this->filters[Indices] =
// 			  new quotient_filter<index_to_reamainder_bits(Indices)>(
//                   ONE << this->quotient_bits,
//                   this->hf),
// 			  true)
// 		 )
// 		 || ...);
//     }
// }

template <class base_filter, class rec_strat>
void
growing_quotient_filter_base<base_filter, rec_strat>::
create_filter(size_t current_index)
{
    // if constexpr (!is_templated)
    // {
    //     this->filters[current_index] =
    //         new quotient_filter(ONE << this->quotient_bits,
    //                             index_to_reamainder_bits(current_index),
    //                             this->hf);
    // }
    // else
    // {
    //     create_filter(current_index, REMAINDER_SEQUENCE);
    // }
    filter_creator::create_filter(*this, current_index);
}









// Handle Stuff ****************************************************************
template <class base_filter, class rec_strat>
growing_quotient_filter_base<base_filter, rec_strat>::handle_type::handle_type(
    growing_quotient_filter_base& qf, rec_handle&& handle)
    : qf(qf), handle(std::move(handle))
{ }

template <class base_filter, class rec_strat>
growing_quotient_filter_base<base_filter, rec_strat>::handle_type::handle_type(
    handle_type&& other)
    : qf(other.qf), handle(std::move(other.handle))
{ }

template <class base_filter, class rec_strat>
typename growing_quotient_filter_base<base_filter, rec_strat>::handle_type&
growing_quotient_filter_base<base_filter, rec_strat>::handle_type::operator=(
    handle_type&& other)
{
    if (this == &other) return *this;

    this->~handle_type();
    new (this) handle_type(std::move(other));
    return *this;
}

template <class base_filter, class rec_strat>
bool
growing_quotient_filter_base<base_filter, rec_strat>::
handle_type::insert(const key_type& key)
{ return qf.insert(key, handle); }

template <class base_filter, class rec_strat>
bool
growing_quotient_filter_base<base_filter, rec_strat>::handle_type::
insert_hash(const hashed_type& hashed)
{ return qf.insert_hash(hashed, handle); }

template <class base_filter, class rec_strat>
bool
growing_quotient_filter_base<base_filter, rec_strat>::handle_type::
quick_insert(const key_type& key)
{ return qf.quick_insert(key); }

template <class base_filter, class rec_strat>
bool
growing_quotient_filter_base<base_filter, rec_strat>::handle_type::
quick_insert_hash(const hashed_type& hashed)
{ return qf.quick_insert_hash(hashed); }

template <class base_filter, class rec_strat>
qf::ContainsResult
growing_quotient_filter_base<base_filter, rec_strat>::handle_type::
contains(const key_type& key)
{ return qf.contains(key, handle); }

template <class base_filter, class rec_strat>
qf::ContainsResult
growing_quotient_filter_base<base_filter, rec_strat>::handle_type::
contains_hash(const hashed_type& hashed)
{ return qf.contains_hash(hashed, handle); }

template <class base_filter, class rec_strat>
qf::ContainsResult
growing_quotient_filter_base<base_filter, rec_strat>::handle_type::
unsafe_contains(const key_type& key)
{ return qf.unsafe_contains(key); }

template <class base_filter, class rec_strat>
qf::ContainsResult
growing_quotient_filter_base<base_filter, rec_strat>::handle_type::
unsafe_contains_hash(const hashed_type& hashed)
{ return qf.unsafe_contains_hash(hashed); }

template <class base_filter, class rec_strat>
size_t
growing_quotient_filter_base<base_filter, rec_strat>::handle_type::
capacity() const
{ return qf.capacity(); }


template <class base_filter>
using growing_quotient_filter           = growing_quotient_filter_base<
                                              base_filter,
                                              default_reclamation_strategy>;

} // namespace qf
