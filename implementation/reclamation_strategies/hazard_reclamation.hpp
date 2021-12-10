#pragma once
/*******************************************************************************
 * implementation/reclamation_strategies/hazard_reclamation.hpp
 *
 * This file implements a hazard reclamation strategy.  When a pointer
 * is freed, its thread checks, weather it is used somewhere.  If it
 * is not used the pointer is freed, otherwise the pointer is marked
 * at the thread who uses it. This thread is now responsible for
 * freeing the pointer.
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/



#include <atomic>
#include <iostream>
#include <tuple>

#include "utils/mark_pointer.hpp"

namespace qf::hazard
{
using namespace utils_tm;

static constexpr size_t MAX_THREADS = 64;
static constexpr size_t MAX_PTR     = 2;
static constexpr size_t DEL_MARK    = 0;

template <class T> using atomic_hazard = std::atomic<T*>;


enum class hazard_op_result
{
    DELETED,
    NOT_DELETED
};

template <class T> class hazard_handle;

template <class T> class hazard_manager
{
  public:
    using handle = hazard_handle<T>;
    friend handle;

  private:
    std::atomic_size_t   next_id;
    std::atomic<handle*> handles[MAX_THREADS];

    inline void safe_delete(atomic_hazard<T>& aptr)
    {
        auto old      = mark::clear(aptr.exchange(nullptr)); //  & PROT_BITMASK
        auto nthreads = next_id.load();

        for (size_t i = 0; i < nthreads; ++i)
        {
            auto h = handles[i].load();
            if (h->delete_check(old)) return;
        }

        delete old;
    }

    inline void process_delete_mark(T* ptr, size_t id)
    {
        auto nthreads = next_id.load();
        for (size_t i = id; i < nthreads; ++i)
        {
            auto h = handles[i].load();
            if (h->delete_check(ptr)) return;
        }

        delete ptr;
    }

  public:
    hazard_manager() : next_id(0)
    {
        for (size_t i = 0; i < MAX_THREADS; ++i) handles[i].store(nullptr);
    }

    hazard_manager(hazard_manager& other) : next_id(other.next_id.load())
    {
        for (size_t i = 0; i < MAX_THREADS; ++i)
            handles[i].store(other.handles[i]);
    }

    hazard_manager(hazard_manager&& other) : next_id(other.next_id.load())
    {
        for (size_t i = 0; i < MAX_THREADS; ++i)
            handles[i].store(other.handles[i]);
    }

    hazard_manager& operator=(hazard_manager& other)
    {
        next_id.store(other.next_id);
        for (size_t i = 0; i < MAX_THREADS; ++i)
            handles[i].store(other.handles[i]);
        return *this;
    }

    hazard_manager& operator=(hazard_manager&& other)
    {
        next_id.store(other.next_id);
        for (size_t i = 0; i < MAX_THREADS; ++i)
            handles[i].store(other.handles[i]);
        return *this;
    }

    inline handle get_handle() { return handle(*this); }
};


template <class T> class hazard_handle
{
  private:
    using manager_type = hazard_manager<T>;
    friend manager_type;

    manager_type&      manager;
    size_t             id;
    std::atomic_size_t nmbr;
    std::atomic<T*>    t[MAX_PTR];

    inline void push(T* ptr)
    {
        auto pos = nmbr.fetch_add(1);
        if (pos > MAX_PTR)
        {
            std::cout << "hzrd overload on thread " << id << std::endl;
            for (size_t i = 0; i < MAX_PTR; i++) { std::cout << t[i] << ", "; }
            std::cout << std::endl;
            exit(2);
        }
        t[pos].store(ptr);
    }

    inline int find(T* ptr)
    {
        for (int i = int(nmbr.load()); i >= 0; --i)
        {
            auto temp = mark::clear(t[i].load());
            if (temp == ptr) return i;
        }
        return -1;
    }

    inline T* pop(T* ptr)
    {
        int pos = find(mark::clear(ptr));
        if (pos < 0)
        {
            std::cout << "unsuccessful pop in hazard pointer not protected!"
                      << std::endl;
            std::cout << "find pointer" << ptr << std::endl;
            for (size_t i = 0; i < MAX_PTR; ++i)
            {
                std::cout << i << ":" << t[i].load() << std::endl;
            }
            exit(3);
        }
        auto lastpos = nmbr.load() - 1;
        auto lastptr = t[lastpos].load();
        auto temp    = t[pos].exchange(lastptr);
        nmbr.fetch_sub(1);
        auto nlastptr = t[lastpos].exchange(nullptr);
        if (nlastptr != lastptr)
        {
            // may need some more complicated unification of marked flags
            // currently only the del flag could be set
            t[pos].store(nlastptr);
            // std::cout << "edge case last is marked during pop in
            // hazard_specific" << std::flush;
        }
        return temp;
    }

    inline T* pop_last()
    {
        size_t pos = nmbr.fetch_sub(1) - 1;
        auto   ptr = t[pos].exchange(nullptr);
        return ptr;
    }

    inline bool check(T* ptr)
    {
        size_t n = nmbr.load();

        for (int i = n - 1; i >= 0; --i)
        {
            auto temp = mark::clear(t[i].load());
            if (ptr == temp) return true;
        }
        return false;
    }

    inline bool delete_check(T* ptr)
    {
        size_t n = nmbr.load();

        for (int i = int(n) - 1; i >= 0; --i)
        {
            auto temp = t[i].load();
            if (ptr == temp) // cannot be marked <-> otherwise double free
            {
                // flag the item for deletion at a later date!
                // protecting[i].compare_exchange_strong(temp, ptr | DEL_FLAG);
                if (mark::atomic_mark<DEL_MARK>(t[i], temp))
                    return true;
                else
                    i++;
            }
        }
        return false;
    }

  public:
    hazard_handle(hazard_manager<T>& manager)
        : manager(manager), id(manager.next_id.fetch_add(1)), nmbr(0)
    {
        for (size_t i = 0; i < MAX_PTR; ++i) t[i].store(nullptr);
        manager.handles[id].store(this);
    }

    hazard_handle(hazard_handle&& other)
        : manager(other.manager), id(other.id), nmbr(other.nmbr.load())
    {
        for (size_t i = 0; i < MAX_PTR; ++i) t[i].store(other.t[i].load());
        manager.handles[id].store(this);
    }

    hazard_handle& operator=(hazard_handle&& other)
    {
        manager = other.manager;
        id      = other.id;
        nmbr    = other.nmbr.load();

        for (size_t i = 0; i < MAX_PTR; ++i) t[i].store(other.t[i].load());

        manager.handles[id].store(this);
        return *this;
    }


    inline T* protect(atomic_hazard<T>& aptr)
    {
        T* read = mark::clear(aptr.load());
        if (read == nullptr) return nullptr;

        push(read);

        T* read2 = mark::clear(aptr.load());
        if (read2 == nullptr) pop_last();

        return read2;
    }

    inline void unprotect_last(int count = 1)
    {
        while (count-- > 0)
        {
            auto temp = pop_last();
            if (mark::get_mark<DEL_MARK>(temp))
                return manager.process_delete_mark(mark::clear(temp), id + 1);
        }
    }

    inline void unprotect(T* ptr)
    {
        auto temp = pop(ptr);

        if (mark::get_mark<DEL_MARK>(temp))
            return manager.process_delete_mark(ptr, id + 1);
    }

    inline void safe_delete(atomic_hazard<T>& aptr)
    {
        manager.safe_delete(aptr);
    }
};

class hazard_reclamation_strategy
{
  public:
    template <class T> using atomic_reclamation  = atomic_hazard<T>;
    template <class T> using reclamation_manager = hazard_manager<T>;
    template <class T> using reclamation_handle  = hazard_handle<T>;

    hazard_reclamation_strategy() = delete;
};

} // namespace qf::hazard
