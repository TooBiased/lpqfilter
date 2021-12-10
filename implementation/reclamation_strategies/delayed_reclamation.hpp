#pragma once
/*******************************************************************************
 * implementation/reclamation_strategies/delayed_reclamation.hpp
 *
 * This delayed strategy is mostly used for our tests.  It delays freeing of
 * memory, until the handle is destroyed (usually at the end of our tests).
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/



#include <atomic>
#include <vector>

#include "utils/mark_pointer.hpp"

namespace qf::delayed
{

using namespace utils_tm;

template <class T> using atomic_delayed = std::atomic<T*>;


template <class T> class delayed_handle;


template <class T> class delayed_manager
{
  public:
    using handle = delayed_handle<T>;
    friend handle;

  public:
    delayed_manager()                       = default;
    delayed_manager(const delayed_manager&) = default;
    delayed_manager& operator=(const delayed_manager&) = default;

    inline handle get_handle() { return handle(); }
};

template <class T> class delayed_handle
{
  public:
    delayed_handle()                      = default;
    delayed_handle(const delayed_handle&) = delete;
    delayed_handle& operator=(const delayed_handle&) = delete;
    delayed_handle(delayed_handle&& src) { std::swap(src.freelist, freelist); }
    delayed_handle& operator=(delayed_handle&& src)
    {
        std::swap(src.freelist, freelist);
        return *this;
    }
    ~delayed_handle()
    {
        for (auto curr : freelist) delete curr;
    }

    inline T* protect(atomic_delayed<T>& aptr)
    {
        T* read = mark::clear(aptr.load());
        return read;
    }

    inline void unprotect_last(int = 1) { return; }

    inline void unprotect(T*) { return; }

    inline void safe_delete(atomic_delayed<T>& aptr)
    {
        auto old = mark::clear(aptr.exchange(nullptr));
        if (!old) return;
        freelist.push_back(old);
    }

  private:
    std::vector<T*> freelist;
};

class delayed_reclamation_strategy
{
  public:
    template <class T> using atomic_reclamation  = atomic_delayed<T>;
    template <class T> using reclamation_manager = delayed_manager<T>;
    template <class T> using reclamation_handle  = delayed_handle<T>;

    delayed_reclamation_strategy() = delete;
};

}; // namespace qf::delayed
