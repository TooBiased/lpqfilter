#pragma once

#include <atomic>

using size_t = long unsigned int;

namespace proj
{

namespace mark
{
/* FLAGS & BITMASK DEFINITIONS ************************************************/
    template<size_t i=0>
    inline constexpr size_t flag()
    {
        return (1ull << (63-i));
    }

    template<size_t i=0>
    inline constexpr size_t mask()
    {
        return ~flag<i>();
    }

    template<size_t i=0>
    inline constexpr size_t lower()
    {
        return flag<i>()-1;
    }


/* MARK ***********************************************************************/
    template<size_t i=0, class T=void >
    inline bool atomic_mark(std::atomic<T*>& tar, T*& exp)
    {
        auto temp = (T*)(size_t(exp) | flag<i>());
        return tar.compare_exchange_strong(exp, temp);
    }

    template<size_t i=0, class T=void >
    inline constexpr T* mark(T* ptr)
    {
        return (T*)(size_t(ptr) | flag<i>());
    }


/* UNMARK *********************************************************************/
    template<size_t i=0, class T=void >
    inline bool atomic_unmark_cas(std::atomic<T*>& tar, T* exp)
    {
        return tar.compare_exchange_strong(exp, exp & mask<i>());
    }

    template<size_t i=0, class T=void >
    inline bool atomic_unmark(std::atomic<T*>& tar)
    {
        return tar.fetch_and(mask<i>()) & flag<i>();
    }

    template<size_t i=0, class T=void >
    inline constexpr T* unmark(T* ptr)
    {
        return (T*)(size_t(ptr) & mask<i>());
    }


/* CLEAR **********************************************************************/
    template<size_t i=0, class T=void >
    inline bool atomic_clear(std::atomic<T*>& tar)
    {
        return tar.fetch_and(lower<15>());
    }

    template<class T=void>
    inline constexpr T* clear(T* ptr)
    {
        return (T*)(size_t(ptr) & lower<15>());
    }


/******************************************************************************/
    template<size_t i=0, class T=void>
    inline constexpr bool get_mark(T* ptr)
    {
        return size_t(ptr) & flag<i>();
    }

}; // namespace proj

}; // namespace mark
