#pragma once

#include <vector>

namespace qf::seq_reclamation
{

template <class T>
class atomic_seq_reclamation
{
private:
    T* ptr;
public:

    T* load() const    { return ptr; }
    void store(T* ptr) { this->ptr = ptr; }
    atomic_seq_reclamation& operator=(T* new_ptr)
    {
        ptr = new_ptr;
        return *this;
    }
};



template<class T>
class seq_reclamation_handle;


template<class T>
class seq_reclamation_manager
{
public:
    using handle = seq_reclamation_handle<T>;

public:
    seq_reclamation_manager() = default;
    seq_reclamation_manager(const seq_reclamation_manager&) = default;
    seq_reclamation_manager& operator=(const seq_reclamation_manager&) = default;

    handle get_handle() { return handle(); }

};

template<class T>
class seq_reclamation_handle
{
public:
    seq_reclamation_handle()                                         = default;
    seq_reclamation_handle(const seq_reclamation_handle&)            = default;
    seq_reclamation_handle& operator=(const seq_reclamation_handle&) = default;
    seq_reclamation_handle(seq_reclamation_handle&& src)             = default;
    seq_reclamation_handle& operator=(seq_reclamation_handle&& src)  = default;
    ~seq_reclamation_handle() = default;

    T* protect(atomic_seq_reclamation<T>& aptr)
    {
        return aptr.load();
    }

    void unprotect_last(int = 1)
    {
        return;
    }

    void unprotect(T*)
    {
        return;
    }

    void safe_delete(atomic_seq_reclamation<T>& aptr)
    {
        delete aptr.load();
        aptr.store(nullptr);
    }
};



class seq_reclamation_strategy
{
public:
    template <class T>
    using atomic_reclamation = atomic_seq_reclamation<T>;
    template <class T>
    using reclamation_handle = seq_reclamation_handle<T>;
    template <class T>
    using reclamation_manager = seq_reclamation_manager<T>;
};

};
