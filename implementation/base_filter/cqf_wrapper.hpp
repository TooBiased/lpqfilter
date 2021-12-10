#pragma once
/*******************************************************************************
 * implementation/base_filter/cqf_wrapper.hpp
 *
 * this is a wrapper for the counting quotient filter positioned in misc/cqf/...
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/



#include "utils/default_hash.hpp"
namespace htm = utils_tm::hash_tm;

#include "implementation/definitions.hpp"
#include "implementation/handle_wrapper.hpp"
#include "implementation/utilities.hpp"

extern "C"
{
#include "misc/cqf/include/gqf.h"
#include "misc/cqf/include/gqf_file.h"
#include "misc/cqf/include/gqf_int.h"
}

template <class K, size_t Remainder, class HF = htm::default_hash>
class cqf_wrapper
{
  private:
    static constexpr size_t _remainder = Remainder;
    QF                      _cqf;

  public:
    static constexpr bool is_growing_compatible = false;
    static constexpr bool is_templated          = false;
    static constexpr bool is_sequential         = false;
    static constexpr bool is_growing            = true;
    static constexpr bool is_dynamic            = false;
    static constexpr bool uses_handle           = false;

    using key_type = K;
    using hasher   = HF;

    cqf_wrapper(size_t                         min_capacity = 2048,
                [[maybe_unused]] const hasher& hf           = hasher())
    {
        size_t s   = 2048;
        size_t log = 11;
        while (s < min_capacity)
        {
            s <<= 1;
            ++log;
        }

        if (!qf_malloc(&_cqf, s, log + _remainder, 0, QF_HASH_DEFAULT, 0))
        {
            exit(666);
        }
    }

    cqf_wrapper(const cqf_wrapper&) = default;
    cqf_wrapper& operator=(const cqf_wrapper&) = default;
    cqf_wrapper(cqf_wrapper&&)                 = default;
    cqf_wrapper& operator=(cqf_wrapper&&) = default;

    inline bool insert(const key_type& key)
    {
        int ret = 0;

        // do
        // {
        //     ret = qf_insert(&_cqf, key, 0, 1, QF_TRY_ONCE_LOCK);
        // } while (ret == -2);

        ret = qf_insert(&_cqf, key, 0, 1, QF_WAIT_FOR_LOCK);

        return ret >= 0;
    }

    inline bool contains(const key_type& key)
    {
        size_t count = qf_count_key_value(&_cqf, key, 0, 0);
        return (count > 0);
    }

    inline size_t capacity() const { return _cqf.metadata->nslots; }

    inline size_t memory_usage_bytes() const
    {
        return _cqf.metadata->total_size_in_bytes;
    }

    inline size_t unused_memory_bits() const
    {
        return (_cqf.metadata->total_size_in_bytes << 3) -
               (double(_cqf.metadata->nslots) * (2.125 + _remainder));
    }

  private:
};
