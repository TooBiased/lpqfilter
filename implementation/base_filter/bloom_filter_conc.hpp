#pragma once
/*******************************************************************************
 * implementation/base_filter/bloom_filter_conc.hpp
 *
 * concurrent Bloom filter implemtentation
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/



#include <algorithm>
#include <atomic>
#include <limits>
#include <memory>

#include "utils/default_hash.hpp"
namespace htm = utils_tm::hash_tm;
#include "implementation/utilities.hpp"
#include "utils/fastrange.hpp"

class atomic_bit_array
{
  public:
    atomic_bit_array(size_t capacity = 1000)
    {
        size_t cap = capacity >> 3;
        if (capacity & 7) cap += 1;

        char* table = static_cast<char*>(operator new(cap));
        std::fill(table, table + cap, 0);
        _table = std::unique_ptr<std::atomic_char[]>(
            reinterpret_cast<std::atomic_char*>(table));
    }

    bool get(size_t i) const
    {
        char temp      = _table[i >> 3].load();
        char selection = 1 << (i & 7);
        return temp & selection;
    }

    bool set(size_t i)
    {
        char temp      = _table[i >> 3].load();
        char selection = 1 << (i & 7);
        while (!(temp & selection))
        {
            char goal = temp | selection;
            if (_table[i >> 3].compare_exchange_weak(temp, goal)) return true;
        }
        return false;
    }

  private:
    std::unique_ptr<std::atomic_char[]> _table;
};


template <class E, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS,
          class hasher = htm::default_hash>
class bloom_filter_conc
{
  public:
    static constexpr bool is_growing_compatible = false;
    static constexpr bool is_templated          = true;
    static constexpr bool is_sequential         = false;
    static constexpr bool is_growing            = false;
    static constexpr bool is_dynamic            = false;
    static constexpr bool uses_handle           = false;
    template <size_t rem_bits>
    using instanciated = bloom_filter_conc<E, rem_bits, hasher>;

    using key_type = E;

    bloom_filter_conc(size_t        min_capacity = qf::DEFAULT_MIN_CAPACITY,
                      const hasher& hf           = hasher())
    {
        _capacity = (1ull << qf::capacity_to_quotient_bits(min_capacity));
        _nbits    = _capacity * (3 + remainder_bits);
        //_factor   = double(cap)/double(std::numeric_limits<size_t>::max());
        _table = atomic_bit_array(_nbits);

        for (size_t i = 0; i < H; ++i)
        {
            _hashfcts[i] = hasher(hf.seed + i * 230498129012931ull);
        }
    }

    bool insert(const E& e)
    {
        for (size_t i = 0; i < H; ++i) { _table.set(hash(i, e)); }
        return true;
    }

    bool contains(const E& e) const
    {
        for (size_t i = 0; i < H; ++i)
        {
            if (!_table.get(hash(i, e))) return false;
        }
        return true;
    }


    size_t capacity() const { return _capacity; }



  private:
    static constexpr size_t get_cap(size_t capacity)
    {
        auto cap = (1ull << qf::capacity_to_quotient_bits(capacity));
        cap *= remainder_bits + 3;
        return cap;
    }

    static constexpr size_t H = std::max<int>(1, (remainder_bits >> 1));
    hasher                  _hashfcts[H];
    // double           _factor;
    size_t           _nbits;
    size_t           _capacity;
    atomic_bit_array _table;

    size_t hash(size_t i, const E& e) const
    {
        return utils_tm::fastrange64(_nbits, _hashfcts[i](e)); // * _factor;
    }
};


template <class E, class hasher = htm::default_hash>
class classic_bloom_filter_conc
{
  public:
    static constexpr bool is_growing_compatible = false;
    static constexpr bool is_templated          = false;
    static constexpr bool is_sequential         = false;
    static constexpr bool is_growing            = false;
    static constexpr bool is_dynamic            = false;
    static constexpr bool uses_handle           = false;

    using key_type = E;

    classic_bloom_filter_conc(
        size_t        min_capacity   = qf::DEFAULT_MIN_CAPACITY,
        size_t        remainder_bits = qf::DEFAULT_REMAINDER_BITS,
        const hasher& hf             = hasher())
    {
        double p = (1 << remainder_bits);

        _nbits =
            double(min_capacity) * std::log(p) / (std::log(2) * std::log(2));
        _nhash = double(_nbits) / double(min_capacity) * std::log(2);

        _capacity = min_capacity;

        _table = atomic_bit_array(_nbits);

        _hashfcts = std::make_unique<hasher[]>(_nhash);
        for (size_t i = 0; i < _nhash; ++i)
        {
            _hashfcts[i] = hasher(hf.seed + i * 230498129012931ull);
        }
    }

    bool insert(const E& e)
    {
        for (size_t i = 0; i < _nhash; ++i) { _table.set(hash(i, e)); }
        return true;
    }

    bool contains(const E& e) const
    {
        for (size_t i = 0; i < _nhash; ++i)
        {
            if (!_table.get(hash(i, e))) return false;
        }
        return true;
    }


    size_t capacity() const { return _capacity; }



  private:
    static constexpr size_t get_cap(size_t capacity)
    {
        // Where do I use this function? is it necessary?
        // in this case we would need the intended fp-rate
        return capacity;
    }

    size_t                    _nbits;
    size_t                    _nhash;
    size_t                    _capacity;
    std::unique_ptr<hasher[]> _hashfcts;
    atomic_bit_array          _table;

    size_t hash(size_t i, const E& e) const
    {
        return utils_tm::fastrange64(_nbits, _hashfcts[i](e)); // * _factor;
    }
};
