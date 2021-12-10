#pragma once
/*******************************************************************************
 * implementation/utilities.hpp
 *
 * some common functions
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/


#include <atomic>
#include <bitset>
#include <cmath>

#include "definitions.hpp"
#include "handle_wrapper.hpp"
#include "utils/default_hash.hpp"

namespace qf
{
struct failed_operation_t
{
};
constexpr failed_operation_t FAILED_OPERATION;

class InsertResult
{
  private:
    static constexpr std::uint64_t failed_op_mask      = ONE << 63;
    static constexpr std::uint64_t cluster_length_mask = ~failed_op_mask;

    std::uint64_t data = failed_op_mask;

  public:
    InsertResult() = default;
    InsertResult(std::uint64_t cluster_length)
        : data(cluster_length & cluster_length_mask)
    {
    }
    InsertResult(failed_operation_t) : data(failed_op_mask) {}

    operator bool() const { return successful(); }

    bool successful() const
    {
        return !static_cast<bool>(data & failed_op_mask);
    }
    std::uint64_t cluster_length() const { return data & cluster_length_mask; }

    void set_cluster_length(std::uint64_t cluster_length)
    {
        data = cluster_length & cluster_length_mask;
    }
    void set_failed_operation() { data = failed_op_mask; }
};



enum class ContainsResultEnum
{
    FAILED_OPERATION,
    ENTRY_NOT_PRESENT,
    ENTRY_PRESENT,
    CANONICAL_SLOT_EMPTY
};

class ContainsResult
{
  private:
    ContainsResultEnum data = ContainsResultEnum::FAILED_OPERATION;

  public:
    ContainsResult() = default;
    ContainsResult(ContainsResultEnum result) : data(result) {}
    ContainsResult(failed_operation_t)
        : data(ContainsResultEnum::FAILED_OPERATION)
    {
    }

    operator bool() const { return data == ContainsResultEnum::ENTRY_PRESENT; }

    ContainsResultEnum result() { return data; }

    void set_result(ContainsResultEnum result) { data = result; }
};



template <class remainder_type>
std::pair<size_t, remainder_type>
get_quotient_and_remainder(std::uint64_t hash, size_t quotient_bits,
                           size_t remainder_bits)
{
    const std::uint64_t fingerprint =
        hash >> (64 - quotient_bits - remainder_bits);
    const size_t         q = fingerprint >> remainder_bits;
    const remainder_type r = fingerprint & ((ONE << remainder_bits) - ONE);
    return {q, r};
}



template <class remainder_type, size_t remainder_bits>
std::pair<size_t, remainder_type>
get_quotient_and_remainder(std::uint64_t hash, size_t quotient_bits)
{
    const std::uint64_t fingerprint =
        hash >> (64 - quotient_bits - remainder_bits);
    const size_t         q = fingerprint >> remainder_bits;
    const remainder_type r = fingerprint & ((ONE << remainder_bits) - ONE);
    return {q, r};
}



// This is the type which holds sequences
template <size_t... Ns> struct Sequence
{
};

// First define the template signature
template <size_t min, size_t... Ns> struct seq_gen;

// Recursion case
template <size_t min, size_t I, size_t... Ns> struct seq_gen<min, I, Ns...>
{
    using type = typename seq_gen<min, I - 1, I, Ns...>::type;
};

// Recursion abort
template <size_t min, size_t... Ns> struct seq_gen<min, min, Ns...>
{
    using type = Sequence<min, Ns...>;
};

template <size_t min, size_t max>
using sequence_t = typename seq_gen<min, max>::type;



size_t capacity_to_quotient_bits(size_t min_capacity)
{
    return static_cast<size_t>(std::ceil(std::log2(min_capacity)));
}

size_t next_power_of_two(size_t n)
{
    return ONE << (static_cast<size_t>(std::ceil(std::log2(n))));
}

size_t filter_capacity(size_t min_capacity)
{
    return next_power_of_two(min_capacity);
}

size_t false_positive_rate_to_remainder_bits(double fp_rate)
{
    return static_cast<size_t>(std::ceil(std::log2(1.0 / fp_rate)));
}

double remainder_bits_to_false_positive_rate(size_t remainder_bits)
{
    return 1. / (ONE << remainder_bits);
}

constexpr size_t int_log2(size_t n)
{
    if (n == 1) return 0;

    size_t count = 0;
    for (; n > 1; count++) { n >>= 1; }

    return count;
}

template <class T> std::string bitstring(T val, size_t bits = sizeof(T) * 8)
{
    std::string str;
    for (T bit = ONE << (bits - 1); bit != 0; bit >>= 1)
    {
        if (val & bit)
            str.push_back('1');
        else
            str.push_back('0');
    }

    return str;
}

template <class QuotientFilter> auto create_handle(QuotientFilter& filter)
{
    if constexpr (QuotientFilter::uses_handle) { return filter.get_handle(); }
    else if constexpr (!QuotientFilter::uses_handle)
    {
        return QF_Handle_Wrapper(filter);
    }
}

} // namespace qf
