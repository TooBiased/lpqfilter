#pragma once
/*******************************************************************************
 * implementation/definitions.hpp
 *
 * This file holds some of the standard definitions of our quotient
 * filter classes.
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/


#include <cstdint>
#include <cstddef>

namespace qf
{

enum class TableStatus : std::uint8_t { active, growing, inactive };

using status_type                = std::uint8_t;
using grouped_cell_base_type     = std::uint64_t;

#ifndef NO_QI
    constexpr bool QF_USE_QUICK_INSERT = true;
#else
    constexpr bool QF_USE_QUICK_INSERT = false;
#endif

constexpr size_t        BASE_INTEGER_BITS = 64;
constexpr std::uint64_t ONE               = 1;
constexpr size_t        status_bits       = 3;
constexpr size_t        shift_buffer      = 100;

constexpr size_t DEFAULT_MIN_CAPACITY    = 1000;
constexpr size_t DEFAULT_REMAINDER_BITS  = 10;
constexpr size_t MIN_REMAINDER_BITS      = 1;
constexpr size_t MAX_REMAINDER_BITS      = BASE_INTEGER_BITS - status_bits;
constexpr size_t REMAINDER_COUNT         = MAX_REMAINDER_BITS - MIN_REMAINDER_BITS + 1;
constexpr size_t MAX_RUN_SIZE            = 50;
constexpr size_t GROWING_BLOCK_SIZE      = 1000000; // replace by table_capacity / "threads" -> get active thread count for table from handle

constexpr auto REMAINDER_SEQUENCE = std::make_index_sequence<REMAINDER_COUNT>{};

} // namespace qf
