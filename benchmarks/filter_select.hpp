#pragma once
/*******************************************************************************
 * benchmarks/filter_select.hpp
 *
 * This is a chooser for that switches between different filter versions at
 * compile time
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/


#include "utils/default_hash.hpp"
namespace htm = utils_tm::hash_tm;

//* NON TEMPLATED COMPACT FILTER VARIANTS **************************************

#if defined(QFILTER)
#include "implementation/base_filter/standard/standard_qfilter_conc.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::standard_qfilter_conc<Key, htm::default_hash>;

#elif defined(QFILTER_SEQ)
#include "implementation/base_filter/standard/standard_qfilter_seq.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::standard_qfilter_seq<Key, htm::default_hash>;

#elif defined(QFILTER_LOCKING)
#include "implementation/base_filter/standard/standard_qfilter_locking.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::standard_qfilter_locking<Key, htm::default_hash>;

#elif defined(LPFILTER)
#include "implementation/base_filter/standard/standard_lpfilter_conc.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::standard_lpfilter_conc<Key, htm::default_hash>;

#elif defined(LPFILTER_SEQ)
#include "implementation/base_filter/standard/standard_lpfilter_seq.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::standard_lpfilter_seq<Key, htm::default_hash>;






//* NONCOMPACT FILTER VARIANTS *************************************************

#elif defined(NONGROUPED_QFILTER)
#include "implementation/base_filter/nongrouped/nongrouped_qfilter_conc.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::nongrouped_qfilter_conc<Key, htm::default_hash>;

#elif defined(NONGROUPED_QFILTER_SEQ)
#include "implementation/base_filter/nongrouped/nongrouped_qfilter_seq.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::nongrouped_qfilter_seq<Key, htm::default_hash>;

#elif defined(NONGROUPED_QFILTER_LOCKING)
#include "implementation/base_filter/nongrouped/nongrouped_qfilter_locking.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::nongrouped_qfilter_locking<Key, htm::default_hash>;

#elif defined(NONGROUPED_LPFILTER)
#include "implementation/base_filter/nongrouped/nongrouped_lpfilter_conc.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::nongrouped_lpqfilter_conc<Key, htm::default_hash>;

#elif defined(NONGROUPED_LPQFILTER_SEQ)
#include "implementation/base_filter/nongrouped/nongrouped_lpfilter_seq.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::nongrouped_lpfilter_seq<Key, htm::default_hash>;






//* TEMPLATED FILTER VARIANTS **************************************************

#elif defined(TEMPLATED_QFILTER)
#include "implementation/base_filter/templated/templated_qfilter_conc.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::templated_qfilter_conc_base<Key, htm::default_hash>;

#elif defined(TEMPLATED_QFILTER_SEQ)
#include "implementation/base_filter/templated/templated_qfilter_seq.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::templated_qfilter_seq_base <Key, htm::default_hash>;

#elif defined(TEMPLATED_QFILTER_LOCKING)
#include "implementation/base_filter/templated/templated_qfilter_locking.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::templated_qfilter_locking<Key, remainder_bits, htm::default_hash>;

#elif defined(TEMPLATED_LPFILTER)
#include "implementation/base_filter/templated/templated_lpfilter_conc.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::templated_lpfilter_conc<Key, remainder_bits, htm::default_hash>;

#elif defined(TEMPLATED_LPFILTER_SEQ)
#include "implementation/base_filter/templated/templated_lpfilter_seq.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = qf::templated_lpfilter_seq<Key, remainder_bits, htm::default_hash>;






//* COMPETITOR VARIANTS ********************************************************

#elif defined(BLOOM)
#include "implementation/base_filter/bloom_filter_conc.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = bloom_filter_conc<Key, remainder_bits, htm::default_hash>;

#elif defined(CLASSIC_BLOOM)
#include "implementation/base_filter/bloom_filter_conc.hpp"
template <class Key>
using QF_type = classic_bloom_filter_conc<Key, htm::default_hash>;


// #elif defined(NOLOSS_QFILTER)
// #include "implementation/base_filter/noloss/noloss_qfilter_wrapper.hpp"
// template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
// using QF_type = noloss_qfilter_wrapper<Key, remainder_bits, htm::default_hash>;
// constexpr QF_concurrency_type QF_concurrency_variant = QF_concurrency_type::concurrent;
// constexpr bool QF_compact = true;

// #elif defined(NOLOSS_TEMPL_QFILTER)
// #include "implementation/base_filter/noloss/noloss_templated_qfilter.hpp"
// template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
// using QF_type = noloss_templated_qfilter<Key, remainder_bits, htm::default_hash>;
// constexpr QF_concurrency_type QF_concurrency_variant = QF_concurrency_type::concurrent;
// constexpr bool QF_compact = true;

#elif defined(CQF_WRAPPER)
#include "implementation/base_filter/cqf_wrapper.hpp"
template <class Key, size_t remainder_bits = qf::DEFAULT_REMAINDER_BITS>
using QF_type = cqf_wrapper<Key, remainder_bits, htm::default_hash>;

#else
#pragma message ( "warning: no qfilter selected" )
#endif // QFILTERs
