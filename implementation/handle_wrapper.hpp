#pragma once
/*******************************************************************************
 * implementation/handle_wrapper.hpp
 *
 * Some of our quotient filters use handles, this wrapper enables us
 * to use the same interface for implementations that do not use
 * handles.
 *
 * Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
 *
 * Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
 *
 * All rights reserved. Published under the BSD-2 license in the LICENSE file.
 ******************************************************************************/



#include <cstddef>
#include <utility>

namespace qf
{

template<class QuotientFilter>
class QF_Handle_Wrapper
{
private:
	QuotientFilter& qf;

public:
	using key_type = typename QuotientFilter::key_type;

	QF_Handle_Wrapper(QuotientFilter& filter) : qf(filter) {};

	auto insert(const key_type& key)
	{ return qf.insert(key); }

	auto contains(const key_type& key)
	{ return qf.contains(key); }

	size_t capacity() const
	{ return qf.capacity(); }

	size_t memory_usage_bytes() const
	{ return qf.memory_usage_bytes(); }

	size_t unused_memory_bits() const
	{ return qf.unused_memory_bits(); }

	double fill_level() const
	{ return qf.fill_level(); }
};



template <typename QuotientFilter, bool UseHazard>
struct HandleTypeHelper;

template<class QuotientFilter>
struct HandleTypeHelper<QuotientFilter, false>
{
	using type = QF_Handle_Wrapper<QuotientFilter>;
};

template<class QuotientFilter>
struct HandleTypeHelper<QuotientFilter, true>
{
	using type = decltype(std::declval<QuotientFilter>().get_handle());
};

template<class QuotientFilter>
using HandleType = typename HandleTypeHelper<QuotientFilter,QuotientFilter::uses_handle>::type;

} // namespace qf
