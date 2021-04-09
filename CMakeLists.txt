cmake_minimum_required (VERSION 2.8)

################################################################################
# CMakeLists.txt
#
# cmake build script for our benchmark executables.  Should be fairly obvious
#
# Part of Project lpqfilter - https://github.com/TooBiased/lpqfilter.git
#
# Copyright (C) 2019-2020 Tobias Maier <t.maier@kit.edu>
#
# All rights reserved. Published under the BSD-2 license in the LICENSE file.
################################################################################



#### USER DEFINED ##############################################################

project("quotient_filter")

set(QFILTER_HASHFCT XXH3 CACHE STRING
  "Changes the used hash function if XXHASH is not available, MURMUR2 is used as backoff!")
set_property(CACHE QFILTER_HASHFCT PROPERTY STRINGS XXH3 XXHASH MURMUR2 MURMUR3)

#### BASIC SETTINGS ############################################################

include_directories(${PROJECT_SOURCE_DIR})

if (MSVC)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /std:c++17 /permissive- /MP")
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD /Ox")
  set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd /bigobj")
endif (MSVC)

if (UNIX)
  set (CMAKE_CXX_FLAGS "-std=c++17 -Wall -Wextra -O3 -march=native -g -pthread -static-libstdc++ -Wno-parentheses -Wno-unused-but-set-parameter")

endif (UNIX)

#### FIND PACKAGES #############################################################

find_package(Threads REQUIRED)
set(LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

if (QFILTER_HASHFCT STREQUAL XXHASH)
  message(STATUS "Looking for xxhash")

  find_path(XXHASH_DIR xxHash/xxhash.h
    PATHS ${CMAKE_CURRENT_SOURCE_DIR}/misc/)

  if (XXHASH_DIR)
    include_directories(${XXHASH_DIR}/xxHash)
    message(STATUS "Looking for xxhash -- found")
  else()
    message(STATUS "Looking for xxhash -- not found")
    message("Cannot find xxHash. Use MURMUR2 instead.")
    set(QFILTER_HASHFCT MURMUR2)
  endif()
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/misc)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/misc/cqf/include)

#### MAKE TESTS ################################################################

set(cqf_sources_basic
  benchmarks/basic_test.cpp
  misc/cqf/src/gqf.c
  misc/cqf/src/gqf_file.c
  misc/cqf/src/hashutil.c
  misc/cqf/src/partitioned_counter.c
  )

set(cqf_sources_fill
  benchmarks/fill_ratio_test.cpp
  misc/cqf/src/gqf.c
  misc/cqf/src/gqf_file.c
  misc/cqf/src/hashutil.c
  misc/cqf/src/partitioned_counter.c
)

function( FilterTest variant sources name )
  add_executable(${name} ${sources})
  target_compile_definitions(${name} PRIVATE -D ${variant} -D ${HASH_FCT})
  target_link_libraries(${name} ${LIBRARIES})
endfunction( FilterTest )

FilterTest( BLOOM             benchmarks/basic_test.cpp  basic_bloom_concurrent )
FilterTest( CQF_WRAPPER       "${cqf_sources_basic}"     basic_cqf_wrapper )
FilterTest( QFILTER           benchmarks/basic_test.cpp  basic_compact_qfilter_concurrent )
FilterTest( QFILTER_SEQ       benchmarks/basic_test.cpp  basic_compact_qfilter_seq )
FilterTest( QFILTER_LOCKING   benchmarks/basic_test.cpp  basic_compact_qfilter_locked )
FilterTest( LPFILTER          benchmarks/basic_test.cpp  basic_compact_lpfilter_concurrent )
FilterTest( LPFILTER_SEQ      benchmarks/basic_test.cpp  basic_compact_lpfilter_seq )
# FilterTest( TEMPLATED_QFILTER         benchmarks/basic_test.cpp  basic_templated_qfilter_conc )
# FilterTest( TEMPLATED_QFILTER_SEQ     benchmarks/basic_test.cpp  basic_templated_qfilter_seq )
# FilterTest( TEMPLATED_QFILTER_LOCKING benchmarks/basic_test.cpp  basic_templated_qfilter_locking )
# FilterTest( TEMPLATED_LPFILTER        benchmarks/basic_test.cpp  basic_templated_lpfilter_conc )
# FilterTest( TEMPLATED_LPFILTER_SEQ    benchmarks/basic_test.cpp  basic_templated_lpfilter_seq )


FilterTest( BLOOM                     benchmarks/fill_ratio_test.cpp   fill_bloom_conc  )
FilterTest( CQF_WRAPPER               "${cqf_sources_fill}"            fill_cqf_wrapper )
FilterTest( QFILTER                   benchmarks/fill_ratio_test.cpp   fill_qfilter_conc )
FilterTest( QFILTER_SEQ               benchmarks/fill_ratio_test.cpp   fill_qfilter_seq )
FilterTest( QFILTER_LOCKING           benchmarks/fill_ratio_test.cpp   fill_qfilter_locking )
FilterTest( LPFILTER                  benchmarks/fill_ratio_test.cpp   fill_lpfilter_conc )
FilterTest( LPFILTER_SEQ              benchmarks/fill_ratio_test.cpp   fill_lpfilter_seq )
# FilterTest( TEMPLATED_QFILTER         benchmarks/fill_ratio_test.cpp   fill_templated_qfilter_conc )
# FilterTest( TEMPLATED_QFILTER_SEQ     benchmarks/fill_ratio_test.cpp   fill_templated_qfilter_seq )
# FilterTest( TEMPLATED_QFILTER_LOCKING benchmarks/fill_ratio_test.cpp   fill_templated_qfilter_locking )
# FilterTest( TEMPLATED_LPFILTER        benchmarks/fill_ratio_test.cpp   fill_templated_lpfilter_conc )
# FilterTest( TEMPLATED_LPFILTER_SEQ    benchmarks/fill_ratio_test.cpp   fill_templated_lpfilter_seq )

FilterTest( QFILTER  benchmarks/progression_test.cpp  progression_qfilter_conc_qi )
FilterTest( QFILTER  benchmarks/progression_test.cpp  progression_qfilter_conc_no_qi )
target_compile_definitions(progression_qfilter_conc_no_qi PRIVATE -D NO_QI)
FilterTest( QFILTER  benchmarks/progression_test.cpp  progression_qfilter_conc_no_pre )
target_compile_definitions(progression_qfilter_conc_no_pre PRIVATE -D NO_PRE)
FilterTest( QFILTER  benchmarks/progression_test.cpp  progression_qfilter_conc_no_qi_no_pre )
target_compile_definitions(progression_qfilter_conc_no_qi_no_pre PRIVATE -D NO_QI -D NO_PRE)
