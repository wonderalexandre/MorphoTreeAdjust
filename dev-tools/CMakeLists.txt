cmake_minimum_required(VERSION 3.10)

include("${CMAKE_CURRENT_LIST_DIR}/../cmake/AppleSdkFix.cmake")
project(MorphoTreeAdjustTests LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

if(MSVC)
  add_compile_options(/W4)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")

set(MTA_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
set(MTA_INCLUDE_DIR "${MTA_ROOT_DIR}/morphoTreeAdjust/include")
set(MTA_TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(MTA_STB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external/stb")

include_directories(
  "${MTA_ROOT_DIR}"
  "${MTA_INCLUDE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${MTA_STB_DIR}"
)

if(NOT MSVC)
  set_source_files_properties(
    "${MTA_STB_DIR}/stb_image.cpp"
    PROPERTIES COMPILE_OPTIONS "-Wno-deprecated-declarations;-Wno-missing-field-initializers"
  )
endif()

set(MTA_COMMON_SOURCES
  "${MTA_STB_DIR}/stb_image.cpp"
)

option(MTA_BUILD_TESTS "Build MorphoTreeAdjust unit tests." ON)
option(MTA_BUILD_TOOLS "Build MorphoTreeAdjust developer tools." ON)
option(MTA_BUILD_DEBUG_TOOLS "Build MorphoTreeAdjust ad hoc debug tools." OFF)
option(MTA_BUILD_BENCHMARKS "Build MorphoTreeAdjust benchmarks." ON)

message(STATUS
  "MorphoTreeAdjust targets:"
  " tests=${MTA_BUILD_TESTS}"
  " tools=${MTA_BUILD_TOOLS}"
  " debug-tools=${MTA_BUILD_DEBUG_TOOLS}"
  " benchmarks=${MTA_BUILD_BENCHMARKS}"
)
message(STATUS
  "MorphoTreeAdjust scope: core-only DGMM/JMIV codebase"
)

function(add_mta_test_executable target source_file)
  add_executable("${target}"
    "${CMAKE_CURRENT_SOURCE_DIR}/${source_file}"
    ${MTA_COMMON_SOURCES}
  )
  target_include_directories("${target}" PRIVATE
    "${MTA_ROOT_DIR}"
    "${MTA_INCLUDE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${MTA_STB_DIR}"
  )
endfunction()

function(add_mta_tests_root_executable target source_file)
  add_executable("${target}"
    "${MTA_TESTS_DIR}/${source_file}"
    ${MTA_COMMON_SOURCES}
  )
  target_include_directories("${target}" PRIVATE
    "${MTA_ROOT_DIR}"
    "${MTA_INCLUDE_DIR}"
    "${MTA_TESTS_DIR}"
    "${MTA_STB_DIR}"
  )
endfunction()

enable_testing()

add_subdirectory("${MTA_ROOT_DIR}/unit-tests" "${CMAKE_CURRENT_BINARY_DIR}/unit-tests")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/tools" "${CMAKE_CURRENT_BINARY_DIR}/tools")
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/benchmarks" "${CMAKE_CURRENT_BINARY_DIR}/benchmarks")

add_test(
  NAME benchmark_component_tree_casf_smoke
  COMMAND "${CMAKE_COMMAND}"
          -DMTA_ROOT_DIR=${MTA_ROOT_DIR}
          -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/run_benchmark_component_tree_casf_smoke.cmake"
)
