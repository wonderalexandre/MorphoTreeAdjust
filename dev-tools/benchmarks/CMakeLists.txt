function(add_mta_benchmark_target target source_file)
  add_mta_tests_root_executable("${target}" "${source_file}")
endfunction()

function(link_mta_google_benchmark target)
  find_package(benchmark QUIET)
  if(benchmark_FOUND)
    target_link_libraries("${target}" PRIVATE benchmark::benchmark)
    return()
  endif()

  find_path(BENCHMARK_INCLUDE_DIR benchmark/benchmark.h
    PATHS /opt/homebrew/include /usr/local/include
  )
  find_library(BENCHMARK_LIBRARY benchmark
    PATHS /opt/homebrew/lib /usr/local/lib
  )

  if(BENCHMARK_INCLUDE_DIR AND BENCHMARK_LIBRARY)
    target_include_directories("${target}" PRIVATE "${BENCHMARK_INCLUDE_DIR}")
    target_link_libraries("${target}" PRIVATE "${BENCHMARK_LIBRARY}")
    return()
  endif()

  message(FATAL_ERROR "Google Benchmark not found. Install benchmark or provide benchmark headers/libs.")
endfunction()

function(add_mta_google_benchmark_target target)
  set(mta_sources "")
  foreach(source_file IN LISTS ARGN)
    list(APPEND mta_sources "${MTA_TESTS_DIR}/${source_file}")
  endforeach()

  add_executable("${target}"
    ${mta_sources}
    ${MTA_COMMON_SOURCES}
  )
  target_include_directories("${target}" PRIVATE
    "${MTA_ROOT_DIR}"
    "${MTA_INCLUDE_DIR}"
    "${MTA_TESTS_DIR}"
    "${MTA_STB_DIR}"
  )
  link_mta_google_benchmark("${target}")
endfunction()

if(MTA_BUILD_BENCHMARKS)
  add_mta_benchmark_target(jmiv2026_benchmark benchmarks/jmiv2026_benchmark.cpp)
  add_mta_google_benchmark_target(
    benchmark_component_tree_casf
    benchmarks/google_benchmark_main.cpp
    benchmarks/benchmark_component_tree_casf.cpp
  )
endif()
