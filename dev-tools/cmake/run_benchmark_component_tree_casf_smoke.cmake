if(NOT DEFINED MTA_ROOT_DIR)
  message(FATAL_ERROR "MTA_ROOT_DIR is required")
endif()

if(NOT DEFINED MTA_BUILD_ROOT)
  if(DEFINED ENV{MTA_BUILD_ROOT} AND NOT "$ENV{MTA_BUILD_ROOT}" STREQUAL "")
    set(MTA_BUILD_ROOT "$ENV{MTA_BUILD_ROOT}")
  else()
    cmake_path(GET MTA_ROOT_DIR PARENT_PATH MTA_PARENT_DIR)
    set(MTA_BUILD_ROOT "${MTA_PARENT_DIR}/build/MorphoTreeAdjust")
  endif()
endif()

set(BENCHMARK_SOURCE_DIR "${MTA_ROOT_DIR}/dev-tools")
set(BENCHMARK_BUILD_DIR "${MTA_BUILD_ROOT}/dev-tools/benchmark_component_tree_casf_smoke")
set(BENCHMARK_BINARY "${BENCHMARK_BUILD_DIR}/benchmarks/benchmark_component_tree_casf")

execute_process(
  COMMAND "${CMAKE_COMMAND}"
          -S "${BENCHMARK_SOURCE_DIR}"
          -B "${BENCHMARK_BUILD_DIR}"
          -DMTA_BUILD_TESTS=OFF
          -DMTA_BUILD_TOOLS=OFF
          -DMTA_BUILD_BENCHMARKS=ON
  RESULT_VARIABLE configure_result
)
if(NOT configure_result EQUAL 0)
  message(FATAL_ERROR "Failed to configure benchmark_component_tree_casf")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" --build "${BENCHMARK_BUILD_DIR}" --target benchmark_component_tree_casf -j4
  RESULT_VARIABLE build_result
)
if(NOT build_result EQUAL 0)
  message(FATAL_ERROR "Failed to build benchmark_component_tree_casf")
endif()

execute_process(
  COMMAND "${BENCHMARK_BINARY}"
          "--benchmark_filter=BM_component_tree_casf_(area|naive_area)/structured/32/32$"
          "--benchmark_min_time=0.01s"
  RESULT_VARIABLE run_result
)
if(NOT run_result EQUAL 0)
  message(FATAL_ERROR "benchmark_component_tree_casf smoke test failed")
endif()
