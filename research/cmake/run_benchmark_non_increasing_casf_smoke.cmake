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

set(SMOKE_SOURCE_DIR "${MTA_ROOT_DIR}/research")
set(SMOKE_BUILD_DIR "${MTA_BUILD_ROOT}/research/benchmark_non_increasing_casf_smoke")
set(SMOKE_BINARY "${SMOKE_BUILD_DIR}/research_non_increasing_casf_apply")
set(SMOKE_IMAGE "${MTA_ROOT_DIR}/dat/misc256/4.1.01.png")

if(NOT EXISTS "${SMOKE_IMAGE}")
  message(FATAL_ERROR "Non-increasing CASF smoke image not found: ${SMOKE_IMAGE}")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}"
          -S "${SMOKE_SOURCE_DIR}"
          -B "${SMOKE_BUILD_DIR}"
          -DMTA_RESEARCH_BUILD_TESTS=OFF
          -DMTA_RESEARCH_BUILD_TOOLS=ON
  RESULT_VARIABLE configure_result
)
if(NOT configure_result EQUAL 0)
  message(FATAL_ERROR "Failed to configure non_increasing_casf_apply")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" --build "${SMOKE_BUILD_DIR}" --target research_non_increasing_casf_apply -j4
  RESULT_VARIABLE build_result
)
if(NOT build_result EQUAL 0)
  message(FATAL_ERROR "Failed to build non_increasing_casf_apply")
endif()

execute_process(
  COMMAND "${SMOKE_BINARY}"
          "--mode" "compare"
          "--attribute" "box-area-perimeter-ratio"
          "--dual-maintenance" "merge-node-clean-unchecked"
          "--no-output"
          "${SMOKE_IMAGE}"
          "1" "2"
  RESULT_VARIABLE run_result
)
if(NOT run_result EQUAL 0)
  message(FATAL_ERROR "non_increasing_casf_apply smoke test failed")
endif()
