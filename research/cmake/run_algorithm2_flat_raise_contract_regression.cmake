if(NOT DEFINED MTA_ROOT_DIR)
  message(FATAL_ERROR "MTA_ROOT_DIR is required")
endif()

get_filename_component(MTA_ROOT_DIR "${MTA_ROOT_DIR}" ABSOLUTE)

if(NOT DEFINED MTA_BUILD_ROOT)
  if(DEFINED ENV{MTA_BUILD_ROOT} AND NOT "$ENV{MTA_BUILD_ROOT}" STREQUAL "")
    set(MTA_BUILD_ROOT "$ENV{MTA_BUILD_ROOT}")
  else()
    get_filename_component(MTA_PARENT_DIR "${MTA_ROOT_DIR}" DIRECTORY)
    set(MTA_BUILD_ROOT "${MTA_PARENT_DIR}/build/MorphoTreeAdjust")
  endif()
endif()

if(DEFINED MTA_ALGORITHM2_REGRESSION_IMAGE)
  set(REGRESSION_IMAGE "${MTA_ALGORITHM2_REGRESSION_IMAGE}")
elseif(DEFINED ENV{MTA_ALGORITHM2_REGRESSION_IMAGE} AND NOT "$ENV{MTA_ALGORITHM2_REGRESSION_IMAGE}" STREQUAL "")
  set(REGRESSION_IMAGE "$ENV{MTA_ALGORITHM2_REGRESSION_IMAGE}")
else()
  set(REGRESSION_IMAGE "${MTA_ROOT_DIR}/dat/misc256/4.1.01.png")
endif()

if(NOT EXISTS "${REGRESSION_IMAGE}")
  message(FATAL_ERROR "Algorithm 2 regression image not found: ${REGRESSION_IMAGE}")
endif()

set(REGRESSION_SOURCE_DIR "${MTA_ROOT_DIR}/research")
set(REGRESSION_BUILD_DIR "${MTA_BUILD_ROOT}/research/algorithm2_flat_raise_contract_regression")
set(UNIT_TEST_BINARY "${REGRESSION_BUILD_DIR}/research_debug_region_merge_updater_unit_tests")
set(REAL_EVENTS_BINARY "${REGRESSION_BUILD_DIR}/research_algorithm2_flat_raise_real_events")

function(run_checked description)
  execute_process(
    COMMAND ${ARGN}
    RESULT_VARIABLE run_result
    OUTPUT_VARIABLE run_stdout
    ERROR_VARIABLE run_stderr
  )
  if(NOT run_result EQUAL 0)
    message(STATUS "${description} stdout:\n${run_stdout}")
    message(STATUS "${description} stderr:\n${run_stderr}")
    message(FATAL_ERROR "${description} failed with exit code ${run_result}")
  endif()
  if(NOT run_stderr STREQUAL "")
    string(STRIP "${run_stderr}" run_stderr_stripped)
    message(STATUS "${description}: ${run_stderr_stripped}")
  endif()
endfunction()

execute_process(
  COMMAND "${CMAKE_COMMAND}"
          -S "${REGRESSION_SOURCE_DIR}"
          -B "${REGRESSION_BUILD_DIR}"
          -DMTA_RESEARCH_BUILD_TESTS=ON
          -DMTA_RESEARCH_BUILD_TOOLS=ON
  RESULT_VARIABLE configure_result
)
if(NOT configure_result EQUAL 0)
  message(FATAL_ERROR "Failed to configure Algorithm 2 contract regression")
endif()

execute_process(
  COMMAND "${CMAKE_COMMAND}" --build "${REGRESSION_BUILD_DIR}"
          --target research_debug_region_merge_updater_unit_tests
          --target research_algorithm2_flat_raise_real_events
          -j4
  RESULT_VARIABLE build_result
)
if(NOT build_result EQUAL 0)
  message(FATAL_ERROR "Failed to build Algorithm 2 contract regression targets")
endif()

run_checked(
  "Algorithm 2 flat-raise unit contract regression"
  "${UNIT_TEST_BINARY}"
)

run_checked(
  "Algorithm 2 real-image preserve-Q_a regression"
  "${REAL_EVENTS_BINARY}"
  "--max-events-per-image" "5"
  "--max-events-per-node" "2"
  "--min-region-size" "2"
  "--max-region-size" "16"
  "--qa-case" "preserve"
  "--require-events" "5"
  "--fail-on-rejected"
  "--fail-on-mismatch"
  "${REGRESSION_IMAGE}"
)

run_checked(
  "Algorithm 2 real-image release-Q_a regression"
  "${REAL_EVENTS_BINARY}"
  "--max-events-per-image" "5"
  "--max-events-per-node" "2"
  "--min-region-size" "2"
  "--max-region-size" "64"
  "--qa-case" "release"
  "--require-events" "5"
  "--fail-on-rejected"
  "--fail-on-mismatch"
  "${REGRESSION_IMAGE}"
)
