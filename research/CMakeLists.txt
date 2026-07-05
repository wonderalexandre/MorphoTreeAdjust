cmake_minimum_required(VERSION 3.15)

project(MorphoTreeAdjustResearch LANGUAGES CXX)

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

get_filename_component(MTA_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
set(MTA_RESEARCH_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(MTA_RESEARCH_STB_DIR "${MTA_ROOT_DIR}/dev-tools/external/stb")

option(MTA_RESEARCH_BUILD_TESTS "Build MorphoTreeAdjust research tests." ON)
option(MTA_RESEARCH_BUILD_TOOLS "Build MorphoTreeAdjust research tools." ON)
option(MTA_INSTALL_RESEARCH_HEADERS "Install MorphoTreeAdjust research headers." OFF)

message(STATUS
  "MorphoTreeAdjust research:"
  " tests=${MTA_RESEARCH_BUILD_TESTS}"
  " tools=${MTA_RESEARCH_BUILD_TOOLS}"
  " install-headers=${MTA_INSTALL_RESEARCH_HEADERS}"
)

if(NOT MSVC)
  set_source_files_properties(
    "${MTA_RESEARCH_STB_DIR}/stb_image.cpp"
    PROPERTIES COMPILE_OPTIONS "-Wno-deprecated-declarations;-Wno-missing-field-initializers"
  )
endif()

set(MTA_RESEARCH_COMMON_SOURCES
  "${MTA_RESEARCH_STB_DIR}/stb_image.cpp"
)

function(add_mta_research_executable target source_file)
  add_executable("${target}"
    "${CMAKE_CURRENT_SOURCE_DIR}/${source_file}"
    ${MTA_RESEARCH_COMMON_SOURCES}
  )
  target_include_directories("${target}" PRIVATE
    "${MTA_ROOT_DIR}"
    "${MTA_RESEARCH_INCLUDE_DIR}"
    "${MTA_RESEARCH_STB_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/tools"
  )
endfunction()

function(add_mta_research_test target source_file)
  add_mta_research_executable("${target}" "tests/${source_file}")
  add_test(NAME "${target}" COMMAND "${target}")
  set_tests_properties("${target}" PROPERTIES LABELS research)
endfunction()

enable_testing()

if(MTA_RESEARCH_BUILD_TESTS)
  add_mta_research_test(research_debug_non_increasing_local_oracle_unit_tests debug_non_increasing_local_oracle_unit_tests.cpp)
  add_mta_research_test(research_debug_problem1_contract_unit_tests debug_problem1_contract_unit_tests.cpp)
  add_mta_research_test(research_debug_problem1_level_materialization_unit_tests debug_problem1_level_materialization_unit_tests.cpp)
  add_mta_research_test(research_debug_problem1_oracle_unit_tests debug_problem1_oracle_unit_tests.cpp)
  add_mta_research_test(research_debug_problem1_survey_unit_tests debug_problem1_survey_unit_tests.cpp)
  add_mta_research_test(research_debug_problem1_updater_unit_tests debug_problem1_updater_unit_tests.cpp)
  add_mta_research_test(research_debug_problem2_subtree_update_contract_unit_tests debug_problem2_subtree_update_contract_unit_tests.cpp)
  add_mta_research_test(research_debug_region_level_change_unit_tests debug_region_level_change_unit_tests.cpp)
  add_mta_research_test(research_debug_region_materialization_unit_tests debug_region_materialization_unit_tests.cpp)
  add_mta_research_test(research_debug_region_merge_updater_unit_tests debug_region_merge_updater_unit_tests.cpp)
  add_mta_research_test(research_debug_region_merge_workspace_unit_tests debug_region_merge_workspace_unit_tests.cpp)
  add_mta_research_test(research_debug_region_rebuild_updater_unit_tests debug_region_rebuild_updater_unit_tests.cpp)
  add_mta_research_test(research_debug_region_split_updater_unit_tests debug_region_split_updater_unit_tests.cpp)
  add_mta_research_test(research_debug_region_update_oracle_unit_tests debug_region_update_oracle_unit_tests.cpp)
  add_mta_research_test(research_debug_structural_comparison_unit_tests debug_structural_comparison_unit_tests.cpp)
  add_mta_research_test(research_details_non_increasing_casf_kernel_unit_tests details_non_increasing_casf_kernel_unit_tests.cpp)
  add_mta_research_test(research_details_non_increasing_local_partition_unit_tests details_non_increasing_local_partition_unit_tests.cpp)
endif()

if(MTA_RESEARCH_BUILD_TOOLS)
  add_mta_research_executable(research_algorithm2_flat_raise_real_events tools/algorithm2_flat_raise_real_events.cpp)
  add_mta_research_executable(research_compare_update_backends tools/compare_update_backends.cpp)
  add_mta_research_executable(research_non_increasing_casf_apply tools/non_increasing_casf_apply.cpp)
endif()

if(MTA_INSTALL_RESEARCH_HEADERS)
  install(DIRECTORY "${MTA_RESEARCH_INCLUDE_DIR}/"
    COMPONENT research_headers
    DESTINATION include
  )
endif()
