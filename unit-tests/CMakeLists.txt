set(MTA_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
set(MTA_INCLUDE_DIR "${MTA_ROOT_DIR}/morphoTreeAdjust/include")
set(MTA_SRC_DIR "${MTA_ROOT_DIR}/morphoTreeAdjust/src")
set(MTA_TESTS_DIR "${MTA_ROOT_DIR}/dev-tools")
set(MTA_STB_DIR "${MTA_TESTS_DIR}/external/stb")

include_directories(
  "${MTA_ROOT_DIR}"
  "${MTA_INCLUDE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${MTA_TESTS_DIR}"
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

function(add_mta_unit_test_executable target source_file)
  add_executable("${target}"
    "${CMAKE_CURRENT_SOURCE_DIR}/${source_file}"
    ${MTA_COMMON_SOURCES}
  )
  target_include_directories("${target}" PRIVATE
    "${MTA_ROOT_DIR}"
    "${MTA_INCLUDE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${MTA_TESTS_DIR}"
    "${MTA_STB_DIR}"
  )
endfunction()

function(add_mta_unit_test target label source_file)
  add_mta_unit_test_executable("${target}" "${source_file}")
  add_test(NAME "${target}" COMMAND "${target}")
  set_tests_properties("${target}" PROPERTIES LABELS "${label}")
endfunction()

if(MTA_BUILD_TESTS)
  add_mta_unit_test(adjacency_relation_unit_tests core adjacency_relation_unit_tests.cpp)
  add_mta_unit_test(dynamic_component_tree_adjustment_leaf_unit_tests core dynamic_component_tree_adjustment_leaf_unit_tests.cpp)
  add_mta_unit_test(dynamic_component_tree_adjustment_unit_tests core dynamic_component_tree_adjustment_unit_tests.cpp)
  add_mta_unit_test(dynamic_component_tree_casf_unit_tests core dynamic_component_tree_casf_unit_tests.cpp)
  add_mta_unit_test(dynamic_component_tree_unit_tests core dynamic_component_tree_unit_tests.cpp)
  add_mta_unit_test(public_core_api_unit_tests core public_core_api_unit_tests.cpp)
endif()
