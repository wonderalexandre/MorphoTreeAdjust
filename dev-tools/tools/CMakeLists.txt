function(add_mta_tool_target target source_file)
  add_mta_tests_root_executable("${target}" "${source_file}")
endfunction()

if(MTA_BUILD_TOOLS)
  add_mta_tool_target(compare_subtree_vs_instrumented tools/compare_subtree_vs_instrumented.cpp)
  add_mta_tool_target(dynamic_casf_apply tools/dynamic_casf_apply.cpp)
  add_mta_tool_target(dynamic_component_tree_adjustment_demo tools/dynamic_component_tree_adjustment_demo.cpp)
endif()
