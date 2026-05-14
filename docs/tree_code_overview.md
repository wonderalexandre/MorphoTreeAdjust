# Core Technical Guide

This document summarizes only the published main line of
`MorphoTreeAdjust`.

For public use, read first:

- [README.md](../README.md)
- [docs/README.md](./README.md)
- [public_core_api.md](./public_core_api.md)

## Code Map

The central points of the repository are:

- [MorphoTreeAdjust.hpp](../morphoTreeAdjust/include/MorphoTreeAdjust.hpp)
  recommended public entry point for C++
- [AdjacencyRelation.hpp](../morphoTreeAdjust/include/AdjacencyRelation.hpp)
  adjacency relation used on the image
- [DynamicComponentTree.hpp](../morphoTreeAdjust/include/DynamicComponentTree.hpp)
  dynamic max-tree and min-tree representation
- [AttributeComputer.hpp](../morphoTreeAdjust/include/AttributeComputer.hpp)
  incremental attributes associated with each node
- [DualMinMaxTreeIncrementalFilter.hpp](../morphoTreeAdjust/include/DualMinMaxTreeIncrementalFilter.hpp)
  subtree-based dual min/max incremental filter
- [DualMinMaxTreeIncrementalFilterLeaf.hpp](../morphoTreeAdjust/include/DualMinMaxTreeIncrementalFilterLeaf.hpp)
  leaf-based dual min/max incremental filter
- [ComponentTreeCasf.hpp](../morphoTreeAdjust/include/ComponentTreeCasf.hpp)
  CASF over the dynamic-tree core
- [morphoTreeAdjust.cpp](../morphoTreeAdjust/morphoTreeAdjust.cpp)
  bindings Python via pybind11

## Structural Reading

The main code flow can be read as follows:

1. `AdjacencyRelation` defines the domain and connectivity.
2. `DynamicComponentTree` materializes the base dynamic tree.
3. `AttributeComputer` and the dynamic attributes track the nodes.
4. `DualMinMaxTreeIncrementalFilter` and
   `DualMinMaxTreeIncrementalFilterLeaf` update the dual tree after pruning in
   the primal tree.
5. `ComponentTreeCasf` organizes that mechanism to apply CASF.
6. `morphoTreeAdjust.cpp` exposes the same layer to Python.

## Cost Hotspots

The most important costs in the repository today are in:

- initial construction of `DynamicComponentTree`;
- incremental update executed by
  `DualMinMaxTreeIncrementalFilter`;
- the leaf-based variant in
  `DualMinMaxTreeIncrementalFilterLeaf`;
- repeated filter application in `ComponentTreeCasf`.

In practical terms:

- the `leaf` variant usually increases the number of processed events;
- the `subtree` variant usually concentrates more work per event;
- the cost difference depends on the image, the polarity, and the CASF
  threshold schedule.

## Where To Inspect Behavior

To study or validate the core, the most useful places are:

- [unit-tests](../unit-tests)
- [dev-tools/tools](../dev-tools/tools)
- [dev-tools/benchmarks/jmiv2026_benchmark.cpp](../dev-tools/benchmarks/jmiv2026_benchmark.cpp)
- [examples/core_cpp_api_example.cpp](../examples/core_cpp_api_example.cpp)
- [examples/core_python_api_example.py](../examples/core_python_api_example.py)

## Suggested Reading Order

If you want to understand the repository core in sequence:

1. [MorphoTreeAdjust.hpp](../morphoTreeAdjust/include/MorphoTreeAdjust.hpp)
2. [DynamicComponentTree.hpp](../morphoTreeAdjust/include/DynamicComponentTree.hpp)
3. [AttributeComputer.hpp](../morphoTreeAdjust/include/AttributeComputer.hpp)
4. [DualMinMaxTreeIncrementalFilter.hpp](../morphoTreeAdjust/include/DualMinMaxTreeIncrementalFilter.hpp)
5. [DualMinMaxTreeIncrementalFilterLeaf.hpp](../morphoTreeAdjust/include/DualMinMaxTreeIncrementalFilterLeaf.hpp)
6. [ComponentTreeCasf.hpp](../morphoTreeAdjust/include/ComponentTreeCasf.hpp)
7. [morphoTreeAdjust.cpp](../morphoTreeAdjust/morphoTreeAdjust.cpp)
