# Core Technical Guide

This document summarizes the public `MorphoTreeAdjust` core and the two
production CASF paths built on it.

For public use, read first:

- [README.md](../README.md)
- [docs/README.md](./README.md)
- [public_core_api.md](./public_core_api.md)

## Code Map

The central points of the repository are:

- [MorphoTreeAdjust.hpp](../morphoTreeAdjust/include/MorphoTreeAdjust.hpp)
  is the recommended public C++ entry point and includes both CASF variants;
- [AdjacencyRelation.hpp](../morphoTreeAdjust/include/AdjacencyRelation.hpp)
  defines the image domain connectivity shared by the paired trees;
- [DynamicComponentTree.hpp](../morphoTreeAdjust/include/DynamicComponentTree.hpp)
  implements the mutable per-pixel max-tree and min-tree representation;
- [AttributeComputer.hpp](../morphoTreeAdjust/include/AttributeComputer.hpp)
  implements incremental area and bounding-box-derived attributes;
- [DualMinMaxTreeIncrementalFilter.hpp](../morphoTreeAdjust/include/DualMinMaxTreeIncrementalFilter.hpp)
  keeps a fixed max-tree/min-tree pair synchronized after subtree pruning or a
  local node-to-parent merge;
- [DualMinMaxTreeIncrementalFilterLeaf.hpp](../morphoTreeAdjust/include/DualMinMaxTreeIncrementalFilterLeaf.hpp)
  implements the leaf-event variant of the dual-tree update;
- [ComponentTreeCasf.hpp](../morphoTreeAdjust/include/ComponentTreeCasf.hpp)
  implements conventional CASF with integer thresholds and the `Updating`,
  `Naive`, and `Hybrid` modes;
- [ComponentTreeCasfNonIncreasing.hpp](../morphoTreeAdjust/include/ComponentTreeCasfNonIncreasing.hpp)
  implements local-merge CASF for five supported non-order-preserving
  bounding-box attributes and floating-point thresholds;
- [morphoTreeAdjust.cpp](../morphoTreeAdjust/morphoTreeAdjust.cpp)
  exposes the supported high-level layer through pybind11.

## Attribute Paths

The CASF classes have intentionally different contracts.

`ComponentTreeCasf` is the conventional subtree-pruning path. It defaults to
`AREA`, consumes integer thresholds, and supports three execution modes:

- `Updating` applies all thresholds through the incremental dual-tree engine;
- `Naive` reconstructs and rebuilds trees for every threshold;
- `Hybrid` rebuilds for the first threshold and updates incrementally after it.

`ComponentTreeCasfNonIncreasing` is the local node-merge path. It accepts only:

- `AREA_PER_BOX_PERIMETER`;
- `RECTANGULARITY`;
- `BOX_ASPECT_RATIO`;
- `BOX_ELONGATION`;
- `AREA_PER_BOX_DIAGONAL_SQUARED`.

It defaults to `AREA_PER_BOX_PERIMETER`, consumes `float` thresholds, and has no
mode switch. A node is selected when its current attribute value is at or below
the threshold. Each threshold continues from the current mutable tree state.

In Python, the corresponding canonical strings are
`"area_per_box_perimeter"`, `"rectangularity"`,
`"box_aspect_ratio"`, `"box_elongation"`, and
`"area_per_box_diagonal_squared"`. The default is
`"area_per_box_perimeter"`.

## Structural Flow

The shared foundation can be read as follows:

1. `AdjacencyRelation` defines which image pixels are neighbors.
2. `DynamicComponentTree` materializes a mutable max-tree or min-tree.
3. `DynamicAttributeComputer` implementations attach incremental state to tree
   nodes and write values to caller-owned buffers.
4. `DualMinMaxTreeIncrementalFilter` associates one fixed min-tree/max-tree pair
   with the shared adjacency and registered attribute computers.
5. The selected CASF class schedules mutations on both tree directions.
6. Reconstruction of the current min-tree produces the filtered image.
7. `morphoTreeAdjust.cpp` adapts the public high-level objects and image buffers
   to Python.

The conventional `ComponentTreeCasf` path is:

1. compute the current attribute values;
2. collect subtree roots selected by an integer threshold;
3. apply max-tree pruning and synchronize the min-tree;
4. apply min-tree pruning and synchronize the max-tree;
5. either keep updating the pair, rebuild it, or combine both strategies
   according to `Updating`, `Naive`, or `Hybrid`.

The `ComponentTreeCasfNonIncreasing` path is:

1. compute bounding-box attribute values for both trees;
2. collect nodes whose current value is at or below a `float` threshold;
3. sort candidates deterministically by altitude and support signature;
4. merge each still-live candidate into its primal parent with
   `mergeNodeIntoParentAndUpdateDual`;
5. update the dual tree and both registered attribute buffers incrementally;
6. repeat in the opposite tree direction and then for the next threshold.

## Public Node-Level Update Operations

The public C++ adjuster has two complementary mutation units:

```cpp
adjust.pruneSubtreeAndUpdateDual(primalTree, dualTree, subtreeRoot);
adjust.mergeNodeIntoParentAndUpdateDual(primalTree, dualTree, nodeId);
```

`pruneSubtreeAndUpdateDual` builds the explicit region from an entire primal
subtree, updates the paired dual tree, and finally prunes that subtree from the
primal tree. It requires a strict parent/subtree level edge.

`mergeNodeIntoParentAndUpdateDual` builds explicit regions from the connected
components of a node's direct proper parts, updates the dual tree for each
component, and then merges the primal node into its parent. This is the
primitive used by `ComponentTreeCasfNonIncreasing`.

Both calls require:

- the same min-tree/max-tree pair given to the adjuster constructor;
- a live, non-root node in the primal tree;
- a matching shared adjacency relation;
- registered attribute computers and adequately sized buffers when incremental
  attribute values must remain synchronized.

These two methods are public C++ APIs but are not currently bound as Python
methods. The Python CASF classes invoke their respective C++ implementations
internally.

## Small API Examples

Conventional CASF in C++:

```cpp
ComponentTreeCasf<> casf(image, adjacency, AREA);
auto output = casf.filter(
    std::vector<int>{10, 20},
    ComponentTreeCasf<>::Mode::Hybrid);
```

Non-order-preserving CASF in C++:

```cpp
ComponentTreeCasfNonIncreasing<> casf(
    image, adjacency, RECTANGULARITY);
auto output = casf.filter(std::vector<float>{0.50f, 0.75f});
```

The equivalent high-level distinction in Python is:

```python
area_output = mta.ComponentTreeCasf(image, "area", adj).filter(
    [10, 20], mode="updating"
)
shape_output = mta.ComponentTreeCasfNonIncreasing(
    image, "rectangularity", adj
).filter([0.50, 0.75])
```

## Cost Hotspots

The most important costs in the repository are:

- initial construction of the paired `DynamicComponentTree` objects;
- subtree support collection and explicit dual-tree updates in
  `pruneSubtreeAndUpdateDual`;
- connected-component discovery over direct proper parts and repeated explicit
  dual updates in `mergeNodeIntoParentAndUpdateDual`;
- incremental attribute maintenance after structural changes;
- candidate traversal and deterministic sorting at every
  `ComponentTreeCasfNonIncreasing` threshold;
- repeated reconstruction and tree construction in `ComponentTreeCasf`'s
  `Naive` path;
- the number and size of update events in the conventional `Updating` and leaf
  paths.

In practical terms:

- `Naive` trades simpler independent threshold steps for repeated full tree
  construction;
- `Updating` avoids those rebuilds but pays for dual-tree synchronization;
- `Hybrid` pays for one rebuild before switching to incremental updates;
- the leaf variant usually increases the number of processed events;
- subtree pruning concentrates work in larger explicit regions;
- local merging may issue several dual updates for one primal node when its
  direct proper part is disconnected;
- all differences depend on image structure, polarity, attribute, adjacency,
  and threshold schedule.

## Where To Inspect Behavior

To study or validate the core, the most useful places are:

- [unit-tests](../unit-tests)
- [dual_min_max_tree_incremental_filter_unit_tests.cpp](../unit-tests/dual_min_max_tree_incremental_filter_unit_tests.cpp)
- [dynamic_component_tree_casf_unit_tests.cpp](../unit-tests/dynamic_component_tree_casf_unit_tests.cpp)
- [dev-tools/tools](../dev-tools/tools)
- [dev-tools/benchmarks/jmiv2026_benchmark.cpp](../dev-tools/benchmarks/jmiv2026_benchmark.cpp)
- [dev-tools/benchmarks/full_casf_icdar_benchmark.cpp](../dev-tools/benchmarks/full_casf_icdar_benchmark.cpp)
- [examples/core_cpp_api_example.cpp](../examples/core_cpp_api_example.cpp)
- [examples/core_python_api_example.py](../examples/core_python_api_example.py)

## Suggested Reading Order

For the public architecture:

1. [MorphoTreeAdjust.hpp](../morphoTreeAdjust/include/MorphoTreeAdjust.hpp)
2. [AdjacencyRelation.hpp](../morphoTreeAdjust/include/AdjacencyRelation.hpp)
3. [DynamicComponentTree.hpp](../morphoTreeAdjust/include/DynamicComponentTree.hpp)
4. [AttributeComputer.hpp](../morphoTreeAdjust/include/AttributeComputer.hpp)
5. [DualMinMaxTreeIncrementalFilter.hpp](../morphoTreeAdjust/include/DualMinMaxTreeIncrementalFilter.hpp)
6. [ComponentTreeCasf.hpp](../morphoTreeAdjust/include/ComponentTreeCasf.hpp)
7. [ComponentTreeCasfNonIncreasing.hpp](../morphoTreeAdjust/include/ComponentTreeCasfNonIncreasing.hpp)
8. [morphoTreeAdjust.cpp](../morphoTreeAdjust/morphoTreeAdjust.cpp)

For the non-order-preserving path specifically, start at
`ComponentTreeCasfNonIncreasing::filter`, follow
`applyThreshold` and `applyMergePhase`, then read
`DualMinMaxTreeIncrementalFilter::mergeNodeIntoParentAndUpdateDual` and the
explicit-region update helpers it calls. For the conventional subtree path,
start at `ComponentTreeCasf::filter`, compare its three mode branches, and then
read `DualMinMaxTreeIncrementalFilter::pruneSubtreeAndUpdateDual` alongside the
existing batch-pruning wrappers.
