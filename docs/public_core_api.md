# Public Core API

This document records the recommended public surface of
`MorphoTreeAdjust` after the repository reorganization.

For general documentation navigation, first see:

- [README.md](../README.md)
- [docs/README.md](./README.md)

## Official Scope

The official API of the main repository is the `core` line, namely:

- per-pixel dynamic max-trees and min-trees;
- per-pixel dual min/max tree incremental filtering;
- CASF for order-preserving attributes;
- CASF for the supported non-order-preserving bounding-box attributes;
- Python bindings for the high-level core layer.

Modules outside the standard public surface:

- none

## C++ Entry Point

Recommended header:

```cpp
#include <MorphoTreeAdjust/MorphoTreeAdjust.hpp>
```

When used directly from a local checkout, the equivalent include is:

```cpp
#include "MorphoTreeAdjust.hpp"
```

Main types in the C++ surface:

- `AdjacencyRelation`
- `DynamicComponentTree`
- `DualMinMaxTreeIncrementalFilter`
- `DualMinMaxTreeIncrementalFilterLeaf`
- `ComponentTreeCasf`
- `ComponentTreeCasfNonIncreasing`
- `DynamicAreaComputer`
- `DynamicBoundingBoxComputer`

Official example:

- [examples/core_cpp_api_example.cpp](../examples/core_cpp_api_example.cpp)

## CASF Variants

### `ComponentTreeCasf`

Use `ComponentTreeCasf` for the conventional subtree-pruning CASF path with
order-preserving attributes. Its default attribute is `AREA`, and its
thresholds are integers.

The C++ constructors accept either an adjacency radius or an existing shared
adjacency relation:

```cpp
ComponentTreeCasf(ImageUInt8Ptr image,
                  double radiusAdj,
                  Attribute attribute = AREA);

ComponentTreeCasf(ImageUInt8Ptr image,
                  AdjacencyRelationPtr adjacency,
                  Attribute attribute = AREA);
```

Filtering accepts a sequence of integer thresholds and one of these modes:

- `Mode::Updating`: apply every threshold by incrementally updating the paired
  trees;
- `Mode::Naive`: reconstruct and rebuild the component trees at every
  threshold;
- `Mode::Hybrid`: use the naive path for the first threshold and incremental
  updates for the remaining thresholds.

```cpp
auto output = casf.filter(
    std::vector<int>{10, 20, 30},
    ComponentTreeCasf<>::Mode::Updating);
```

The string overload accepts `"updating"`, `"naive"`, or `"hybrid"`.
`Updating` is the default mode.

### `ComponentTreeCasfNonIncreasing`

Use `ComponentTreeCasfNonIncreasing` for the supported flat
non-order-preserving bounding-box attributes. This path merges selected nodes
locally into their parents and keeps the paired tree and registered attributes
synchronized incrementally.

The accepted C++ attributes and their canonical Python names are:

| C++ `Attribute` | Python string |
| --- | --- |
| `AREA_PER_BOX_PERIMETER` | `"area_per_box_perimeter"` |
| `RECTANGULARITY` | `"rectangularity"` |
| `BOX_ASPECT_RATIO` | `"box_aspect_ratio"` |
| `BOX_ELONGATION` | `"box_elongation"` |
| `AREA_PER_BOX_DIAGONAL_SQUARED` | `"area_per_box_diagonal_squared"` |

`AREA_PER_BOX_PERIMETER` is the C++ default, and
`"area_per_box_perimeter"` is the Python default. Other attributes are
rejected by this class.

The C++ constructors mirror the conventional CASF constructors:

```cpp
ComponentTreeCasfNonIncreasing(
    ImageUInt8Ptr image,
    double radiusAdj,
    Attribute attribute = AREA_PER_BOX_PERIMETER);

ComponentTreeCasfNonIncreasing(
    ImageUInt8Ptr image,
    AdjacencyRelationPtr adjacency,
    Attribute attribute = AREA_PER_BOX_PERIMETER);
```

Its thresholds are floating-point values. Nodes whose current attribute value
is at or below a threshold are merge candidates. Thresholds are applied in
order to the same mutable max-tree/min-tree pair, so a later threshold
continues from the state produced by the earlier ones.

```cpp
auto output = shapeCasf.filter(std::vector<float>{0.50f, 0.75f});
```

The C++ overloads accept either `const std::vector<float>&` or
`std::span<const float>`. This class has no CASF mode argument.

## Small C++ Example

```cpp
#include <MorphoTreeAdjust/MorphoTreeAdjust.hpp>

#include <memory>
#include <vector>

int main() {
    auto image = ImageUInt8::create(2, 3);
    const std::vector<uint8_t> pixels{4, 4, 1, 3, 2, 0};
    for (std::size_t i = 0; i < pixels.size(); ++i) {
        (*image)[static_cast<int>(i)] = pixels[i];
    }

    auto adjacency = std::make_shared<AdjacencyRelation>(2, 3, 1.5);

    ComponentTreeCasf<> areaCasf(image->clone(), adjacency, AREA);
    auto areaOutput = areaCasf.filter(
        std::vector<int>{1, 2},
        ComponentTreeCasf<>::Mode::Updating);

    ComponentTreeCasfNonIncreasing<> shapeCasf(
        image->clone(), adjacency, RECTANGULARITY);
    auto shapeOutput = shapeCasf.filter(std::vector<float>{0.75f});

    return areaOutput && shapeOutput ? 0 : 1;
}
```

## Low-Level C++ Dual-Tree Updates

`DualMinMaxTreeIncrementalFilter` also exposes two node-level operations for
the fixed min-tree/max-tree pair supplied to its constructor:

```cpp
adjust.pruneSubtreeAndUpdateDual(primalTree, dualTree, subtreeRoot);
adjust.mergeNodeIntoParentAndUpdateDual(primalTree, dualTree, nodeId);
```

- `pruneSubtreeAndUpdateDual` removes a complete primal subtree and updates the
  paired dual tree from that subtree's full support;
- `mergeNodeIntoParentAndUpdateDual` merges one primal node into its parent,
  updating the dual tree for every connected component in the node's direct
  proper parts.

Both operations require the exact paired trees used to construct `adjust` and
a live, non-root primal node. If incremental attributes are being tracked,
register their computers and buffers with `setAttributeComputer` before the
operation. These two node-level methods are C++ APIs; they are not exposed by
the current Python binding.

## Python Entry Point

Recommended entry point:

```python
import morphoTreeAdjust as mta
```

Main types in the Python surface:

- `AdjacencyRelation`
- `DynamicComponentTree`
- `DualMinMaxTreeIncrementalFilter`
- `ComponentTreeCasf`
- `ComponentTreeCasfNonIncreasing`

The high-level CASF constructors accept either the default adjacency radius or
an `AdjacencyRelation` object:

```python
area_casf = mta.ComponentTreeCasf(
    image, attribute="area", radiusAdj=1.5
)
area_output = area_casf.filter([10, 20], mode="hybrid")

adj = mta.AdjacencyRelation(image.shape[0], image.shape[1], 1.5)
shape_casf = mta.ComponentTreeCasfNonIncreasing(
    image, "rectangularity", adj
)
shape_output = shape_casf.filter([0.50, 0.75])
```

For `ComponentTreeCasf`, Python thresholds are integer values and `mode`
accepts `"updating"`, `"naive"`, or `"hybrid"`. For
`ComponentTreeCasfNonIncreasing`, Python thresholds are converted to `float`
and there is no `mode` parameter. Both classes expose `getMinTree()`,
`getMaxTree()`, `minTree`, and `maxTree`.

Official example:

- [examples/core_python_api_example.py](../examples/core_python_api_example.py)

Official notebooks:

- [morphoTreeAdjust_example_leaf.ipynb](../notebooks/morphoTreeAdjust_example_leaf.ipynb)
- [morphoTreeAdjust_example_subtree.ipynb](../notebooks/morphoTreeAdjust_example_subtree.ipynb)

## Usage Rules

- prefer `MorphoTreeAdjust.hpp` in new C++ code;
- choose `ComponentTreeCasf` for the conventional order-preserving path and
  `ComponentTreeCasfNonIncreasing` only for its five supported attributes;
- keep integer threshold schedules separate from floating-point shape-attribute
  schedules;
- use `sdist` or `wheel` for a complete Python distribution, including generated
  version metadata;
- use direct `cmake --install` only to stage the native extension and public C++
  headers;
- use the examples above whenever you need a small and stable starting point.
