# Public Core API

This document records the recommended public surface of
`MorphoTreeAdjust` after the repository reorganization.

For general documentation navigation, first see:

- [README.md](../README.md)
- [docs/README.md](./README.md)

## Official Scope

Today, the official API of the main repository is the `core` line, namely:

- per-pixel dynamic trees;
- per-pixel dual min/max tree incremental filtering;
- CASF over that main line;
- Python bindings for the same layer.

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
- `DynamicAreaComputer`
- `DynamicBoundingBoxComputer`

Official example:

- [examples/core_cpp_api_example.cpp](../examples/core_cpp_api_example.cpp)

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

Official example:

- [examples/core_python_api_example.py](../examples/core_python_api_example.py)

Official notebooks:

- [morphoTreeAdjust_example_leaf.ipynb](../notebooks/morphoTreeAdjust_example_leaf.ipynb)
- [morphoTreeAdjust_example_subtree.ipynb](../notebooks/morphoTreeAdjust_example_subtree.ipynb)

## Usage Rules

- prefer `MorphoTreeAdjust.hpp` in new C++ code;
- when documenting the main repository, point first to the `core` API;
- treat `sdist`, `wheel`, and `cmake --install` as distributions of the same
  `core` surface;
- use the examples above whenever you need a small and stable starting point.
