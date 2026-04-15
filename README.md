# MorphoTreeAdjust

[![CI](https://github.com/wonderalexandre/MorphoTreeAdjust/actions/workflows/unit-tests.yml/badge.svg?branch=main)](https://github.com/wonderalexandre/MorphoTreeAdjust/actions/workflows/unit-tests.yml)

MorphoTreeAdjust provides the core C++ and Python implementation for dynamic
updates of component trees and connected alternating sequential filters (CASF).
The repository is organized around a small public API, example programs, and
developer tools for experiments, validation, and benchmarks.

## What You Can Do With It

- build min-trees and max-trees with a shared adjacency model;
- update primal/dual trees after local structural changes instead of rebuilding
  everything from scratch;
- run CASF on grayscale images with dynamic update strategies;
- use the same core concepts from C++ and Python;
- inspect benchmarks and development tools used in the DGMM and JMIV work.

## Quick Start

If you want the Python surface first, install it from PyPI:

```sh
python -m pip install morphoTreeAdjust
```

If you want the current repository version instead, install it from the
repository root:

```sh
python -m pip install .
```

Then use the public module:

```python
import numpy as np
import morphoTreeAdjust as mta

image = np.array(
    [
        [4, 4, 2, 1],
        [4, 3, 2, 1],
        [5, 5, 2, 0],
        [5, 6, 6, 0],
    ],
    dtype=np.uint8,
)

adj = mta.AdjacencyRelation(image.shape[0], image.shape[1], 1.5)
maxtree = mta.DynamicComponentTree(image, True, adj)
mintree = mta.DynamicComponentTree(image, False, adj)

adjust = mta.DynamicComponentTreeAdjustment(mintree, maxtree)
adjust.refreshAttributes()

casf = mta.ComponentTreeCasf(image, "area", adj)
filtered = casf.filter([1, 2])
```

For a complete runnable script, see
[examples/core_python_api_example.py](./examples/core_python_api_example.py).

If you want a developer build with benchmarks and tools:

```sh
cmake -S dev-tools -B ../build/MorphoTreeAdjust/dev-tools
cmake --build ../build/MorphoTreeAdjust/dev-tools
```

If you want the install surface from the root project:

```sh
cmake -S . -B ../build/MorphoTreeAdjust/install \
  -DPYTHON_EXECUTABLE="$(python -c 'import sys; print(sys.executable)')"
cmake --build ../build/MorphoTreeAdjust/install
cmake --install ../build/MorphoTreeAdjust/install --prefix /tmp/mta-install
```

`cmake --install` publishes:

- the Python package under `lib/morphoTreeAdjust`;
- the public C++ headers under `include/MorphoTreeAdjust` when
  `MTA_INSTALL_PUBLIC_CORE_HEADERS=ON`.

## Public API

Recommended public C++ include:

```cpp
#include <MorphoTreeAdjust/MorphoTreeAdjust.hpp>
```

When working directly from a local checkout, the equivalent include is:

```cpp
#include "MorphoTreeAdjust.hpp"
```

Main C++ types:

- `AdjacencyRelation`
- `DynamicComponentTree`
- `DynamicComponentTreeAdjustment`
- `DynamicComponentTreeAdjustmentLeaf`
- `ComponentTreeCasf`
- `DynamicAreaComputer`
- `DynamicBoundingBoxComputer`

Main Python types:

- `AdjacencyRelation`
- `DynamicComponentTree`
- `DynamicComponentTreeAdjustment`
- `ComponentTreeCasf`

Official public examples:

- [examples/core_cpp_api_example.cpp](./examples/core_cpp_api_example.cpp)
- [examples/core_python_api_example.py](./examples/core_python_api_example.py)

Notebooks:

- [notebooks/morphoTreeAdjust_example_leaf.ipynb](./notebooks/morphoTreeAdjust_example_leaf.ipynb)
- [notebooks/morphoTreeAdjust_example_subtree.ipynb](./notebooks/morphoTreeAdjust_example_subtree.ipynb)

For the maintained API contract, see
[docs/public_core_api.md](./docs/public_core_api.md).

## Repository Layout

- Public C++ entrypoint:
  [morphoTreeAdjust/include/MorphoTreeAdjust.hpp](./morphoTreeAdjust/include/MorphoTreeAdjust.hpp)
- Core headers:
  [morphoTreeAdjust/include](./morphoTreeAdjust/include)
- Python bindings:
  [morphoTreeAdjust/morphoTreeAdjust.cpp](./morphoTreeAdjust/morphoTreeAdjust.cpp)
- Examples:
  [examples](./examples)
- Documentation:
  [docs](./docs)
- Developer tools and benchmarks:
  [dev-tools](./dev-tools)

Recommended documentation entry points:

- [docs/README.md](./docs/README.md)
- [docs/public_core_api.md](./docs/public_core_api.md)
- [docs/tree_code_overview.md](./docs/tree_code_overview.md)
- [CONTRIBUTING.md](./CONTRIBUTING.md)

## Contribution and Releases

The repository uses Conventional Commit style pull request titles together with
`release-please` automation.

- merge regular pull requests with **Squash and merge**;
- keep the PR title and the final squash commit title in Conventional Commit
  form, such as `feat: ...`, `fix: ...`, or `docs: ...`;
- when `main` accumulates releasable changes, `release-please` opens or updates
  a release PR;
- the Python package version is derived from Git tags via `setuptools-scm`;
- merging that release PR creates the version tag, updates `CHANGELOG.md`, and
  creates the GitHub release;
- publication to PyPI remains a local manual step.

For the expected contribution workflow, see [CONTRIBUTING.md](./CONTRIBUTING.md).

## Research Context

The repository is motivated by dynamic update problems on component trees.
Given an image `f`, its min-tree and max-tree, and a local change applied on one
tree, the central question is how to update the dual tree efficiently instead of
recomputing it from scratch.

This is the basis for the dynamic primal/dual adjustment work and for CASF
implementations built on top of those structures.

Related papers:

- Wonder Alves, Nicolas Passat, Dênnis José da Silva, Alexandre Morimitsu,
  Ronaldo F. Hashimoto. *Efficient connected alternating sequential filters
  based on component trees*. International Conference on Discrete Geometry and
  Mathematical Morphology (DGMM), November 2025, Groningen, Netherlands.
  [hal-05163556](https://hal.science/hal-05163556/)
- Wonder Alves, Nicolas Passat, Dênnis José da Silva, Alexandre Morimitsu,
  Ronaldo F. Hashimoto. *Component tree: Update rather than rebuild*. Journal
  of Mathematical Imaging and Vision (under review), 2026.

## Developer Tools

Ad hoc debug binaries under `dev-tools/tools/debug_*.cpp` are excluded from the
default developer build. Enable them explicitly when needed:

```sh
cmake -S dev-tools -B ../build/MorphoTreeAdjust/dev-tools \
  -DMTA_BUILD_DEBUG_TOOLS=ON
cmake --build ../build/MorphoTreeAdjust/dev-tools
```

Core benchmark example:

```sh
../build/MorphoTreeAdjust/dev-tools/benchmarks/jmiv2026_benchmark \
  --repeat 5 --warmup 1 --json --no-validate image1.png image2.png
```

### `dynamic_casf_apply`

`dynamic_casf_apply` applies connected alternating sequential filters on
grayscale images and compares the dynamic implementations against the naive
rebuild-based baseline.

After building `dev-tools`, the executable is available at:

```sh
../build/MorphoTreeAdjust/dev-tools/tools/dynamic_casf_apply
```

Usage:

```sh
../build/MorphoTreeAdjust/dev-tools/tools/dynamic_casf_apply \
  [--mode dynamic-subtree|dynamic-leaf|naive|compare] \
  [--attribute area|bbox_width|bbox_height|bbox_diagonal] \
  [--radio-adj <radius>] \
  [--iter-timing] \
  [--no-output] \
  <input.png> [<output.png>] <threshold1> [threshold2 ...]
```

Modes:

- `dynamic-subtree`: CASF with dynamic subtree-based updates;
- `dynamic-leaf`: CASF with dynamic leaf-based updates;
- `naive`: CASF by rebuilding the trees at every threshold;
- `compare`: runs `dynamic-subtree`, `dynamic-leaf`, and `naive`, checks that
  outputs match, and prints a timing summary.

Attributes:

- `area`
- `bbox_width`
- `bbox_height`
- `bbox_diagonal`

Examples:

```sh
../build/MorphoTreeAdjust/dev-tools/tools/dynamic_casf_apply \
  --mode dynamic-subtree \
  --attribute area \
  cameraman.png cameraman_dynamic.png 64 128 256
```

```sh
../build/MorphoTreeAdjust/dev-tools/tools/dynamic_casf_apply \
  --mode compare \
  --attribute bbox_diagonal \
  --radio-adj 1.5 \
  --no-output \
  cameraman.png 64 128 256
```

Timing output:

- `Initialization time`: dynamic setup cost before the threshold loop,
  including adjacency creation, max-tree/min-tree construction, attribute
  computer setup, initial attribute computation, and attribute binding;
- `Iteration time`: accumulated threshold-loop time, excluding initialization
  and final reconstruction;
- `Reconstruction time`: final image reconstruction time;
- `Total time`: `Initialization time + Iteration time + Reconstruction time`.

For `naive`, `Total time` measures the full per-threshold rebuild-and-filter
work. In `compare` mode, the reported dynamic totals include initialization so
they are directly comparable to the naive baseline.
