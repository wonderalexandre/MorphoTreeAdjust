# MorphoTreeAdjust

MorphoTreeAdjust is the canonical C++/Python repository for the core
algorithms behind the DGMM and JMIV papers on dynamic update of component
trees.

The central problem is:

Let $f$ be an image, and $T^\min_f$ and $T^\max_f$ be the min-tree and the
max-tree of $f$. We set:

- $S$, a subtree in the min-tree $T^\min_f$;
- $T^\min_g$, the min-tree obtained after pruning the subtree $S$ from
  $T^\min_f$;
- $g$, the image reconstructed from the pruned min-tree $T^\min_g$;
- $T^\max_g$, the max-tree of the image $g$.

Problem:

> How to modify $T^\max_f$ to obtain $T^\max_g$?

## Repository focus

This repository is the paper-first codebase for:

- dynamic component trees;
- dynamic primal/dual adjustment;
- CASF filtering over those structures;
- Python bindings for the same core API.

The Python module exposes:

- `AdjacencyRelation`
- `DynamicComponentTree`
- `DynamicComponentTreeAdjustment`
- `ComponentTreeCasf`

## Code map

- Public C++ entrypoint:
  [morphoTreeAdjust/include/MorphoTreeAdjust.hpp](./morphoTreeAdjust/include/MorphoTreeAdjust.hpp)
- Core tree/update headers:
  [morphoTreeAdjust/include](./morphoTreeAdjust/include)
- Python bindings:
  [morphoTreeAdjust/morphoTreeAdjust.cpp](./morphoTreeAdjust/morphoTreeAdjust.cpp)

Documentation:

- [docs/README.md](./docs/README.md)
- [docs/public_core_api.md](./docs/public_core_api.md)
- [docs/tree_code_overview.md](./docs/tree_code_overview.md)

Recommended public include:

```cpp
#include "MorphoTreeAdjust.hpp"
```

Official examples:

- [examples/core_cpp_api_example.cpp](./examples/core_cpp_api_example.cpp)
- [examples/core_python_api_example.py](./examples/core_python_api_example.py)

Notebooks:

- `notebooks/morphoTreeAdjust_example_leaf.ipynb`
- `notebooks/morphoTreeAdjust_example_subtree.ipynb`

## Papers

- Wonder Alves, Nicolas Passat, Dênnis José da Silva, Alexandre Morimitsu,
  Ronaldo F. Hashimoto. *Efficient connected alternating sequential filters
  based on component trees*. International Conference on Discrete Geometry and
  Mathematical Morphology (DGMM), November 2025, Groningen, Netherlands.
  [hal-05163556](https://hal.science/hal-05163556/)
- Wonder Alves, Nicolas Passat, Dênnis José da Silva, Alexandre Morimitsu,
  Ronaldo F. Hashimoto. *Component tree: Update rather than rebuild*. Journal
  of Mathematical Imaging and Vision (under review), 2026.

## Build

Developer build:

```sh
cmake -S dev-tools -B ../build/MorphoTreeAdjust/dev-tools
cmake --build ../build/MorphoTreeAdjust/dev-tools
```

Ad hoc debug binaries under `dev-tools/tools/debug_*.cpp` are excluded from the
default developer build. Enable them explicitly when needed:

```sh
cmake -S dev-tools -B ../build/MorphoTreeAdjust/dev-tools \
  -DMTA_BUILD_DEBUG_TOOLS=ON
cmake --build ../build/MorphoTreeAdjust/dev-tools
```

Install surface:

- `cmake --install` publishes the Python module under `lib/morphoTreeAdjust`;
- when `MTA_INSTALL_PUBLIC_CORE_HEADERS=ON`, it also installs the public C++
  headers into `include/MorphoTreeAdjust`;
- the Python `sdist/wheel` path packages only the files needed for the same
  core surface.

Example:

```sh
cmake -S . -B ../build/MorphoTreeAdjust/install \
  -DPYTHON_EXECUTABLE="$(python -c 'import sys; print(sys.executable)')"
cmake --build ../build/MorphoTreeAdjust/install
cmake --install ../build/MorphoTreeAdjust/install --prefix /tmp/mta-install
```

Core benchmark example:

```sh
../build/MorphoTreeAdjust/dev-tools/benchmarks/jmiv2026_benchmark --repeat 5 --warmup 1 --json --no-validate image1.png image2.png
```

## Tool: `dynamic_casf_apply`

`dynamic_casf_apply` is a dev tool to apply connected alternating sequential
filters (CASF) on grayscale images and compare the dynamic implementations
against the naive rebuild-based baseline.

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
