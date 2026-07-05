# MorphoTreeAdjust Research

This directory preserves the algorithm-alignment study separately from the
stable public API.

The public product boundary remains:

```text
morphoTreeAdjust/include/
```

The research boundary is:

```text
research/
```

Dependency rule:

```text
morphoTreeAdjust/include/ must not include research/
research/ may include its own header-only study snapshot
```

## Layout

- `include/morphoTreeAdjustResearch/`: header-only snapshot used by the study.
- `tests/`: research and diagnostic regression tests.
- `tools/`: exploratory tools used by the study.
- `docs/`: theory notes, proof sketches, mapping notes, and performance notes.
- `cmake/`: research smoke/regression scripts.
- `scripts/`: auxiliary research scripts.
- `results/`: recorded research outputs.

## Build

Standalone research build:

```bash
cmake -S research -B ../build/MorphoTreeAdjust/research
cmake --build ../build/MorphoTreeAdjust/research -j4
ctest --test-dir ../build/MorphoTreeAdjust/research --output-on-failure
```

Root opt-in build:

```bash
cmake -S . -B ../build/MorphoTreeAdjust/root-research -DMTA_BUILD_RESEARCH=ON
cmake --build ../build/MorphoTreeAdjust/root-research -j4
```

The root build leaves `MTA_BUILD_RESEARCH` off by default.
