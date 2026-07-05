# CASF performance baseline

This file records the CASF benchmark entry points that must remain stable
while general-region update work remains debug/opt-in.

## Branch

Baseline branch for the general-region implementation:

```text
theory/algorithm-alignment
```

## Smoke benchmark

Use the existing CMake script to configure, build, and run the CASF smoke
benchmark:

```bash
cmake -DMTA_ROOT_DIR="$PWD" \
      -P dev-tools/cmake/run_benchmark_component_tree_casf_smoke.cmake
```

The script builds `benchmark_component_tree_casf` and runs:

```text
BM_component_tree_casf_(area|naive_area)/structured/32/32
```

with:

```text
--benchmark_min_time=0.01s
```

## Non-increasing CASF smoke

Use the dedicated CMake script to configure, build, and run the promoted
non-increasing CASF updating path against the naive rebuild oracle:

```bash
cmake -DMTA_ROOT_DIR="$PWD" \
      -P research/cmake/run_benchmark_non_increasing_casf_smoke.cmake
```

The script runs `non_increasing_casf_apply` with the public clean local-update
maintenance mode and the bounding-box area/perimeter attribute on a small
repository image.

## Full benchmark targets

For a broader local run:

```bash
cmake -S dev-tools \
      -B ../build/MorphoTreeAdjust/dev-tools/benchmark_component_tree_casf \
      -DMTA_BUILD_TESTS=OFF \
      -DMTA_BUILD_TOOLS=OFF \
      -DMTA_BUILD_BENCHMARKS=ON

cmake --build ../build/MorphoTreeAdjust/dev-tools/benchmark_component_tree_casf \
      --target benchmark_component_tree_casf -j4

../build/MorphoTreeAdjust/dev-tools/benchmark_component_tree_casf/benchmarks/benchmark_component_tree_casf
```

## CASF update backend comparison

The `benchmark_component_tree_casf` target also contains a paired comparison
for the same CASF threshold workflow using:

- the published subtree backend, `updateTree`;
- the debug explicit-region backend, `updateTreeWithExplicitRegion`.

Build and run a small comparison with:

```bash
cmake -S dev-tools \
      -B ../build/MorphoTreeAdjust/dev-tools/stage8_benchmarks \
      -DMTA_BUILD_TESTS=OFF \
      -DMTA_BUILD_TOOLS=OFF \
      -DMTA_BUILD_BENCHMARKS=ON

cmake --build ../build/MorphoTreeAdjust/dev-tools/stage8_benchmarks \
      --target benchmark_component_tree_casf -j4

../build/MorphoTreeAdjust/dev-tools/stage8_benchmarks/benchmarks/benchmark_component_tree_casf \
      --benchmark_filter='BM_component_tree_casf_(update_tree|explicit_region)_backend_area/structured/32/8$' \
      --benchmark_min_time=0.01s
```

The paired synthetic benchmark is registered for `structured`, `natural_like`,
and `piecewise_scene` image models. To run the full synthetic backend matrix:

```bash
build/dev-tools-benchmarks/benchmarks/benchmark_component_tree_casf \
      --benchmark_filter='BM_component_tree_casf_(update_tree|explicit_region)_backend_area' \
      --benchmark_min_time=0.05s
```

The explicit-region backend reconstructs the same subtree footprint `C` before
each primal prune and then calls `updateTreeWithExplicitRegion`. This measures
the extra cost of externalizing `C` for the CASF workflow; it is still
debug/opt-in and is not the production CASF route.

## Regression rule

The production CASF path must remain on the current subtree-specialized update
until a general-region path is both correct and performance neutral.

Any shared-code refactor must pass:

- the CASF smoke benchmark above;
- the existing unit tests;
- a before/after comparison of the full CASF benchmark target.

Repeated-run regressions above 2% block integration unless the current CASF
path remains available and unchanged.

## Initial smoke run

Local smoke run recorded on 2026-05-14, branch
`theory/algorithm-alignment`.

Command:

```bash
cmake -DMTA_ROOT_DIR="$PWD" \
      -P dev-tools/cmake/run_benchmark_component_tree_casf_smoke.cmake
```

Observed result:

```text
BM_component_tree_casf_area/structured/32/32          178531 ns real, 178175 ns CPU
BM_component_tree_casf_naive_area/structured/32/32   1126663 ns real, 1123833 ns CPU
```

Observed speedup counter:

```text
speedup_vs_naive=6.66764
```

These numbers are a local smoke reference, not a statistically stable
publication benchmark. Use repeated full benchmark runs before accepting any
shared-code refactor into the CASF production path.

## Stage 8 smoke run

Local smoke run recorded on 2026-05-14, branch
`theory/algorithm-alignment`.

CASF smoke:

```text
BM_component_tree_casf_area/structured/32/32           175687 ns real, 175347 ns CPU
BM_component_tree_casf_naive_area/structured/32/32    1132403 ns real, 1131500 ns CPU
speedup_vs_naive=6.76018
```

Decision: keep the generalized region path debug/opt-in. Do not route it
through production CASF before repeated full benchmark runs show a
performance-neutral integration.

## ICDAR 480p backend comparison

Local exploratory run recorded on 2026-05-14, branch
`theory/algorithm-alignment`.

The original ICDAR paths referenced by older local configs were not present on
this machine, so this run used the available 480p grayscale images under:

```text
/Users/wonderalexandre/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_480p
```

Command:

```bash
../build/MorphoTreeAdjust/research/research_compare_update_backends \
  --radio-adj 1.5 \
  /Users/wonderalexandre/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_480p/test_000.png \
  /Users/wonderalexandre/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_480p/test_001.png \
  /Users/wonderalexandre/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_480p/test_002.png \
  /Users/wonderalexandre/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_480p/test_003.png \
  /Users/wonderalexandre/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_480p/test_004.png
```

Configuration:

- image size: 480 x 853;
- attribute: `AREA`;
- adjacency radius: `1.5`;
- thresholds: 50 thresholds from the local 480p CASF configuration;
- output check: every explicit-region run matched the `updateTree` output.

Observed timing:

```text
image        backend                       total_ms  update_ms  phase1_ms  phase2_ms  matches_update_tree
test_000.png updateTree                     75.539    73.843     37.068     36.775     true
test_000.png updateTreeWithExplicitRegion   72.916    55.069     27.855     27.214     true
test_001.png updateTree                     57.551    56.278     29.850     26.428     true
test_001.png updateTreeWithExplicitRegion   64.345    48.426     25.352     23.074     true
test_002.png updateTree                     44.359    43.249     25.526     17.723     true
test_002.png updateTreeWithExplicitRegion   46.240    32.235     19.804     12.431     true
test_003.png updateTree                     49.490    48.289     28.164     20.125     true
test_003.png updateTreeWithExplicitRegion   51.653    37.043     22.278     14.765     true
test_004.png updateTree                     55.001    53.744     33.566     20.178     true
test_004.png updateTreeWithExplicitRegion   59.294    44.056     28.318     15.738     true
```

Average over the five images:

```text
backend                       avg_total_ms  avg_update_ms
updateTree                     56.388        55.081
updateTreeWithExplicitRegion   58.890        43.366
```

Interpretation: on this small single-run sample, the explicit-region backend is
correct for the tested CASF merge cases and reduces the measured update window
by about 21.3%. The end-to-end CASF loop is about 4.4% slower, so this is not
yet evidence to replace the production CASF path. Use repeated runs, randomized
ordering, and a larger ICDAR subset before treating these numbers as stable.

## ICDAR 1080p backend comparison

Local exploratory run recorded on 2026-05-14, branch
`theory/algorithm-alignment`.

Dataset:

```text
/Users/wonderalexandre/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_1080p/test_000.png
...
/Users/wonderalexandre/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_1080p/test_009.png
```

Configuration:

- image size: 1080 x 1920;
- attribute: `AREA`;
- adjacency radius: `1.5`;
- thresholds: `50, 100, ..., 2500`;
- output check: every explicit-region run matched the `updateTree` output.

Observed timing:

```text
image        backend                       total_ms  update_ms  phase1_ms  phase2_ms  matches_update_tree
test_000.png updateTree                    337.079   332.250    183.941    148.309    true
test_000.png updateTreeWithExplicitRegion  354.853   306.094    166.825    139.269    true
test_001.png updateTree                    297.465   293.675    170.357    123.318    true
test_001.png updateTreeWithExplicitRegion  311.779   270.073    159.930    110.143    true
test_002.png updateTree                    236.265   233.358    150.666     82.692    true
test_002.png updateTreeWithExplicitRegion  249.943   209.163    134.081     75.082    true
test_003.png updateTree                    252.820   249.483    152.389     97.094    true
test_003.png updateTreeWithExplicitRegion  263.025   218.991    136.151     82.840    true
test_004.png updateTree                    278.743   275.081    176.069     99.012    true
test_004.png updateTreeWithExplicitRegion  289.594   244.441    161.715     82.726    true
test_005.png updateTree                    287.953   284.596    193.182     91.414    true
test_005.png updateTreeWithExplicitRegion  300.936   256.220    177.767     78.453    true
test_006.png updateTree                    141.831   139.640     78.389     61.251    true
test_006.png updateTreeWithExplicitRegion  146.525   127.404     72.421     54.983    true
test_007.png updateTree                    235.737   232.773    157.719     75.054    true
test_007.png updateTreeWithExplicitRegion  245.264   212.919    148.024     64.895    true
test_008.png updateTree                    169.286   166.950    112.048     54.902    true
test_008.png updateTreeWithExplicitRegion  174.102   141.962     98.035     43.927    true
test_009.png updateTree                    254.780   251.727    163.203     88.524    true
test_009.png updateTreeWithExplicitRegion  274.166   243.590    161.130     82.460    true
```

Average over the ten images:

```text
backend                       avg_total_ms  avg_update_ms  avg_phase1_ms  avg_phase2_ms
updateTree                    249.196       245.953        153.796        92.157
updateTreeWithExplicitRegion  261.019       223.086        141.608        81.478
```

Interpretation: this 1080p run keeps the same correctness signal as 480p. The
explicit-region backend reduces the measured update window by about 9.3%, while
the end-to-end loop is about 4.7% slower. The fixed overhead of collecting and
passing the explicit region is still visible in the total CASF loop, even though
the localized update itself is cheaper.

## ICDAR 1080p backend comparison, 100 images

Local run recorded on 2026-06-04, branch `theory/algorithm-alignment`.

Dataset:

```text
/Volumes/SSD/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_1080p/test_000.png
...
/Volumes/SSD/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_1080p/test_099.png
```

Configuration:

- image size: `1080 x 1920`;
- images: first 100 ICDAR 1080p files;
- attribute: `AREA`;
- adjacency radius: `1.5`;
- thresholds: `50, 100, ..., 2500`;
- output check: every explicit-region run matched the `updateTree` output.

Raw CSV:

```text
research/results/jmiv2026/compare_update_backends_icdar_1080p_100_20260604.csv
```

Command:

```bash
build/dev-tools-tools/tools/compare_update_backends \
  --radio-adj 1.5 \
  --attribute area \
  --thresholds "50,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,1650,1700,1750,1800,1850,1900,1950,2000,2050,2100,2150,2200,2250,2300,2350,2400,2450,2500" \
  /Volumes/SSD/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_1080p/test_000.png \
  ... \
  /Volumes/SSD/GitHub/MorphoTreeAdjust_likeHigra/tests/dat/icdar_1080p/test_099.png
```

Average over the 100 images:

```text
backend                       avg_total_ms  avg_update_ms  avg_phase1_ms  avg_phase2_ms
updateTree                    307.308       303.684        189.842        113.842
updateTreeWithExplicitRegion  316.323       276.073        175.481        100.592
```

Median over the 100 images:

```text
backend                       med_total_ms  med_update_ms  med_phase1_ms  med_phase2_ms
updateTree                    281.495       278.198        171.080        101.721
updateTreeWithExplicitRegion  286.869       245.910        155.541         88.779
```

Pairwise ratios, explicit-region divided by `updateTree`:

```text
metric        avg_ratio  median_ratio  min_ratio  max_ratio
total_ms        1.0304        1.0310    0.9193    1.1343
update_ms       0.9033        0.9028    0.7577    0.9883
phase1_ms       0.9196        0.9185    0.7951    1.0093
phase2_ms       0.8735        0.8760    0.7039    0.9785
```

Counts:

```text
explicit faster in total_ms:   11 / 100 images
explicit faster in update_ms: 100 / 100 images
explicit faster in phase1_ms:  99 / 100 images
explicit faster in phase2_ms: 100 / 100 images
mismatches:                    0 / 100 images
```

Interpretation: this larger run strengthens the earlier signal. The
explicit-region backend consistently reduces the measured update window, but
the whole CASF loop is still about 3.0% slower on average. The likely cost is
outside the dual update itself, mainly the explicit materialization of the
subtree footprint `C` before each call to `updateTreeWithExplicitRegion`.

## Synthetic backend comparison smoke, 2026-06-04

Local smoke run after registering the paired backend benchmark for all three
synthetic image models.

Build:

```bash
cmake -S dev-tools -B build/dev-tools-benchmarks \
      -DMTA_BUILD_TESTS=OFF \
      -DMTA_BUILD_TOOLS=OFF \
      -DMTA_BUILD_BENCHMARKS=ON

cmake --build build/dev-tools-benchmarks \
      --target benchmark_component_tree_casf -j4
```

Commands:

```bash
build/dev-tools-benchmarks/benchmarks/benchmark_component_tree_casf \
      --benchmark_filter='BM_component_tree_casf_(update_tree|explicit_region)_backend_area/(structured|natural_like|piecewise_scene)/128/8$' \
      --benchmark_min_time=0.02s

build/dev-tools-benchmarks/benchmarks/benchmark_component_tree_casf \
      --benchmark_filter='BM_component_tree_casf_(update_tree|explicit_region)_backend_area/(structured|natural_like|piecewise_scene)/256/16$' \
      --benchmark_min_time=0.05s
```

Observed wall-clock benchmark times:

```text
model            size  thresholds  updateTree_ms  explicitRegion_ms
structured        128      8             2.962             2.952
natural_like      128      8             2.454             2.390
piecewise_scene   128      8             2.812             2.406
structured        256     16            12.927            13.105
natural_like      256     16            13.274            12.875
piecewise_scene   256     16            12.104            12.159
```

Each benchmark case validates the backend output against the naive area
baseline before timing. This smoke run indicates that the explicit-region
backend is comparable on these synthetic CASF cases, sometimes slightly faster
and sometimes slightly slower. Treat these values as a regression/smoke signal,
not as a replacement decision for the production CASF path.

## Production CASF regression after Algorithm 2 contract work, 2026-06-04

Local regression after adding the Algorithm 2 flat-raise eligibility API,
contract regression, real-event scanner, and ICDAR contract results. The
production CASF path remains unchanged.

Unit/core regression:

```bash
cmake -S dev-tools \
      -B build/dev-tools-prod-regression \
      -DMTA_BUILD_TESTS=ON \
      -DMTA_BUILD_TOOLS=OFF \
      -DMTA_BUILD_BENCHMARKS=OFF

cmake --build build/dev-tools-prod-regression \
      --target dynamic_component_tree_casf_unit_tests \
               dual_min_max_tree_incremental_filter_unit_tests \
               public_core_api_unit_tests -j4

ctest --test-dir build/dev-tools-prod-regression \
      -R 'dynamic_component_tree_casf_unit_tests|dual_min_max_tree_incremental_filter_unit_tests|public_core_api_unit_tests' \
      --output-on-failure
```

Result:

```text
3/3 tests passed
```

CASF smoke:

```bash
cmake -DMTA_ROOT_DIR="$PWD" \
      -P dev-tools/cmake/run_benchmark_component_tree_casf_smoke.cmake
```

Result:

```text
BM_component_tree_casf_area/structured/32/32           164089 ns real,  163744 ns CPU
BM_component_tree_casf_naive_area/structured/32/32    1073340 ns real, 1071077 ns CPU
```

Interpretation: this verifies the immediate stabilization gate for the current
Algorithm 2 work: the core CASF tests pass and the production CASF smoke still
passes through the existing subtree-specialized path.
