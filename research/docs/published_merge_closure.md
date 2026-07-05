# Published merge closure

This note closes the implementation scope for the merge algorithm used in the
published CASF paper.

## Decision

The published merge path is the current subtree-specialized production path:

```cpp
DualMinMaxTreeIncrementalFilter::updateTree(DynamicComponentTree *dualTree,
                                            NodeId subtreeRoot)
```

and the public CASF wrappers that call it:

```cpp
DualMinMaxTreeIncrementalFilter::pruneMaxTreeAndUpdateMinTree(...)
DualMinMaxTreeIncrementalFilter::pruneMinTreeAndUpdateMaxTree(...)
ComponentTreeCasf::filter(...)
```

This path is considered closed for the published-paper scope. Further work on
arbitrary regions, local splits, local node creation, raise updates, or
multiple-anchor recomposition belongs to the debug/general-region track
and should not be presented as part of the published merge implementation.

## Implemented paper objects

The production implementation realizes the paper merge update in the
subtree-pruning setting:

| Paper object | Production implementation |
| --- | --- |
| `C` | `properPartSetC_`, collected from all proper parts in `subtree(subtreeRoot)` |
| `a` | `altitudeCa`, the altitude of `nodeCa` in the dual tree |
| `b` | altitude of `parent(subtreeRoot)` in the primal tree |
| `nodeCa` | lowest/highest dual anchor selected from the owners of pixels in `C` |
| `MergeNodes(C, -)` | adjacent dual seeds around `C` |
| `MergeNodes([a,b], -)` | nodes collected by climbing the adjacent-seed paths |
| `MergeNodes(v, -)` | `mergeNodesByLevel_.getMergedNodes(v)` |
| `MergedNode(v, +)` | `currentUnionNode` created/reused during the level sweep |
| frontier nodes | `mergeNodesByLevel_.getFrontierNodesAboveB()` |

The diagnostic debug workspace
`debug::region::RegionMergeWorkspace` mirrors these collections with paper names
for inspection and tests, but it is not part of the production CASF path.

## Production invariants

The published merge implementation is accepted only under the existing
subtree-specialized assumptions:

- `C` is the footprint of a live non-root subtree in the primal tree;
- `C` is a union of current proper parts;
- the primal tree is pruned after the dual tree has been updated;
- the target level `b` is obtained from `parent(subtreeRoot)`;
- the dual tree reconstruction must match the pruned primal tree;
- dynamic attributes must remain consistent with full recomputation when
  validation is enabled.

These assumptions are implementation constraints, not general mathematical
limitations of the paper formulation.

## Correctness evidence

The closure relies on these test groups:

- `dual_min_max_tree_incremental_filter_unit_tests`
  - checks production `updateTree` logging for `C`, `F_lambda`, and `nodeCa`;
  - compares dynamic CASF flow against a naive prune/rebuild baseline;
  - covers a recorded subtree regression on `dat/cameraman.pgm`;
  - validates sequential min-tree prunes and reconstruction equality;
  - checks tree connectivity and proper-part consistency after updates;
  - exercises all shared non-root subtree ids on the 12x12 fixture.
- `dynamic_component_tree_casf_unit_tests`
  - compares public `ComponentTreeCasf` area and bounding-box sequences against
    the naive rebuild baseline;
  - includes structured stress images up to `256 x 256`.
- `debug_region_level_change_unit_tests`
  - verifies that the subtree descriptor builder recovers the same `C`, `a`,
    `b`, and `nodeCa` used by the production update.
- `debug_region_merge_updater_unit_tests`
  - checks that a subtree-derived explicit descriptor matches the production
    dual update and the rebuild oracle.
- `debug_region_merge_workspace_unit_tests`
  - verifies that the Stage 7 workspace exposes the merge collections using the
    paper notation.

## Performance evidence

The production CASF benchmark target remains:

```text
benchmark_component_tree_casf
```

The smoke benchmark remains:

```bash
cmake -DMTA_ROOT_DIR="$PWD" \
      -P dev-tools/cmake/run_benchmark_component_tree_casf_smoke.cmake
```

Local smoke run recorded on 2026-05-14, branch
`theory/algorithm-alignment`:

```text
BM_component_tree_casf_area/structured/32/32           175687 ns real, 175347 ns CPU
BM_component_tree_casf_naive_area/structured/32/32    1132403 ns real, 1131500 ns CPU
speedup_vs_naive=6.76018
```

## Explicitly out of scope

The following items are intentionally excluded from the published merge
closure:

- arbitrary connected `C` that is not a subtree footprint;
- `C` cutting a current direct proper part;
- local split/update with node creation;
- multiple `nodeCa` anchors;
- explicit raise updates;
- replacing the production CASF hot path with the debug workspace;
- claiming performance neutrality for generalized-region code.

These items are valid research directions, but they belong to a new paper or a
future general-region implementation plan.

## Maintenance rule

Bug fixes and regression tests may continue on the production merge path. New
general-region behavior should remain under `research/include/morphoTreeAdjustResearch/debug`
or equivalent opt-in code until it has independent correctness and performance
evidence.
