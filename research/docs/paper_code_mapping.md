# Paper notation to code mapping

This note maps the current paper notation for the merge update to the current
implementation and to the new experimental region descriptor.

## Current production scope

The current production CASF path implements the subtree-specialized case:

```cpp
DualMinMaxTreeIncrementalFilter::updateTree(DynamicComponentTree *dualTree,
                                            NodeId subtreeRoot)
```

In that path, `C` is not supplied explicitly. It is collected from the primal
tree as:

```text
C = union of proper parts of all nodes in subtree(subtreeRoot)
```

For the implementation closure of this published merge scope, see
`research/docs/published_merge_closure.md`.

The experimental path introduces an explicit descriptor:

```cpp
debug::region::RegionLevelChange<PixelType>
```

so that the implementation can eventually match the paper formulation directly.
The current experimental execution baseline is:

```cpp
debug::region::RegionLevelChangeRebuildUpdater
```

which applies the explicit region-level change by rebuilding both trees from
the globally modified image. This is a correctness baseline, not the final
local incremental algorithm.

The first local experimental updater is:

```cpp
debug::region::RegionMergeUpdater
```

It currently supports explicit materialized lower updates (`b < a`) on the
dual tree only. Explicit raise updates are intentionally rejected until the
local recomposition for that direction is formalized.

The opt-in method
`RegionMergeUpdater::updateDualTreeOnlyAfterMaterializingSinglePartialOwner`
supports the first non-materialized lower cases: `C` may cut one current
direct owner in the dual tree, optionally together with additional fully
covered owners. The method locally materializes the partial owner and then
runs the same explicit lower-update merge sweep.

The broader opt-in method
`RegionMergeUpdater::updateDualTreeOnlyAfterMaterializingPartialOwners`
materializes every partial owner detected in the current plan before running
the same explicit lower-update merge sweep. The current oracle tests include
multi-owner cuts and a case where one partial owner contributes multiple
inside connected components.

`RegionMergeUpdater::updatePairWithLocalDualAndRebuiltPrimal` closes the
experimental Stage 6 scope by updating the dual tree locally and rebuilding the
primal tree from the updated image, so both tree reconstructions can be checked
against the full-rebuild oracle.

The Stage 7 diagnostic workspace is:

```cpp
debug::region::RegionMergeWorkspace
```

It mirrors the affected-node collection used by the current lower-update sweep
and exposes the intermediate sets using the paper notation. It is read-only and
is not routed through the production CASF path.

## Mapping table

| Paper object | Current production code | Experimental code |
| --- | --- | --- |
| `C` | `properPartSetC_`, collected inside `updateTree(..., subtreeRoot)` | `RegionLevelChange::pixelsC` |
| `a` | `altitudeCa`, derived from `nodeCa` in the dual tree | `RegionLevelChange::a` |
| `b` | `primalTree->getAltitude(parent(subtreeRoot))` | `RegionLevelChange::b` |
| `nodeCa` | local variable `nodeCa` in `updateTree` | `AnchoredRegionLevelChange::nodeCa` |
| `MergeNodes(C, -)` | implicit adjacent seeds around `properPartSetC_` | `RegionMergeWorkspace::mergeNodesCMinus` |
| `MergeNodes([a,b], -)` | `mergeNodesByLevel_` plus climbed paths | `RegionMergeWorkspace::mergeNodesIntervalMinus` |
| `MergeNodes(v, -)` | `mergeNodesByLevel_.getMergedNodes(v)` | `RegionMergeWorkspace::mergeNodesByLevel` |
| `MergedNode(v, +)` | `currentUnionNode` during the level sweep | `RegionMergeWorkspace::mergedNodeVPlusCandidate(...)` |
| frontier nodes | `mergeNodesByLevel_.getFrontierNodesAboveB()` | `RegionMergeWorkspace::frontierNodesAboveB` |
| affected thresholds | `mergeNodesByLevel_.firstMergeLevel()` and `nextMergeLevel()` | `RegionMergeWorkspace::affectedThresholds` |
| affected pixels | `properPartSetC_` and `pixelsInCMarks_` | `RegionLevelChange::pixelsC` plus workspace marks |
| locally affected direct owners | implicit in selected proper-part moves | `RegionMaterializationPlan::affectedOwners` |
| partial owner/proper-part cut | unsupported in production path | detected by `RegionOwnerMaterialization::cutsDirectProperParts()` |
| connected pieces inside/outside `C` | not represented explicitly | `insideComponents` and `outsideComponents` in the materialization plan |
| local node creation | unsupported in production path | `DynamicComponentTree::createDetachedNode` plus `materializeSinglePartialOwner` |
| materialize-then-merge update | unsupported in production path | `RegionMergeUpdater::updateDualTreeOnlyAfterMaterializingSinglePartialOwner` |
| multi-owner materialize-then-merge update | unsupported in production path | `RegionMergeUpdater::updateDualTreeOnlyAfterMaterializingPartialOwners` |
| pair update after materialization | unsupported in production path | `RegionMergeUpdater::updatePairWithLocalDualAndRebuiltPrimal` |

## Current assumptions not present in the paper

The production path currently assumes:

- `C` is generated by a subtree root;
- `C` is a union of existing proper parts;
- `C` does not cut a proper part;
- `b` is the level of `parent(subtreeRoot)`;
- the update is tied to pruning in the primal tree;
- the CASF event is close to a regional-extremum removal.

The paper-general path must eventually remove these assumptions. The first
experimental implementation will relax only the subtree-root assumption while
still requiring `C` to be a union of current proper parts.

## Unsupported cases for the current implementation

The current production path does not directly support:

- an explicit connected region `C` that is not a subtree footprint;
- a region `C` that selects only part of a node's direct proper parts;
- a region update whose target level `b` is not the parent level of a pruned
  subtree;
- multiple `nodeCa` anchors produced by one explicit region;
- local node creation during materialization of `C`.

These cases belong to the experimental implementation stages and must not be
routed through `ComponentTreeCasf` until correctness and performance are
validated.
