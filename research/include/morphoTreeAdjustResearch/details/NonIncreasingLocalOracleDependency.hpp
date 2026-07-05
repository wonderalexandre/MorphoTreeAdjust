#pragma once

#include <sstream>
#include <string>

#include "NonIncreasingLocalPartition.hpp"

namespace morphotree::details::non_increasing {

struct LocalOracleDependencyReport {
    int pixelOwners = 0;
    int pixelOwnersWithOracleDualOwner = 0;
    int preDualSnapshots = 0;
    int oracleDualSnapshots = 0;
    int oracleSnapshotsWithParentLinks = 0;
    int oracleSnapshotsWithChildren = 0;
    int oracleSnapshotsWithDirectProperParts = 0;

    int partitionLevels = 0;
    int partitionComponents = 0;
    int partitionComponentsWithOracleDualOwners = 0;
    int partitionComponentsWithOracleThresholdOwners = 0;
    int partitionComponentsMatchingOracleSnapshot = 0;

    int recompositionNodes = 0;
    int recompositionNodesMatchedToOracle = 0;
    int recompositionNodesWithOracleParent = 0;
    int recompositionNodesWithOracleChildren = 0;

    bool usesOracleDualOwnerLabels() const {
        return pixelOwnersWithOracleDualOwner > 0 ||
               partitionComponentsWithOracleDualOwners > 0 ||
               partitionComponentsWithOracleThresholdOwners > 0;
    }

    bool usesOracleSnapshotSupport() const {
        return oracleDualSnapshots > 0 &&
               partitionComponentsMatchingOracleSnapshot > 0;
    }

    bool usesOraclePartitionGuidance() const {
        return usesOracleDualOwnerLabels() || usesOracleSnapshotSupport();
    }

    bool usesOracleRecompositionGuidance() const {
        return recompositionNodesMatchedToOracle > 0 ||
               recompositionNodesWithOracleParent > 0 ||
               recompositionNodesWithOracleChildren > 0 ||
               oracleSnapshotsWithParentLinks > 0 ||
               oracleSnapshotsWithChildren > 0 ||
               oracleSnapshotsWithDirectProperParts > 0;
    }

    bool isOracleFree() const {
        return !usesOraclePartitionGuidance() &&
               !usesOracleRecompositionGuidance();
    }
};

struct LocalMaterializationOracleDependencyReport {
    int preDualSnapshots = 0;
    int oracleDualSnapshots = 0;
    int pixelOwners = 0;
    int pixelOwnersWithOracleDualOwner = 0;

    int planNodes = 0;
    int planNodesWithOracleNodeId = 0;
    int planNodesWithOracleParent = 0;
    int planNodesWithOracleChildren = 0;
    int planNodesWithStablePreChildren = 0;
    int stablePreChildren = 0;

    int oracleSnapshotsWithParentLinks = 0;
    int oracleSnapshotsWithChildren = 0;
    int oracleSnapshotsWithDirectProperParts = 0;

    bool hasOraclePlanInputs() const {
        return planNodesWithOracleNodeId > 0 ||
               planNodesWithOracleParent > 0 ||
               planNodesWithOracleChildren > 0;
    }

    bool hasOracleSnapshotInputs() const {
        return oracleDualSnapshots > 0 ||
               oracleSnapshotsWithParentLinks > 0 ||
               oracleSnapshotsWithChildren > 0 ||
               oracleSnapshotsWithDirectProperParts > 0;
    }

    bool hasOracleOwnerInputs() const {
        return pixelOwnersWithOracleDualOwner > 0;
    }

    bool usesLocalPreStateInputs() const {
        return preDualSnapshots > 0 ||
               planNodesWithStablePreChildren > 0 ||
               stablePreChildren > 0;
    }

    bool mayUseOracleGuidedMaterialization() const {
        return hasOraclePlanInputs() ||
               hasOracleSnapshotInputs() ||
               hasOracleOwnerInputs();
    }

    bool isOracleFreeMaterializationInput() const {
        return !mayUseOracleGuidedMaterialization();
    }
};

inline LocalOracleDependencyReport auditLocalOracleDependencies(
    const LocalMergeOracle &oracle,
    const LocalPartitionPlan *partitionPlan = nullptr,
    const LocalRecompositionPlan *recompositionPlan = nullptr) {
    LocalOracleDependencyReport report;
    report.pixelOwners = static_cast<int>(oracle.pixelOwners.size());
    for (const auto &owner : oracle.pixelOwners) {
        if (owner.oracleDualOwner != InvalidNode) {
            ++report.pixelOwnersWithOracleDualOwner;
        }
    }

    report.preDualSnapshots = static_cast<int>(oracle.preDualNodes.size());
    report.oracleDualSnapshots = static_cast<int>(oracle.oracleDualNodes.size());
    for (const auto &snapshot : oracle.oracleDualNodes) {
        if (snapshot.parentInSnapshot &&
            snapshot.parentId != InvalidNode &&
            snapshot.parentId != snapshot.nodeId) {
            ++report.oracleSnapshotsWithParentLinks;
        }
        if (!snapshot.childrenInSnapshot.empty()) {
            ++report.oracleSnapshotsWithChildren;
        }
        if (!snapshot.directLocalProperParts.empty()) {
            ++report.oracleSnapshotsWithDirectProperParts;
        }
    }

    if (partitionPlan != nullptr) {
        report.partitionLevels = static_cast<int>(partitionPlan->levels.size());
        for (const auto &levelPartition : partitionPlan->levels) {
            for (const auto &component : levelPartition.components) {
                ++report.partitionComponents;
                if (!component.oracleDualOwners.empty()) {
                    ++report.partitionComponentsWithOracleDualOwners;
                }
                if (!component.oracleThresholdOwners.empty()) {
                    ++report.partitionComponentsWithOracleThresholdOwners;
                }
                if (componentMatchesOracleSnapshot(oracle, component)) {
                    ++report.partitionComponentsMatchingOracleSnapshot;
                }
            }
        }
    }

    if (recompositionPlan != nullptr) {
        report.recompositionNodes = static_cast<int>(recompositionPlan->nodes.size());
        for (const auto &node : recompositionPlan->nodes) {
            if (node.oracleNodeId != InvalidNode) {
                ++report.recompositionNodesMatchedToOracle;
            }
            if (node.oracleParentId != InvalidNode) {
                ++report.recompositionNodesWithOracleParent;
            }
            if (!node.oracleChildrenInSnapshot.empty()) {
                ++report.recompositionNodesWithOracleChildren;
            }
        }
    }

    return report;
}

inline LocalMaterializationOracleDependencyReport auditLocalMaterializationOracleDependencies(
    const LocalMergeOracle &oracle,
    const LocalRecompositionPlan &recompositionPlan) {
    LocalMaterializationOracleDependencyReport report;
    report.preDualSnapshots = static_cast<int>(oracle.preDualNodes.size());
    report.oracleDualSnapshots = static_cast<int>(oracle.oracleDualNodes.size());
    report.pixelOwners = static_cast<int>(oracle.pixelOwners.size());
    for (const auto &owner : oracle.pixelOwners) {
        if (owner.oracleDualOwner != InvalidNode) {
            ++report.pixelOwnersWithOracleDualOwner;
        }
    }

    for (const auto &snapshot : oracle.oracleDualNodes) {
        if (snapshot.parentInSnapshot &&
            snapshot.parentId != InvalidNode &&
            snapshot.parentId != snapshot.nodeId) {
            ++report.oracleSnapshotsWithParentLinks;
        }
        if (!snapshot.childrenInSnapshot.empty()) {
            ++report.oracleSnapshotsWithChildren;
        }
        if (!snapshot.directLocalProperParts.empty()) {
            ++report.oracleSnapshotsWithDirectProperParts;
        }
    }

    report.planNodes = static_cast<int>(recompositionPlan.nodes.size());
    for (const auto &node : recompositionPlan.nodes) {
        if (node.oracleNodeId != InvalidNode) {
            ++report.planNodesWithOracleNodeId;
        }
        if (node.oracleParentId != InvalidNode) {
            ++report.planNodesWithOracleParent;
        }
        if (!node.oracleChildrenInSnapshot.empty()) {
            ++report.planNodesWithOracleChildren;
        }
        if (!node.stableChildPreNodes.empty()) {
            ++report.planNodesWithStablePreChildren;
            report.stablePreChildren += static_cast<int>(node.stableChildPreNodes.size());
        }
    }

    return report;
}

inline LocalOracleDependencyReport auditLocalPlanOracleDependencies(
    const LocalPartitionPlan *partitionPlan = nullptr,
    const LocalRecompositionPlan *recompositionPlan = nullptr) {
    LocalOracleDependencyReport report;

    if (partitionPlan != nullptr) {
        report.partitionLevels = static_cast<int>(partitionPlan->levels.size());
        for (const auto &levelPartition : partitionPlan->levels) {
            for (const auto &component : levelPartition.components) {
                ++report.partitionComponents;
                if (!component.oracleDualOwners.empty()) {
                    ++report.partitionComponentsWithOracleDualOwners;
                }
                if (!component.oracleThresholdOwners.empty()) {
                    ++report.partitionComponentsWithOracleThresholdOwners;
                }
            }
        }
    }

    if (recompositionPlan != nullptr) {
        report.recompositionNodes = static_cast<int>(recompositionPlan->nodes.size());
        for (const auto &node : recompositionPlan->nodes) {
            if (node.oracleNodeId != InvalidNode) {
                ++report.recompositionNodesMatchedToOracle;
            }
            if (node.oracleParentId != InvalidNode) {
                ++report.recompositionNodesWithOracleParent;
            }
            if (!node.oracleChildrenInSnapshot.empty()) {
                ++report.recompositionNodesWithOracleChildren;
            }
        }
    }

    return report;
}

inline std::string debugString(const LocalOracleDependencyReport &report) {
    std::ostringstream out;
    out << std::boolalpha;
    out << "LocalOracleDependencyReport\n";
    out << "pixel_owners=" << report.pixelOwners
        << " pixel_owners_with_oracle_dual_owner=" << report.pixelOwnersWithOracleDualOwner
        << "\n";
    out << "pre_dual_snapshots=" << report.preDualSnapshots
        << " oracle_dual_snapshots=" << report.oracleDualSnapshots
        << " oracle_parent_links=" << report.oracleSnapshotsWithParentLinks
        << " oracle_child_links=" << report.oracleSnapshotsWithChildren
        << " oracle_direct_parts=" << report.oracleSnapshotsWithDirectProperParts
        << "\n";
    out << "partition_levels=" << report.partitionLevels
        << " partition_components=" << report.partitionComponents
        << " partition_components_with_oracle_dual_owners="
        << report.partitionComponentsWithOracleDualOwners
        << " partition_components_with_oracle_threshold_owners="
        << report.partitionComponentsWithOracleThresholdOwners
        << " partition_components_matching_oracle_snapshot="
        << report.partitionComponentsMatchingOracleSnapshot
        << "\n";
    out << "recomposition_nodes=" << report.recompositionNodes
        << " recomposition_nodes_matched_to_oracle="
        << report.recompositionNodesMatchedToOracle
        << " recomposition_nodes_with_oracle_parent="
        << report.recompositionNodesWithOracleParent
        << " recomposition_nodes_with_oracle_children="
        << report.recompositionNodesWithOracleChildren
        << "\n";
    out << "uses_oracle_dual_owner_labels=" << report.usesOracleDualOwnerLabels()
        << " uses_oracle_snapshot_support=" << report.usesOracleSnapshotSupport()
        << " uses_oracle_partition_guidance=" << report.usesOraclePartitionGuidance()
        << " uses_oracle_recomposition_guidance=" << report.usesOracleRecompositionGuidance()
        << " is_oracle_free=" << report.isOracleFree()
        << "\n";
    return out.str();
}

inline std::string debugString(const LocalMaterializationOracleDependencyReport &report) {
    std::ostringstream out;
    out << std::boolalpha;
    out << "LocalMaterializationOracleDependencyReport\n";
    out << "pre_dual_snapshots=" << report.preDualSnapshots
        << " oracle_dual_snapshots=" << report.oracleDualSnapshots
        << " pixel_owners=" << report.pixelOwners
        << " pixel_owners_with_oracle_dual_owner="
        << report.pixelOwnersWithOracleDualOwner
        << "\n";
    out << "plan_nodes=" << report.planNodes
        << " plan_nodes_with_oracle_node_id="
        << report.planNodesWithOracleNodeId
        << " plan_nodes_with_oracle_parent="
        << report.planNodesWithOracleParent
        << " plan_nodes_with_oracle_children="
        << report.planNodesWithOracleChildren
        << " plan_nodes_with_stable_pre_children="
        << report.planNodesWithStablePreChildren
        << " stable_pre_children=" << report.stablePreChildren
        << "\n";
    out << "oracle_parent_links=" << report.oracleSnapshotsWithParentLinks
        << " oracle_child_links=" << report.oracleSnapshotsWithChildren
        << " oracle_direct_parts=" << report.oracleSnapshotsWithDirectProperParts
        << "\n";
    out << "has_oracle_plan_inputs=" << report.hasOraclePlanInputs()
        << " has_oracle_snapshot_inputs=" << report.hasOracleSnapshotInputs()
        << " has_oracle_owner_inputs=" << report.hasOracleOwnerInputs()
        << " uses_local_pre_state_inputs=" << report.usesLocalPreStateInputs()
        << " may_use_oracle_guided_materialization="
        << report.mayUseOracleGuidedMaterialization()
        << " is_oracle_free_materialization_input="
        << report.isOracleFreeMaterializationInput()
        << "\n";
    return out.str();
}

} // namespace morphotree::details::non_increasing
