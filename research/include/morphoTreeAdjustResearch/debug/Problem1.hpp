#pragma once

#include "StructuralComparison.hpp"
#include "details/ExplicitRegionProblem1CasfGate.hpp"
#include "details/ExplicitRegionProblem1CheckedUpdate.hpp"
#include "details/ExplicitRegionProblem1Contract.hpp"
#include "details/ExplicitRegionProblem1LevelMaterialization.hpp"
#include "details/ExplicitRegionProblem1LocalCorrectnessValidator.hpp"
#include "details/ExplicitRegionProblem1MultiAnchor.hpp"
#include "details/ExplicitRegionProblem1Oracle.hpp"
#include "details/ExplicitRegionProblem1Survey.hpp"
#include "details/ExplicitRegionProblem1Updater.hpp"

namespace debug::problem1 {

using ::debug::region::AnchoredRegionLevelChange;
using ::debug::region::RegionLevelChange;
using ::debug::region::RegionMaterializationPlan;
using ::debug::region::RegionOwnerMaterialization;
using ::debug::region::RegionUpdateDirection;
using ::debug::region::analyzeRegionMaterialization;
using ::debug::region::buildExplicitAnchoredRegionChange;
using ::debug::region::buildExplicitRegionChange;
using ::debug::region::buildSubtreeRegionChange;
using ::debug::region::collectLocalConnectedComponents;
using ::debug::region::directionFromLevels;
using ::debug::region::isConnectedRegion;
using ::debug::structural::compareDynamicComponentTreeStructure;
using ::debug::structural::dynamicComponentTreesAreStructurallyEquivalent;
using ::debug::structural::DynamicComponentTreeStructuralComparisonOptions;
using ::debug::structural::DynamicComponentTreeStructuralComparisonResult;

} // namespace debug::problem1
