#pragma once

#include "details/NonIncreasingLocalOracle.hpp"

namespace debug::non_increasing {

using ::morphotree::details::non_increasing::LocalDualNodeSnapshot;
using ::morphotree::details::non_increasing::LocalMergeNodeContractionDescriptor;
using ::morphotree::details::non_increasing::LocalMergeOracle;
using ::morphotree::details::non_increasing::LocalMergePixelOwner;
using ::morphotree::details::non_increasing::LocalMergePixelRole;
using ::morphotree::details::non_increasing::appendNodeVector;
using ::morphotree::details::non_increasing::appendPixelVector;
using ::morphotree::details::non_increasing::appendPixels;
using ::morphotree::details::non_increasing::applyLocalMergeDescriptorToImage;
using ::morphotree::details::non_increasing::buildLocalDualNodeSnapshot;
using ::morphotree::details::non_increasing::buildLocalMergeNodeContractionDescriptor;
using ::morphotree::details::non_increasing::buildLocalMergeOracle;
using ::morphotree::details::non_increasing::classifyLocalPixel;
using ::morphotree::details::non_increasing::collectBoundaryPixels;
using ::morphotree::details::non_increasing::collectDirectProperParts;
using ::morphotree::details::non_increasing::collectNodeSupportPixels;
using ::morphotree::details::non_increasing::debugString;
using ::morphotree::details::non_increasing::findSortedPixelIndex;
using ::morphotree::details::non_increasing::localMergeOracleDomainSize;
using ::morphotree::details::non_increasing::localMergeOracleHasSparseUpdatedTarget;
using ::morphotree::details::non_increasing::localMergeOracleHasUpdatedTarget;
using ::morphotree::details::non_increasing::localMergeOracleUpdatedPixelValue;
using ::morphotree::details::non_increasing::localMergePixelRoleName;
using ::morphotree::details::non_increasing::makePixelMask;
using ::morphotree::details::non_increasing::sortedConnectedComponents;
using ::morphotree::details::non_increasing::sortedUniquePixels;

} // namespace debug::non_increasing
