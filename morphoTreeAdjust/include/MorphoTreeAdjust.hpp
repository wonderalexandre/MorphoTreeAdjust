#pragma once

/**
 * @brief Public entry point for the MorphoTreeAdjust `core`.
 *
 * This header aggregates the C++ surface considered official for the main
 * repository:
 *
 * - per-pixel dynamic trees;
 * - dual min/max tree incremental filtering;
 * - CASF on the main per-pixel line.
 *
 * It defines the published C++ surface of the repository.
 */

#include "AdjacencyRelation.hpp"
#include "AttributeComputer.hpp"
#include "Common.hpp"
#include "ComponentTreeCasf.hpp"
#include "DynamicComponentTree.hpp"
#include "DualMinMaxTreeIncrementalFilter.hpp"
#include "DualMinMaxTreeIncrementalFilterLeaf.hpp"
