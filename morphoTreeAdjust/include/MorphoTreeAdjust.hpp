#pragma once

/**
 * @brief Public entry point for the MorphoTreeAdjust `core`.
 *
 * This header aggregates the C++ surface considered official for the main
 * repository:
 *
 * - per-pixel dynamic trees;
 * - primal/dual dynamic adjustment;
 * - CASF on the main per-pixel line.
 *
 * It defines the published C++ surface of the repository.
 */

#include "AdjacencyRelation.hpp"
#include "AttributeComputer.hpp"
#include "Common.hpp"
#include "ComponentTreeCasf.hpp"
#include "DynamicComponentTree.hpp"
#include "DynamicComponentTreeAdjustment.hpp"
#include "DynamicComponentTreeAdjustmentLeaf.hpp"
