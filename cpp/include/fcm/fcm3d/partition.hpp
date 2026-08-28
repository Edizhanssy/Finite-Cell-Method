#pragma once
#include <vector>

#include "fcm/fcm3d/config.hpp"

namespace fcm {

// Octree partition on cut-cells
std::vector<Box> partition_cell(const Box& cell, const Config3D& cfg);

}  // namespace fcm