#pragma once
#include <vector>

#include "fcm/fcm3d/config.hpp"

namespace fcm {

/// Kesik hucrelerde octree bolme; kesilmeyen bolgeler tek yaprak kalir.
std::vector<Box> partition_cell(const Box& cell, const Config3D& cfg);

}  // namespace fcm