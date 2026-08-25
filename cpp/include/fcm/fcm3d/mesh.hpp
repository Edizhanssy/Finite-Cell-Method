#pragma once
#include <array>
#include <vector>

#include "fcm/fcm3d/config.hpp"
#include "fcm/fcm3d/shape.hpp"

namespace fcm {

struct Cell3D {
    std::array<int, 3> idx;
    Box box;
};


struct Mesh3D {
    std::array<int, 3> n{};
    int p        = 0;
    int n_modes  = 0;      // hucre basina (p+1)^3
    int n_dof    = 0;      // global skaler
    int n_dof_total = 0;   // 3 * n_dof
        std::vector<Vec3> entity_center;   // n_dof: entitenin geometrik merkezi
        std::vector<char> is_node;         // n_dof: dugum modu mu

    std::vector<Cell3D>           cells;
    std::vector<std::vector<int>> ltog;   // [hucre][mod] -> skaler DOF

    int off_node = 0, off_edge = 0, off_face = 0, off_vol = 0;
};

Mesh3D build_mesh(const Config3D& cfg);

}  // namespace fcm