#pragma once
#include <vector>

#include "fcm/fcm3d/config.hpp"

namespace fcm {

struct CellQuadrature3D {
    std::vector<Vec3>   xi;
    std::vector<double> w;
    std::vector<Vec3>   x;
    std::vector<double> mat;
    // Tensor structure. Node index is: k = leaf*ng^3 + (a*ng + b)*ng + c
    int ng = 0;
    // [leap][axis] -> ng times cell-referenced coordinates
    std::vector<std::array<std::vector<double>, 3>> leaf_xi1d;
};

// Octree leafes tensor-multiplying quadrature.
// reference solution: sum(w) = 8
CellQuadrature3D build_cell_quadrature(const Box& cell, const Config3D& cfg);

inline Vec3 map_to_cell(const Vec3& xi, const Box& cell) {
    Vec3 x{};
    for (int d = 0; d < 3; ++d) {
        const std::size_t i = static_cast<std::size_t>(d);
        x[i] = 0.5 * ((1.0 - xi[i]) * cell.lo[i] + (1.0 + xi[i]) * cell.hi[i]);
    }
    return x;
}

}  // namespace fcm