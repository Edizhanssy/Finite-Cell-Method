#pragma once
#include <functional>
#include <vector>

#include "fcm/core/dense_matrix.hpp"
#include "fcm/fcm3d/assembly.hpp"
#include "fcm/fcm3d/config.hpp"
#include "fcm/fcm3d/mesh.hpp"
#include "fcm/fcm3d/quadrature.hpp"

namespace fcm {

struct DirichletBC {
    std::function<bool(const Vec3&)>          pred;
    std::function<double(const Vec3&, int)>   value;
    std::array<bool, 3>                       dirs{true, true, true};
};

struct SolveResult3D {
    Mesh3D                          mesh;
    std::vector<CellQuadrature3D>   quads;
    DenseMatrix                     K;
    std::vector<double>             F;
    std::vector<double>             u;
    int                             n_constrained = 0;
};

SolveResult3D solve(const Config3D& cfg, const std::vector<DirichletBC>& bcs);

std::vector<Vec3> node_coordinates(const Mesh3D& m, const Config3D& cfg);

Vec3 displacement_at(const Config3D& cfg, const Mesh3D& m,
                     const std::vector<double>& u, const Vec3& x);


std::array<double, 6> strain_at(const Config3D& cfg, const Mesh3D& m,
                                const std::vector<double>& u, const Vec3& x);

}  // namespace fcm