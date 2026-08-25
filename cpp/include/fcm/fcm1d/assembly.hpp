#pragma once
#include <cstddef>
#include <vector>
#include "fcm/core/dense_matrix.hpp"
#include "fcm/fcm1d/config.hpp"
#include "fcm/fcm1d/mesh.hpp"
#include "fcm/fcm1d/quadrature.hpp"

namespace fcm {

std::vector<double> element_stiffness(const Element1D& el, const Config& cfg,
                                      const ElementQuadrature& q);
std::vector<double> element_force(const Config& cfg, const ElementQuadrature& q);

DenseMatrix         assemble_stiffness(const Mesh& m, const Config& cfg,
                                       const std::vector<ElementQuadrature>& quads);
std::vector<double> assemble_force(const Mesh& m, const Config& cfg,
                                   const std::vector<ElementQuadrature>& quads);

}