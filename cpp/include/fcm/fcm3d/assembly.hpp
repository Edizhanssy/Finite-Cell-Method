#pragma once
#include <vector>

#include "fcm/fcm3d/config.hpp"
#include "fcm/fcm3d/quadrature.hpp"
#include "fcm/fcm3d/shape.hpp"

namespace fcm {

void strain_displacement(const ShapeValues3D& s, const Vec3& jac,
                         std::vector<double>& B);

std::vector<double> element_stiffness(const Box& cell, const Config3D& cfg,
                                      const std::vector<Mode3D>& modes,
                                      const CellQuadrature3D& q);

std::vector<double> element_force(const Config3D& cfg,
                                  const std::vector<Mode3D>& modes,
                                  const CellQuadrature3D& q);

}  // namespace fcm