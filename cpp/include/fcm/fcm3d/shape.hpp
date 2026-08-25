#pragma once
#include <vector>

namespace fcm {

enum class ModeKind { Node = 0, Edge = 1, Face = 2, Volume = 3 };

struct Mode3D {
    int i, j, k;
    ModeKind kind;
};

std::vector<Mode3D> enumerate_modes(int p);

struct ShapeValues3D {
    std::vector<double> N, dNdxi, dNdeta, dNdzeta;
};

ShapeValues3D shape_3d(int p, const std::vector<Mode3D>& modes,
                       double xi, double eta, double zeta);

}  // namespace fcm