#pragma once
#include <array>
#include <vector>

namespace fcm {

using Vec3 = std::array<double, 3>;

struct Box {
    Vec3 lo;
    Vec3 hi;
};

struct Config3D {
    Vec3 lo{0.0, 0.0, 0.0};
    Vec3 hi{3.0, 1.0, 1.0};
    std::array<int, 3> n_elements{2, 1, 1};

    double E     = 1.0;
    double nu    = 0.0;
    double alpha = 1e-8;

    int    p         = 4;
    int    max_depth = 3;
    double penalty   = 1e5;

    std::vector<Box> fictitious{};

    double load_span0 = 0.0;
    double load_span1 = 0.0;
    double load_amp   = 0.0;
    double load_freq  = 0.0;

    int n_gauss()   const { return p + 1; }
    int n_modes()   const { return (p + 1) * (p + 1) * (p + 1); }
    int n_dof_elem() const { return 3 * n_modes(); }

    bool   inside_fictitious(const Vec3& x) const;
    double material_factor(const Vec3& x) const;
    Vec3   body_load(const Vec3& x) const;

    std::array<double, 36> elasticity() const;
};

}  // namespace fcm