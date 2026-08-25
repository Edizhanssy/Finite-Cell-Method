#include "fcm/fcm3d/config.hpp"

#include <cmath>

namespace fcm {

bool Config3D::inside_fictitious(const Vec3& x) const {
    for (const Box& b : fictitious) {
        bool in = true;
        for (int d = 0; d < 3; ++d)
            if (x[static_cast<std::size_t>(d)] < b.lo[static_cast<std::size_t>(d)] ||
                x[static_cast<std::size_t>(d)] > b.hi[static_cast<std::size_t>(d)]) {
                in = false;
                break;
            }
        if (in) return true;
    }
    return false;
}

double Config3D::material_factor(const Vec3& x) const {
    return inside_fictitious(x) ? alpha : 1.0;
}

Vec3 Config3D::body_load(const Vec3& x) const {
    if (load_amp != 0.0 && load_span0 <= x[0] && x[0] <= load_span1)
        return Vec3{load_amp * std::sin(load_freq * x[0]), 0.0, 0.0};
    return Vec3{0.0, 0.0, 0.0};
}

std::array<double, 36> Config3D::elasticity() const {
    std::array<double, 36> C{};
    const double lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu  = E / (2.0 * (1.0 + nu));

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            C[static_cast<std::size_t>(i) * 6 + j] = lam + (i == j ? 2.0 * mu : 0.0);
    for (int i = 3; i < 6; ++i)
        C[static_cast<std::size_t>(i) * 6 + i] = mu;
    return C;
}

}  // namespace fcm