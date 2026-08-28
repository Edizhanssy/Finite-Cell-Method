#include "fcm/fcm3d/shape.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {
int checked = 0, failed = 0;

void check(bool ok, const char* what, int a = 0, int b = 0) {
    ++checked;
    if (!ok && ++failed <= 8) std::printf("FAIL %s (%d, %d)\n", what, a, b);
}
}  // namespace

int main() {
    for (int p = 2; p <= 5; ++p) {
        const std::vector<fcm::Mode3D> modes = fcm::enumerate_modes(p);
        const int n = (p + 1) * (p + 1) * (p + 1);
        check(static_cast<int>(modes.size()) == n, "total mode number",
              static_cast<int>(modes.size()), n);

        int cnt[4] = {0, 0, 0, 0};
        for (const fcm::Mode3D& m : modes) ++cnt[static_cast<int>(m.kind)];
        const int q = p - 1;
        check(cnt[0] == 8,           "mesh mode",  cnt[0], 8);
        check(cnt[1] == 12 * q,      "line mode",  cnt[1], 12 * q);
        check(cnt[2] == 6 * q * q,   "face mode",    cnt[2], 6 * q * q);
        check(cnt[3] == q * q * q,   "volume mode",  cnt[3], q * q * q);


        const double pts[4][3] = {{0.3, -0.7, 0.1}, {-0.9, 0.2, 0.55},
                                  {0.0, 0.0, 0.0}, {1.0, -1.0, 0.4}};
        for (const auto& P : pts) {
            const fcm::ShapeValues3D s = fcm::shape_3d(p, modes, P[0], P[1], P[2]);
            double sum = 0.0;
            for (std::size_t m = 0; m < modes.size(); ++m)
                if (modes[m].kind == fcm::ModeKind::Node) sum += s.N[m];
            check(std::fabs(sum - 1.0) < 1e-13, "unit partitioning");

            // 2) Ic modlar ilgili yuzde sifirlanir
            for (std::size_t m = 0; m < modes.size(); ++m) {
                if (modes[m].i < 2) continue;
                const fcm::ShapeValues3D e =
                    fcm::shape_3d(p, modes, 1.0, P[1], P[2]);
                check(std::fabs(e.N[m]) < 1e-13, "internal mode is zero at xi=+1");
                break;
            }
        }

        // 3) Turevler: merkezi fark
        const double h = 1e-5, x = 0.23, y = -0.41, z = 0.62;
        const fcm::ShapeValues3D c  = fcm::shape_3d(p, modes, x, y, z);
        const fcm::ShapeValues3D xp = fcm::shape_3d(p, modes, x + h, y, z);
        const fcm::ShapeValues3D xm = fcm::shape_3d(p, modes, x - h, y, z);
        double worst = 0.0;
        for (std::size_t m = 0; m < modes.size(); ++m) {
            const double fd = (xp.N[m] - xm.N[m]) / (2 * h);
            worst = std::fmax(worst, std::fabs(fd - c.dNdxi[m]));
        }
        check(worst < 1e-6, "dNdxi central difference");
    }

    std::printf("%s  shape3d: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}