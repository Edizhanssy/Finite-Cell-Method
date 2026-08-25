#include "fcm/fcm3d/solver.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {
int checked = 0, failed = 0;
void check(bool ok, const char* what, double a = 0, double b = 0) {
    ++checked;
    if (!ok && ++failed <= 10) std::printf("FAIL %-38s  %.10e vs %.10e\n", what, a, b);
}
}  // namespace

int main() {
    // Patch testi: u = (a*x, b*y, c*z), tum sinir dugumlerinde dayatilir.
    // nu = 0 icin sabit gerinim tam uretilmeli, ic modlar sifir olmali.
    const double A = 0.013, Bc = -0.007, Cc = 0.021;

    fcm::Config3D cfg;
    cfg.p  = 3;
    cfg.nu = 0.0;
    cfg.E  = 2.5;
    cfg.lo = {0.0, 0.0, 0.0};
    cfg.hi = {2.0, 1.5, 0.5};
    cfg.n_elements = {2, 2, 2};
    cfg.penalty = 1e8;

    fcm::DirichletBC bc;
    bc.pred = [&](const fcm::Vec3& x) {
        for (int d = 0; d < 3; ++d) {
            const std::size_t s = static_cast<std::size_t>(d);
            if (std::fabs(x[s] - cfg.lo[s]) < 1e-12 ||
                std::fabs(x[s] - cfg.hi[s]) < 1e-12) return true;
        }
        return false;
    };
    bc.value = [&](const fcm::Vec3& x, int d) {
        const double k[3] = {A, Bc, Cc};
        return k[d] * x[static_cast<std::size_t>(d)];
    };

    const fcm::SolveResult3D r = fcm::solve(cfg, {bc});
    check(r.n_constrained > 0, "kisitli DOF var", r.n_constrained, 0);

    {
        int n_nodal = 0, n_high = 0, n_bnd_nodal = 0, n_bnd_high = 0;
        for (int nd = 0; nd < r.mesh.n_dof; ++nd) {
            const bool nodal = r.mesh.is_node[static_cast<std::size_t>(nd)] != 0;
            nodal ? ++n_nodal : ++n_high;
            if (bc.pred(r.mesh.entity_center[static_cast<std::size_t>(nd)]))
                nodal ? ++n_bnd_nodal : ++n_bnd_high;
        }
        std::printf("n_dof=%d  nodal=%d (bnd %d)  high=%d (bnd %d)  constrained=%d\n",
                    r.mesh.n_dof, n_nodal, n_bnd_nodal, n_high, n_bnd_high, r.n_constrained);
    }

    // 1) Ic DOF'lar (kenar/yuz/hacim) sifir olmali
    double worst_int = 0.0, scale = 0.0;
    for (double v : r.u) scale = std::fmax(scale, std::fabs(v));
    for (int nd = r.mesh.off_edge; nd < r.mesh.n_dof; ++nd)
        for (int d = 0; d < 3; ++d)
            worst_int = std::fmax(worst_int,
                std::fabs(r.u[static_cast<std::size_t>(3 * nd + d)]));
    check(worst_int < 1e-8, "ic mod katsayilari sifir", worst_int, 0.0);

    // 2) Rastgele noktalarda tam lineer alan
    unsigned seed = 987u;
    auto rnd = [&seed]() {
        seed = seed * 1103515245u + 12345u;
        return static_cast<double>((seed >> 16) & 0x7fff) / 32768.0;
    };
    double worst_u = 0.0;
    for (int t = 0; t < 50; ++t) {
        const fcm::Vec3 x{cfg.lo[0] + rnd() * (cfg.hi[0] - cfg.lo[0]),
                          cfg.lo[1] + rnd() * (cfg.hi[1] - cfg.lo[1]),
                          cfg.lo[2] + rnd() * (cfg.hi[2] - cfg.lo[2])};
        const fcm::Vec3 uu = fcm::displacement_at(cfg, r.mesh, r.u, x);
        const double ex[3] = {A * x[0], Bc * x[1], Cc * x[2]};
        for (int d = 0; d < 3; ++d)
            worst_u = std::fmax(worst_u,
                std::fabs(uu[static_cast<std::size_t>(d)] - ex[d]));
    }
    check(worst_int < 1e-8, "ic mod katsayilari sifir", worst_int, 0.0);

    std::printf("worst interior coeff: %.3e   worst u error: %.3e\n", worst_int, worst_u);
    if (checked < 3) { std::printf("FAILED  solve3d: only %d checks\n", checked); return 2; }
    std::printf("%s  solve3d: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}