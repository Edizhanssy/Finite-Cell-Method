#include "fcm/fcm3d/assembly.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

int checked = 0, failed = 0;

void check(bool ok, const char* what, double a = 0.0, double b = 0.0) {
    ++checked;
    if (!ok && ++failed <= 10)
        std::printf("FAIL %-34s  %.6e vs %.6e\n", what, a, b);
}

// Ke * v
std::vector<double> matvec(const std::vector<double>& Ke, const std::vector<double>& v) {
    const std::size_t n = v.size();
    std::vector<double> r(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < n; ++j) s += Ke[i * n + j] * v[j];
        r[i] = s;
    }
    return r;
}

double norm_inf(const std::vector<double>& v) {
    double m = 0.0;
    for (double x : v) m = std::fmax(m, std::fabs(x));
    return m;
}

}  // namespace

int main() {
    fcm::Config3D cfg;
    cfg.p  = 3;
    cfg.E  = 2.5;
    cfg.nu = 0.0;

    const fcm::Box cell{{0.0, 0.0, 0.0}, {2.0, 1.5, 0.5}};
    const std::vector<fcm::Mode3D> modes = fcm::enumerate_modes(cfg.p);
    const fcm::CellQuadrature3D q = fcm::build_cell_quadrature(cell, cfg);
    const std::vector<double> Ke = fcm::element_stiffness(cell, cfg, modes, q);

    const std::size_t nm = modes.size();
    const std::size_t n  = 3 * nm;

    // ---- 1) sizes ------------------------------------------------------
    check(q.xi.size() == static_cast<std::size_t>(cfg.n_gauss() * cfg.n_gauss() * cfg.n_gauss()),
          "quadrature nokta sayisi",
          static_cast<double>(q.xi.size()),
          static_cast<double>(cfg.n_gauss() * cfg.n_gauss() * cfg.n_gauss()));
    check(Ke.size() == n * n, "Ke boyutu",
          static_cast<double>(Ke.size()), static_cast<double>(n * n));

    // ---- 2) total weight = volume of the cell ---------------------------------
    {
        double vol = 0.0, detJ = 1.0;
        for (double w : q.w) vol += w;
        for (int d = 0; d < 3; ++d)
            detJ *= 0.5 * (cell.hi[static_cast<std::size_t>(d)] -
                           cell.lo[static_cast<std::size_t>(d)]);
        const double exact = (cell.hi[0]-cell.lo[0]) * (cell.hi[1]-cell.lo[1]) *
                             (cell.hi[2]-cell.lo[2]);
        check(std::fabs(vol * detJ - exact) < 1e-12 * exact, "sum(w)*detJ = hacim",
              vol * detJ, exact);
    }

    double scale = 0.0;
    for (double v : Ke) scale = std::fmax(scale, std::fabs(v));

    // ---- 3) symmetry ------------------
    {
        double worst = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = i + 1; j < n; ++j)
                worst = std::fmax(worst, std::fabs(Ke[i*n+j] - Ke[j*n+i]));
        check(worst < 1e-11 * scale, "Ke simetrik", worst, 0.0);
    }

    // ---- 4) six rigid body modes ------------

    {
        std::vector<fcm::Vec3> node_x;
        for (const fcm::Mode3D& m : modes) {
            if (m.kind != fcm::ModeKind::Node) continue;
            const fcm::Vec3 xi{m.i == 0 ? -1.0 : 1.0,
                               m.j == 0 ? -1.0 : 1.0,
                               m.k == 0 ? -1.0 : 1.0};
            node_x.push_back(fcm::map_to_cell(xi, cell));
        }
        check(node_x.size() == 8, "number of modes:",
              static_cast<double>(node_x.size()), 8.0);

        const char* names[6] = {"displacement x", "displacement y", "displacement z",
                                "torsion x", "torsion y", "torsion z"};
        for (int mode = 0; mode < 6; ++mode) {
            std::vector<double> u(n, 0.0);
            std::size_t nd = 0;
            for (std::size_t m = 0; m < nm; ++m) {
                if (modes[m].kind != fcm::ModeKind::Node) continue;
                const fcm::Vec3& X = node_x[nd++];
                double c[3] = {0.0, 0.0, 0.0};
                switch (mode) {
                    case 0: c[0] = 1.0; break;
                    case 1: c[1] = 1.0; break;
                    case 2: c[2] = 1.0; break;
                    case 3: c[1] = -X[2]; c[2] =  X[1]; break;   // x axis
                    case 4: c[0] =  X[2]; c[2] = -X[0]; break;   // y axis
                    case 5: c[0] = -X[1]; c[1] =  X[0]; break;   // z axis
                }
                for (int d = 0; d < 3; ++d) u[3*m + static_cast<std::size_t>(d)] = c[d];
            }
            const double r = norm_inf(matvec(Ke, u));
            check(r < 1e-10 * scale * norm_inf(u), names[mode], r, 0.0);
        }
    }

    // ---- 5) patch test !! ----------
    // u = (a*x, b*y, c*z). Internal modes are constant
    {
        const double a = 0.013, b = -0.007, c = 0.021;
        std::vector<double> u(n, 0.0);
        std::size_t nd = 0;
        for (std::size_t m = 0; m < nm; ++m) {
            if (modes[m].kind != fcm::ModeKind::Node) continue;
            const fcm::Vec3 xi{modes[m].i == 0 ? -1.0 : 1.0,
                               modes[m].j == 0 ? -1.0 : 1.0,
                               modes[m].k == 0 ? -1.0 : 1.0};
            const fcm::Vec3 X = fcm::map_to_cell(xi, cell);
            u[3*m + 0] = a * X[0];
            u[3*m + 1] = b * X[1];
            u[3*m + 2] = c * X[2];
            ++nd;
        }
        check(nd == 8, "patch: mesh size", static_cast<double>(nd), 8.0);

        fcm::Vec3 jac{};
        for (int d = 0; d < 3; ++d)
            jac[static_cast<std::size_t>(d)] =
                0.5 * (cell.hi[static_cast<std::size_t>(d)] - cell.lo[static_cast<std::size_t>(d)]);

        const double eps_exact[6] = {a, b, c, 0.0, 0.0, 0.0};
        std::vector<double> B;
        double worst_eps = 0.0;
        for (std::size_t k = 0; k < q.xi.size(); ++k) {
            const fcm::ShapeValues3D s =
                fcm::shape_3d(cfg.p, modes, q.xi[k][0], q.xi[k][1], q.xi[k][2]);
            fcm::strain_displacement(s, jac, B);
            for (std::size_t r = 0; r < 6; ++r) {
                double e = 0.0;
                for (std::size_t j = 0; j < n; ++j) e += B[r * n + j] * u[j];
                worst_eps = std::fmax(worst_eps, std::fabs(e - eps_exact[r]));
            }
        }
        check(worst_eps < 1e-12, "patch: sabit gerinim", worst_eps, 0.0);

        const std::vector<double> f = matvec(Ke, u);

        double worst_own = 0.0, surf = 0.0;
        for (std::size_t m = 0; m < nm; ++m) {
            const int idx[3] = {modes[m].i, modes[m].j, modes[m].k};
            for (int d = 0; d < 3; ++d) {
                const double v = std::fabs(f[3*m + static_cast<std::size_t>(d)]);
                if (idx[d] >= 2) worst_own = std::fmax(worst_own, v);
                if (modes[m].kind == fcm::ModeKind::Node) surf = std::fmax(surf, v);
            }
        }
        check(worst_own < 1e-10 * scale, "patch: does not have effect on its own axis, worst_own, 0.0);
        check(surf > 1e-6 * scale,       "patch: effect on the surface",       surf, 0.0);
        //check(worst_int < 1e-10 * scale, "patch: ic mod tepkisi sifir", worst_int, 0.0);
        //check(surf > 1e-6 * scale, "patch: yuzeyde tepki var", surf, 0.0);

        // (c) toplam kuvvet dengesi: sum(f) = 0 her yonde
        for (int d = 0; d < 3; ++d) {
            double s = 0.0;
            for (std::size_t m = 0; m < nm; ++m) s += f[3*m + static_cast<std::size_t>(d)];
            check(std::fabs(s) < 1e-10 * scale, "patch: stability of the force", s, 0.0);
        }
    }

    if (checked < 15) { std::printf("FAILED  assembly3d: only %d checks\n", checked); return 2; }
    std::printf("%s  assembly3d: %d checks, %d failures\n", failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}