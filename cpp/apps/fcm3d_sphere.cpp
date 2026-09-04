#include <cmath>
#include <cstdio>
#include <vector>
#include <array>
#include "fcm/fcm3d/solver.hpp"

int main(int argc, char** argv) {
    const int    p    = argc > 1 ? std::atoi(argv[1]) : 4;
    const int    dep  = argc > 2 ? std::atoi(argv[2]) : 3;
    const int    ne   = argc > 3 ? std::atoi(argv[3]) : 2;

    fcm::Config3D cfg;
    cfg.lo = {0.0, 0.0, 0.0};
    cfg.hi = {1.0, 1.0, 1.0};
    cfg.n_elements = {ne, ne, ne};
    cfg.E     = 1.0;
    cfg.nu    = 0.3;
    cfg.alpha = 1e-8;
    cfg.p         = p;
    cfg.max_depth = dep;
    cfg.penalty   = 1e6;

    const double R = argc > 5 ? std::atof(argv[5]) : 0.15;
    const fcm::Vec3 C{0.5, 0.5, 0.5};
    cfg.fictitious_spheres = {{C, R}};


    const double eps0 = 0.001;

    // Bottom Dirichlet Boundary Condition
    fcm::DirichletBC bottom;
    bottom.pred  = [](const fcm::Vec3& x) { return std::fabs(x[2]) < 1e-12; };
    bottom.value = [](const fcm::Vec3&, int) { return 0.0; };
    bottom.dirs  = {false, false, true};

    // Top Dirichlet Boundary Condition
    fcm::DirichletBC top;
    top.pred  = [](const fcm::Vec3& x) { return std::fabs(x[2] - 1.0) < 1e-12; };
    top.value = [&](const fcm::Vec3&, int d) { return d == 2 ? eps0 : 0.0; };
    top.dirs  = {false, false, true};

    // Rigid body at x=0,y=0; fix u_x, u_y!
    fcm::DirichletBC pin;
    pin.pred  = [](const fcm::Vec3& x) {
        return std::fabs(x[0]) < 1e-12 && std::fabs(x[1]) < 1e-12;
    };
    pin.value = [](const fcm::Vec3&, int) { return 0.0; };
    pin.dirs  = {true, true, false};

    const fcm::SolveResult3D r = fcm::solve(cfg, {bottom, top, pin});

    std::size_t npts = 0;
    for (const fcm::CellQuadrature3D& q : r.quads) npts += q.xi.size();

    // At the center (z = zc) which is close the surface of the sphere sigma_zz, along theta
    const double nu = cfg.nu;
    const double lam = cfg.E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu  = cfg.E / (2.0 * (1.0 + nu));
    const double s_far = cfg.E * eps0;
    const double off = argc > 4 ? std::atof(argv[4]) : 1.02;

    std::printf("p=%d depth=%d cells=%d^3  off=%.2f  R=%.3f  DOF=%d  quad=%zu\n",
            p, dep, ne, off, R, r.mesh.n_dof_total, npts);

    std::printf("# theta_deg  sigma_zz/far\n");
    double smax = 0.0, smin = 1e300, ssum = 0.0;
    const int NT = 72;
    for (int t = 0; t < NT; ++t) {
        const double th = t * 2.0 * M_PI / NT;
        const fcm::Vec3 x{C[0] + off * R * std::cos(th),
                          C[1] + off * R * std::sin(th),
                          C[2]};
        const std::array<double, 6> e = fcm::strain_at(cfg, r.mesh, r.u, x);
        const double tr  = e[0] + e[1] + e[2];
        const double szz = lam * tr + 2.0 * mu * e[2];
        if (t % 6 == 0) std::printf("%6.1f  %10.4f\n", th * 180.0 / M_PI, szz / s_far);
        smax  = std::fmax(smax, szz);
        smin  = std::fmin(smin, szz);
        ssum += szz;
    }
    const double savg = ssum / NT;

    // Radial cut!!
    std::printf("# r/R  sigma_zz/far\n");
    for (int i = 0; i <= 20; ++i) {
        const double rr = 1.0 + i * 0.15;
        const fcm::Vec3 x{C[0] + rr * R, C[1], C[2]};
        if (x[0] > 0.99) break;
        const std::array<double, 6> e = fcm::strain_at(cfg, r.mesh, r.u, x);
        const double tr = e[0] + e[1] + e[2];
        std::printf("%5.2f  %10.4f\n", rr, (lam * tr + 2.0 * mu * e[2]) / s_far);
    }

    const double scf_exact = (27.0 - 15.0 * nu) / (2.0 * (7.0 - 5.0 * nu));
    const double aniso = 100.0 * (smax - smin) / savg;

    std::printf("SCF: avg %.4f (%+.2f%%)   max %.4f (%+.2f%%)   "
                "anisotropy %.1f%%   analytic %.4f\n",
                savg / s_far, 100.0 * (savg / s_far - scf_exact) / scf_exact,
                smax / s_far, 100.0 * (smax / s_far - scf_exact) / scf_exact,
                aniso, scf_exact);

    std::printf("CSV,%d,%d,%d,%.2f,%.3f,%d,%zu,%.6f,%.6f,%.6f,%.6f\n",
            p, dep, ne, off, R, r.mesh.n_dof_total, npts,
            savg / s_far, smax / s_far, smin / s_far, scf_exact);

    return 0;
}