#include <cmath>
#include <cstdio>
#include <vector>

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

    const double R = 0.15;
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

    const fcm::SolveResult3D r = fcm::solve(cfg, {bot, top, pin});

    std::size_t npts = 0;
    for (const fcm::CellQuadrature3D& q : r.quads) npts += q.xi.size();
    std::printf("p=%d depth=%d cells=%d^3  DOF=%d  quad points=%zu\n",
                p, dep, ne, r.mesh.n_dof_total, npts);

    // At the center (z = zc) which is close the surface of the sphere sigma_zz, along theta
    const double nu = cfg.nu;
    const double lam = cfg.E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu  = cfg.E / (2.0 * (1.0 + nu));

    const double off = 1.02;
    double smax = 0.0;
    for (int t = 0; t < 72; ++t) {
        const double th = t * 2.0 * M_PI / 72.0;
        const fcm::Vec3 x{C[0] + off * R * std::cos(th),
                          C[1] + off * R * std::sin(th),
                          C[2]};
        const fcm::Vec3 e = fcm::strain_at(cfg, r.mesh, r.u, x);
        const double tr = e[0] + e[1] + e[2];
        const double szz = lam * tr + 2.0 * mu * e[2];
        smax = std::fmax(smax, szz);
    }

    const double s_far = cfg.E * eps0;
    const double scf_exact = (27.0 - 15.0 * nu) / (2.0 * (7.0 - 5.0 * nu));
    std::printf("sigma_zz max = %.6e   far field = %.6e\n", smax, s_far);
    std::printf("SCF = %.4f   analytic = %.4f   error = %.2f%%\n",
                smax / s_far, scf_exact,
                100.0 * (smax / s_far - scf_exact) / scf_exact);
    return 0;
}