#include <array>
#include <cmath>
#include <cstdio>
#include <vector>

#include "fcm/fcm3d/solver.hpp"

int main() {
    fcm::Config3D cfg;
    cfg.lo = {0.0, 0.0, 0.0};
    cfg.hi = {3.0, 1.0, 1.0};
    cfg.n_elements = {2, 1, 1};
    cfg.E     = 1.0;
    cfg.nu    = 0.0;
    cfg.alpha = 1e-8;
    cfg.p     = 4;
    cfg.max_depth = 4;
    cfg.penalty   = 1e5;

    // Fictitious slice: x in [1, 7/3]
    cfg.fictitious = {{{1.0, -1.0, -1.0}, {7.0 / 3.0, 2.0, 2.0}}};

    // f_sin = (1/20) * sin(4 pi x)
    cfg.load_span0 = 0.0;
    cfg.load_span1 = 1.0;
    cfg.load_amp   = 1.0 / 20.0;
    cfg.load_freq  = 4.0 * 3.14159265358979323846;

    // Dirichlet Boundary Conditions
    fcm::DirichletBC left;
    left.pred  = [](const fcm::Vec3& x) { return std::fabs(x[0]) < 1e-12; };
    left.value = [](const fcm::Vec3&, int) { return 0.0; };

    fcm::DirichletBC right;
    right.pred  = [&](const fcm::Vec3& x) { return std::fabs(x[0] - cfg.hi[0]) < 1e-12; };
    right.value = [](const fcm::Vec3&, int d) { return d == 0 ? -1.0 : 0.0; };
    right.dirs  = {true, false, false};

    // solve!!
    const fcm::SolveResult3D r = fcm::solve(cfg, {left, right});

    // spreading quadrature nodes to interpolate!!
    std::size_t npts = 0;
    for (const fcm::CellQuadrature3D& q : r.quads) npts += q.xi.size();
    std::printf("cells=%zu  quad points=%zu  DOF=%d\n",
                r.quads.size(), npts, r.mesh.n_dof_total);

    const double y = 0.5, z = 0.5;
    auto ux = [&](double x) {
        return fcm::displacement_at(cfg, r.mesh, r.u, fcm::Vec3{x, y, z})[0];
    };

    const double u0 = ux(0.0), u1 = ux(1.0), u2 = ux(7.0 / 3.0), u3 = ux(3.0);
    const double mean = (u2 - u1) / (7.0 / 3.0 - 1.0);

    std::printf("the total Degree of Freedom of the whole domain:  %d\n", r.mesh.n_dof_total);
    std::printf("u(0)   = %12.6e   expected  0\n", u0);
    std::printf("u(1)   = %12.6e   expected -3.9789e-03\n", u1);
    std::printf("u(7/3) = %12.6e   expected -1\n", u2);
    std::printf("u(3)   = %12.6e   expected -1\n", u3);
    std::printf("mean strain in fictitious = %8.5f   expected -0.74702\n", mean);
    return 0;
}