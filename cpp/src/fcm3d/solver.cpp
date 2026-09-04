#include "fcm/fcm3d/solver.hpp"

#include <cmath>
#include <stdexcept>

#include "fcm/core/linalg.hpp"
#include "fcm/fcm3d/shape.hpp"

namespace fcm {

std::vector<Vec3> node_coordinates(const Mesh3D& m, const Config3D& cfg) {
    const std::vector<Mode3D> modes = enumerate_modes(cfg.p);
    std::vector<Vec3> xn(static_cast<std::size_t>(m.n_dof),
                         Vec3{std::nan(""), std::nan(""), std::nan("")});

    for (std::size_t e = 0; e < m.cells.size(); ++e) {
        const Box& b = m.cells[e].box;
        for (std::size_t a = 0; a < modes.size(); ++a) {
            if (modes[a].i >= 2 || modes[a].j >= 2 || modes[a].k >= 2) continue;
            const int idx[3] = {modes[a].i, modes[a].j, modes[a].k};
            Vec3 x{};
            for (int d = 0; d < 3; ++d) {
                const std::size_t s = static_cast<std::size_t>(d);
                x[s] = idx[d] == 0 ? b.lo[s] : b.hi[s];
            }
            xn[static_cast<std::size_t>(m.ltog[e][a])] = x;
        }
    }
    return xn;
}

SolveResult3D solve(const Config3D& cfg, const std::vector<DirichletBC>& bcs) {
    SolveResult3D r;
    r.mesh = build_mesh(cfg);
    for (const Cell3D& c : r.mesh.cells)
        r.quads.push_back(build_cell_quadrature(c.box, cfg));

    r.K = assemble_stiffness(r.mesh, cfg, r.quads);
    r.F = assemble_force(r.mesh, cfg, r.quads);

    const std::vector<Vec3> xn = node_coordinates(r.mesh, cfg);
    for (int nd = 0; nd < r.mesh.n_dof; ++nd) {
        const Vec3& x = r.mesh.entity_center[static_cast<std::size_t>(nd)];
        const bool nodal = r.mesh.is_node[static_cast<std::size_t>(nd)] != 0;
        for (const DirichletBC& bc : bcs) {
            if (!bc.pred(x)) continue;
            for (int d = 0; d < 3; ++d) {
                if (!bc.dirs[static_cast<std::size_t>(d)]) continue;
                const int g = 3 * nd + d;

                const double val = nodal ? bc.value(x, d) : 0.0;
                r.K(g, g) += cfg.penalty;
                r.F[static_cast<std::size_t>(g)] += cfg.penalty * val;
                ++r.n_constrained;
            }
        }
    }

    r.u = solve_dense(r.K, r.F);
    return r;
}

Vec3 displacement_at(const Config3D& cfg, const Mesh3D& m,
                     const std::vector<double>& u, const Vec3& x) {
    const std::vector<Mode3D> modes = enumerate_modes(cfg.p);

    for (std::size_t e = 0; e < m.cells.size(); ++e) {
        const Box& b = m.cells[e].box;
        bool in = true;
        for (int d = 0; d < 3; ++d) {
            const std::size_t s = static_cast<std::size_t>(d);
            if (x[s] < b.lo[s] - 1e-12 || x[s] > b.hi[s] + 1e-12) { in = false; break; }
        }
        if (!in) continue;

        Vec3 xi{};
        for (int d = 0; d < 3; ++d) {
            const std::size_t s = static_cast<std::size_t>(d);
            xi[s] = 2.0 * (x[s] - b.lo[s]) / (b.hi[s] - b.lo[s]) - 1.0;
        }
        const ShapeValues3D sv = shape_3d(cfg.p, modes, xi[0], xi[1], xi[2]);

        Vec3 res{0.0, 0.0, 0.0};
        for (std::size_t a = 0; a < modes.size(); ++a)
            for (int d = 0; d < 3; ++d)
                res[static_cast<std::size_t>(d)] +=
                    sv.N[a] * u[static_cast<std::size_t>(3 * m.ltog[e][a] + d)];
        return res;
    }
    throw std::runtime_error("displacement_at: point is not in the domain");
}

    std::array<double, 6> strain_at(const Config3D& cfg, const Mesh3D& m,
                                    const std::vector<double>& u, const Vec3& x) {
    const std::vector<Mode3D> modes = enumerate_modes(cfg.p);

    for (std::size_t e = 0; e < m.cells.size(); ++e) {
        const Box& b = m.cells[e].box;
        bool in = true;
        for (int d = 0; d < 3; ++d) {
            const std::size_t s = static_cast<std::size_t>(d);
            if (x[s] < b.lo[s] - 1e-12 || x[s] > b.hi[s] + 1e-12) { in = false; break; }
        }
        if (!in) continue;

        Vec3 xi{}, jac{};
        for (int d = 0; d < 3; ++d) {
            const std::size_t s = static_cast<std::size_t>(d);
            xi[s]  = 2.0 * (x[s] - b.lo[s]) / (b.hi[s] - b.lo[s]) - 1.0;
            jac[s] = 0.5 * (b.hi[s] - b.lo[s]);
        }

        const ShapeValues3D sv = shape_3d(cfg.p, modes, xi[0], xi[1], xi[2]);
        std::vector<double> B;
        strain_displacement(sv, jac, B);

        const std::size_t ncol = 3 * modes.size();
        std::array<double, 6> eps{};
        for (std::size_t row = 0; row < 6; ++row) {
            double v = 0.0;
            for (std::size_t a = 0; a < modes.size(); ++a)
                for (int d = 0; d < 3; ++d)
                    v += B[row * ncol + 3 * a + static_cast<std::size_t>(d)]
                       * u[static_cast<std::size_t>(3 * m.ltog[e][a] + d)];
            eps[row] = v;
        }
        return eps;
    }
    throw std::runtime_error("strain_at: node is not on an cell!");
}


}  // namespace fcm