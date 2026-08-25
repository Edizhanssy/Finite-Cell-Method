#include "fcm/fcm3d/assembly.hpp"

#include <cstddef>

namespace fcm {

void strain_displacement(const ShapeValues3D& s, const Vec3& jac,
                         std::vector<double>& B) {
    const std::size_t nm  = s.N.size();
    const std::size_t ncol = 3 * nm;
    B.assign(6 * ncol, 0.0);

    for (std::size_t m = 0; m < nm; ++m) {
        const double dx = s.dNdxi[m]   / jac[0];
        const double dy = s.dNdeta[m]  / jac[1];
        const double dz = s.dNdzeta[m] / jac[2];
        const std::size_t c = 3 * m;

        B[0 * ncol + c + 0] = dx;                 // eps_xx
        B[1 * ncol + c + 1] = dy;                 // eps_yy
        B[2 * ncol + c + 2] = dz;                 // eps_zz
        B[3 * ncol + c + 1] = dz;                 // gamma_yz
        B[3 * ncol + c + 2] = dy;
        B[4 * ncol + c + 0] = dz;                 // gamma_xz
        B[4 * ncol + c + 2] = dx;
        B[5 * ncol + c + 0] = dy;                 // gamma_xy
        B[5 * ncol + c + 1] = dx;
    }
}

std::vector<double> element_stiffness(const Box& cell, const Config3D& cfg,
                                      const std::vector<Mode3D>& modes,
                                      const CellQuadrature3D& q) {
    const std::size_t nm   = modes.size();
    const std::size_t ncol = 3 * nm;
    const std::array<double, 36> C = cfg.elasticity();

    Vec3 jac{};
    double detJ = 1.0;
    for (int d = 0; d < 3; ++d) {
        const std::size_t i = static_cast<std::size_t>(d);
        jac[i] = 0.5 * (cell.hi[i] - cell.lo[i]);
        detJ  *= jac[i];
    }

    std::vector<double> Ke(ncol * ncol, 0.0);
    std::vector<double> B, CB(6 * ncol);

    for (std::size_t k = 0; k < q.xi.size(); ++k) {
        const ShapeValues3D s =
            shape_3d(cfg.p, modes, q.xi[k][0], q.xi[k][1], q.xi[k][2]);
        strain_displacement(s, jac, B);

        const double scale = q.w[k] * detJ * q.mat[k];

        for (std::size_t r = 0; r < 6; ++r)
            for (std::size_t c = 0; c < ncol; ++c) {
                double v = 0.0;
                for (std::size_t t = 0; t < 6; ++t)
                    v += C[r * 6 + t] * B[t * ncol + c];
                CB[r * ncol + c] = v * scale;
            }

        for (std::size_t i = 0; i < ncol; ++i)
            for (std::size_t j = 0; j < ncol; ++j) {
                double v = 0.0;
                for (std::size_t r = 0; r < 6; ++r)
                    v += B[r * ncol + i] * CB[r * ncol + j];
                Ke[i * ncol + j] += v;
            }
    }
    return Ke;
}

std::vector<double> element_force(const Box& cell, const Config3D& cfg,
                                  const std::vector<Mode3D>& modes,
                                  const CellQuadrature3D& q) {
    const std::size_t nm = modes.size();
    std::vector<double> Fe(3 * nm, 0.0);

    double detJ = 1.0;
    for (int d = 0; d < 3; ++d) {
        const std::size_t s = static_cast<std::size_t>(d);
        detJ *= 0.5 * (cell.hi[s] - cell.lo[s]);
    }

    for (std::size_t k = 0; k < q.xi.size(); ++k) {
        const Vec3 f = cfg.body_load(q.x[k]);
        if (f[0] == 0.0 && f[1] == 0.0 && f[2] == 0.0) continue;

        const ShapeValues3D s =
            shape_3d(cfg.p, modes, q.xi[k][0], q.xi[k][1], q.xi[k][2]);
        const double scale = q.w[k] * detJ * q.mat[k];

        for (std::size_t m = 0; m < nm; ++m)
            for (int d = 0; d < 3; ++d)
                Fe[3 * m + static_cast<std::size_t>(d)] +=
                    s.N[m] * f[static_cast<std::size_t>(d)] * scale;
    }
    return Fe;
}

DenseMatrix assemble_stiffness(const Mesh3D& m, const Config3D& cfg,
                               const std::vector<CellQuadrature3D>& quads) {
    const std::vector<Mode3D> modes = enumerate_modes(cfg.p);
    const std::size_t nm = modes.size();
    const std::size_t ncol = 3 * nm;

    DenseMatrix K(m.n_dof_total);
    for (std::size_t e = 0; e < m.cells.size(); ++e) {
        const std::vector<double> Ke =
            element_stiffness(m.cells[e].box, cfg, modes, quads[e]);
        const std::vector<int>& L = m.ltog[e];

        for (std::size_t a = 0; a < nm; ++a)
            for (int da = 0; da < 3; ++da) {
                const int ga = 3 * L[a] + da;
                for (std::size_t b = 0; b < nm; ++b)
                    for (int db = 0; db < 3; ++db) {
                        const int gb = 3 * L[b] + db;
                        K(ga, gb) += Ke[(3 * a + static_cast<std::size_t>(da)) * ncol
                                        + 3 * b + static_cast<std::size_t>(db)];
                    }
            }
    }
    return K;
}

std::vector<double> assemble_force(const Mesh3D& m, const Config3D& cfg,
                                   const std::vector<CellQuadrature3D>& quads) {
    const std::vector<Mode3D> modes = enumerate_modes(cfg.p);
    const std::size_t nm = modes.size();

    std::vector<double> F(static_cast<std::size_t>(m.n_dof_total), 0.0);
    for (std::size_t e = 0; e < m.cells.size(); ++e) {
        const std::vector<double> Fe =
            element_force(m.cells[e].box, cfg, modes, quads[e]);
        const std::vector<int>& L = m.ltog[e];
        for (std::size_t a = 0; a < nm; ++a)
            for (int d = 0; d < 3; ++d)
                F[static_cast<std::size_t>(3 * L[a] + d)] +=
                    Fe[3 * a + static_cast<std::size_t>(d)];
    }
    return F;
}

}  // namespace fcm