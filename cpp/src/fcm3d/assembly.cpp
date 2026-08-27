#include "fcm/fcm3d/assembly.hpp"

#include <cstddef>
#include "fcm/core/legendre.hpp"
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
    const int p  = cfg.p;
    const int np = p + 1;
    const std::size_t ng  = static_cast<std::size_t>(q.ng);
    const std::size_t ng3 = ng * ng * ng;

    const double lam = cfg.E * cfg.nu / ((1.0 + cfg.nu) * (1.0 - 2.0 * cfg.nu));
    const double mu  = cfg.E / (2.0 * (1.0 + cfg.nu));

    Vec3 jac{}, inv{};
    double detJ = 1.0;
    for (int d = 0; d < 3; ++d) {
        const std::size_t i = static_cast<std::size_t>(d);
        jac[i] = 0.5 * (cell.hi[i] - cell.lo[i]);
        inv[i] = 1.0 / jac[i];
        detJ  *= jac[i];
    }

    std::vector<double> Ke(ncol * ncol, 0.0);
    std::vector<double> gx(nm), gy(nm), gz(nm);

    std::array<std::vector<double>, 3> Nt, Dt;
    for (int d = 0; d < 3; ++d) {
        const std::size_t s = static_cast<std::size_t>(d);
        Nt[s].resize(ng * static_cast<std::size_t>(np));
        Dt[s].resize(ng * static_cast<std::size_t>(np));
    }

    std::vector<int> mi(nm), mj(nm), mk(nm);
    for (std::size_t m = 0; m < nm; ++m) {
        mi[m] = modes[m].i; mj[m] = modes[m].j; mk[m] = modes[m].k;
    }

    for (std::size_t leaf = 0; leaf < q.leaf_xi1d.size(); ++leaf) {
        // yaprak basina 1D tablolar
        for (int d = 0; d < 3; ++d) {
            const std::size_t s = static_cast<std::size_t>(d);
            for (std::size_t a = 0; a < ng; ++a) {
                const double t = q.leaf_xi1d[leaf][s][a];
                const std::vector<double> N1 = shape_functions(p, t);
                const std::vector<double> D1 = shape_function_derivs(p, t);
                for (int i = 0; i < np; ++i) {
                    const std::size_t o = a * static_cast<std::size_t>(np)
                                        + static_cast<std::size_t>(i);
                    Nt[s][o] = N1[static_cast<std::size_t>(i)];
                    Dt[s][o] = D1[static_cast<std::size_t>(i)];
                }
            }
        }

        for (std::size_t a = 0; a < ng; ++a)
            for (std::size_t b = 0; b < ng; ++b)
                for (std::size_t c = 0; c < ng; ++c) {
                    const std::size_t k = leaf * ng3 + (a * ng + b) * ng + c;

                    const double* Na = &Nt[0][a * static_cast<std::size_t>(np)];
                    const double* Da = &Dt[0][a * static_cast<std::size_t>(np)];
                    const double* Nb = &Nt[1][b * static_cast<std::size_t>(np)];
                    const double* Db = &Dt[1][b * static_cast<std::size_t>(np)];
                    const double* Nc = &Nt[2][c * static_cast<std::size_t>(np)];
                    const double* Dc = &Dt[2][c * static_cast<std::size_t>(np)];

                    for (std::size_t m = 0; m < nm; ++m) {
                        const double nx = Na[mi[m]], ny = Nb[mj[m]], nz = Nc[mk[m]];
                        const double dx = Da[mi[m]], dy = Db[mj[m]], dz = Dc[mk[m]];
                        gx[m] = dx * ny * nz * inv[0];
                        gy[m] = nx * dy * nz * inv[1];
                        gz[m] = nx * ny * dz * inv[2];
                    }

                    const double w  = q.w[k] * detJ * q.mat[k];
                    const double wl = w * lam;
                    const double wm = w * mu;

                    for (std::size_t ia = 0; ia < nm; ++ia) {
                        const double ax = gx[ia], ay = gy[ia], az = gz[ia];
                        const std::size_t ra = 3 * ia;

                        for (std::size_t ib = ia; ib < nm; ++ib) {
                            const double bx = gx[ib], by = gy[ib], bz = gz[ib];
                            const double dot = wm * (ax * bx + ay * by + az * bz);
                            const std::size_t rb = 3 * ib;

                            const double blk[3][3] = {
                                {wl*ax*bx + wm*ax*bx + dot, wl*ax*by + wm*ay*bx, wl*ax*bz + wm*az*bx},
                                {wl*ay*bx + wm*ax*by, wl*ay*by + wm*ay*by + dot, wl*ay*bz + wm*az*by},
                                {wl*az*bx + wm*ax*bz, wl*az*by + wm*ay*bz, wl*az*bz + wm*az*bz + dot}};

                            for (int i = 0; i < 3; ++i)
                                for (int j = 0; j < 3; ++j) {
                                    const double v = blk[i][j];
                                    Ke[(ra + static_cast<std::size_t>(i)) * ncol
                                       + rb + static_cast<std::size_t>(j)] += v;
                                    if (ia != ib)
                                        Ke[(rb + static_cast<std::size_t>(j)) * ncol
                                           + ra + static_cast<std::size_t>(i)] += v;
                                }
                        }
                    }
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