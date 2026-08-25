#include "fcm/fcm3d/mesh.hpp"

#include <stdexcept>

namespace fcm {

Mesh3D build_mesh(const Config3D& cfg) {
    const int nx = cfg.n_elements[0];
    const int ny = cfg.n_elements[1];
    const int nz = cfg.n_elements[2];
    if (nx < 1 || ny < 1 || nz < 1)
        throw std::invalid_argument("build_mesh: n_elements < 1");
    if (cfg.p < 1) throw std::invalid_argument("build_mesh: p < 1");

    const int q = cfg.p - 1;                 // kenar basina mod
    const int qq = q * q, qqq = q * q * q;

    const int c_node = (nx + 1) * (ny + 1) * (nz + 1);
    const int c_ex   = nx * (ny + 1) * (nz + 1);
    const int c_ey   = (nx + 1) * ny * (nz + 1);
    const int c_ez   = (nx + 1) * (ny + 1) * nz;
    const int c_fx   = (nx + 1) * ny * nz;
    const int c_fy   = nx * (ny + 1) * nz;
    const int c_fz   = nx * ny * (nz + 1);
    const int c_vol  = nx * ny * nz;

    Mesh3D msh;
    msh.n       = cfg.n_elements;
    msh.p       = cfg.p;
    msh.n_modes = (cfg.p + 1) * (cfg.p + 1) * (cfg.p + 1);

    msh.off_node = 0;
    msh.off_edge = msh.off_node + c_node;
    msh.off_face = msh.off_edge + (c_ex + c_ey + c_ez) * q;
    msh.off_vol  = msh.off_face + (c_fx + c_fy + c_fz) * qq;
    msh.n_dof    = msh.off_vol  + c_vol * qqq;
    msh.entity_center.assign(static_cast<std::size_t>(msh.n_dof), Vec3{});
    msh.is_node.assign(static_cast<std::size_t>(msh.n_dof), 0);
    msh.n_dof_total = 3 * msh.n_dof;

    auto id_node = [&](int i, int j, int k) { return (i * (ny + 1) + j) * (nz + 1) + k; };
    auto id_ex   = [&](int i, int j, int k) { return (i * (ny + 1) + j) * (nz + 1) + k; };
    auto id_ey   = [&](int i, int j, int k) { return (i * ny + j) * (nz + 1) + k; };
    auto id_ez   = [&](int i, int j, int k) { return (i * (ny + 1) + j) * nz + k; };
    auto id_fx   = [&](int i, int j, int k) { return (i * ny + j) * nz + k; };
    auto id_fy   = [&](int i, int j, int k) { return (i * (ny + 1) + j) * nz + k; };
    auto id_fz   = [&](int i, int j, int k) { return (i * ny + j) * (nz + 1) + k; };
    auto id_vol  = [&](int i, int j, int k) { return (i * ny + j) * nz + k; };

    Vec3 h{};
    for (int d = 0; d < 3; ++d) {
        const std::size_t s = static_cast<std::size_t>(d);
        h[s] = (cfg.hi[s] - cfg.lo[s]) / cfg.n_elements[s];
    }

    const std::vector<Mode3D> modes = enumerate_modes(cfg.p);

    msh.cells.reserve(static_cast<std::size_t>(c_vol));
    msh.ltog.reserve(static_cast<std::size_t>(c_vol));

    for (int ci = 0; ci < nx; ++ci)
        for (int cj = 0; cj < ny; ++cj)
            for (int ck = 0; ck < nz; ++ck) {
                Cell3D cell{};
                cell.idx = {ci, cj, ck};
                const int gi[3] = {ci, cj, ck};
                for (int d = 0; d < 3; ++d) {
                    const std::size_t s = static_cast<std::size_t>(d);
                    cell.box.lo[s] = cfg.lo[s] + gi[d] * h[s];
                    cell.box.hi[s] = cell.box.lo[s] + h[s];
                }
                msh.cells.push_back(cell);

                std::vector<int> L(static_cast<std::size_t>(msh.n_modes), -1);
                for (std::size_t m = 0; m < modes.size(); ++m) {
                    const int a = modes[m].i, b = modes[m].j, c = modes[m].k;
                    const bool ia = a >= 2, ib = b >= 2, ic = c >= 2;
                    const int internal = ia + ib + ic;
                    int dof = -1;

                    if (internal == 0) {
                        dof = msh.off_node + id_node(ci + a, cj + b, ck + c);
                    } else if (internal == 1) {
                        if (ia)
                            dof = msh.off_edge + (id_ex(ci, cj + b, ck + c)) * q + (a - 2);
                        else if (ib)
                            dof = msh.off_edge + (c_ex + id_ey(ci + a, cj, ck + c)) * q + (b - 2);
                        else
                            dof = msh.off_edge + (c_ex + c_ey + id_ez(ci + a, cj + b, ck)) * q + (c - 2);
                    } else if (internal == 2) {
                        if (!ia)
                            dof = msh.off_face + (id_fx(ci + a, cj, ck)) * qq + (b - 2) * q + (c - 2);
                        else if (!ib)
                            dof = msh.off_face + (c_fx + id_fy(ci, cj + b, ck)) * qq + (a - 2) * q + (c - 2);
                        else
                            dof = msh.off_face + (c_fx + c_fy + id_fz(ci, cj, ck + c)) * qq + (a - 2) * q + (b - 2);
                    } else {
                        dof = msh.off_vol + id_vol(ci, cj, ck) * qqq
                            + ((a - 2) * q + (b - 2)) * q + (c - 2);
                    }
                    const int idx[3] = {a, b, c};
                    Vec3 ctr{};
                    bool nodal = true;
                    for (int d = 0; d < 3; ++d) {
                        const std::size_t s = static_cast<std::size_t>(d);
                        if (idx[d] < 2)
                            ctr[s] = (idx[d] == 0) ? cell.box.lo[s] : cell.box.hi[s];
                        else {
                            ctr[s] = 0.5 * (cell.box.lo[s] + cell.box.hi[s]);
                            nodal = false;
                        }
                    }
                    msh.entity_center[static_cast<std::size_t>(dof)] = ctr;
                    msh.is_node[static_cast<std::size_t>(dof)] = nodal ? 1 : 0;
                    L[m] = dof;
                }
                msh.ltog.push_back(std::move(L));
            }

    return msh;
}

}  // namespace fcm