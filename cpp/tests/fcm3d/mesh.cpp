#include "fcm/fcm3d/mesh.hpp"

#include <cstdio>
#include <map>
#include <set>
#include <vector>

namespace {

int checked = 0, failed = 0;

void check(bool ok, const char* what, long a = 0, long b = 0) {
    ++checked;
    if (!ok && ++failed <= 10) std::printf("FAIL %-40s  %ld vs %ld\n", what, a, b);
}

std::map<std::array<int,3>, std::size_t> mode_index(const std::vector<fcm::Mode3D>& modes) {
    std::map<std::array<int,3>, std::size_t> mi;
    for (std::size_t m = 0; m < modes.size(); ++m)
        mi[{modes[m].i, modes[m].j, modes[m].k}] = m;
    return mi;
}

/// d ekseninde komsu iki hucrenin paylastigi DOF'lari dogrular.
void check_sharing(const fcm::Mesh3D& msh, const std::vector<fcm::Mode3D>& modes,
                   int cellA, int cellB, int d, const char* label) {
    const auto mi = mode_index(modes);
    const int p = msh.p, q = p - 1;
    const std::vector<int>& A = msh.ltog[static_cast<std::size_t>(cellA)];
    const std::vector<int>& B = msh.ltog[static_cast<std::size_t>(cellB)];

    // A'nin +1 yuzu (d ekseninde indeks 1) <-> B'nin -1 yuzu (indeks 0)
    int matched = 0;
    for (int u = 0; u <= p; ++u)
        for (int v = 0; v <= p; ++v) {
            std::array<int,3> ka{}, kb{};
            const int o1 = (d + 1) % 3, o2 = (d + 2) % 3;
            ka[static_cast<std::size_t>(d)]  = 1;
            kb[static_cast<std::size_t>(d)]  = 0;
            ka[static_cast<std::size_t>(o1)] = kb[static_cast<std::size_t>(o1)] = u;
            ka[static_cast<std::size_t>(o2)] = kb[static_cast<std::size_t>(o2)] = v;
            const std::size_t ma = mi.at(ka), mb = mi.at(kb);
            ++checked;
            if (A[ma] != B[mb] && ++failed <= 10)
                std::printf("FAIL %s paylasim (%d,%d,%d)|(%d,%d,%d)  %d vs %d\n",
                            label, ka[0],ka[1],ka[2], kb[0],kb[1],kb[2], A[ma], B[mb]);
            ++matched;
        }
    check(matched == (p + 1) * (p + 1), "paylasilan yuz modu sayisi", matched, (p+1)*(p+1));

    // Toplam paylasim: 4 dugum + 4q kenar + q^2 yuz
    std::set<int> sa(A.begin(), A.end()), sb(B.begin(), B.end()), inter;
    for (int x : sa) if (sb.count(x)) inter.insert(x);
    const long expect = 4 + 4L * q + static_cast<long>(q) * q;
    check(static_cast<long>(inter.size()) == expect, "paylasilan toplam DOF",
          static_cast<long>(inter.size()), expect);
}

}  // namespace

int main() {
    // ---- 1) tek hucre: n_dof == (p+1)^3 ----------------------------------
    for (int p = 2; p <= 5; ++p) {
        fcm::Config3D cfg;
        cfg.p = p;
        cfg.n_elements = {1, 1, 1};
        const fcm::Mesh3D msh = fcm::build_mesh(cfg);
        check(msh.n_dof == msh.n_modes, "tek hucre n_dof = n_modes",
              msh.n_dof, msh.n_modes);
    }

    // ---- 2) elle hesaplanmis kucuk ornek ----------------------------------
    {
        fcm::Config3D cfg;
        cfg.p = 2;
        cfg.n_elements = {2, 1, 1};
        const fcm::Mesh3D msh = fcm::build_mesh(cfg);
        // 12 dugum + 20 kenar + 11 yuz + 2 hacim = 45
        check(msh.n_dof == 45, "2x1x1 p=2 n_dof", msh.n_dof, 45);
        check(msh.n_dof_total == 135, "vektor DOF", msh.n_dof_total, 135);
    }

    // ---- 3) genel izgara: aralik, teklik, ortusme, paylasim ----------------
    for (int p = 2; p <= 4; ++p) {
        fcm::Config3D cfg;
        cfg.p = p;
        cfg.n_elements = {2, 2, 2};
        const fcm::Mesh3D msh = fcm::build_mesh(cfg);
        const std::vector<fcm::Mode3D> modes = fcm::enumerate_modes(p);

        check(msh.cells.size() == 8, "hucre sayisi",
              static_cast<long>(msh.cells.size()), 8);

        std::set<int> all;
        bool in_range = true, unique_in_cell = true;
        for (const std::vector<int>& L : msh.ltog) {
            std::set<int> s;
            for (int d : L) {
                if (d < 0 || d >= msh.n_dof) in_range = false;
                if (!s.insert(d).second) unique_in_cell = false;
                all.insert(d);
            }
            if (static_cast<int>(s.size()) != msh.n_modes) unique_in_cell = false;
        }
        check(in_range,       "tum DOF aralikta");
        check(unique_in_cell, "hucre icinde tekrar yok");
        check(static_cast<long>(all.size()) == msh.n_dof, "her DOF en az bir hucrede",
              static_cast<long>(all.size()), msh.n_dof);

        // hucre sirasi (i*ny+j)*nz+k, ny=nz=2
        check_sharing(msh, modes, 0, 4, 0, "x");   // (0,0,0)-(1,0,0)
        check_sharing(msh, modes, 0, 2, 1, "y");   // (0,0,0)-(0,1,0)
        check_sharing(msh, modes, 0, 1, 2, "z");   // (0,0,0)-(0,0,1)
    }

    if (checked < 100) { std::printf("FAILED  mesh3d: only %d checks\n", checked); return 2; }
    std::printf("%s  mesh3d: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}