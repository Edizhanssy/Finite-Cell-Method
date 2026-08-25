#include "fcm/fcm3d/shape.hpp"

#include <stdexcept>

#include "fcm/core/legendre.hpp"

namespace fcm {

    std::vector<Mode3D> enumerate_modes(int p) {
        if (p < 1) throw std::invalid_argument("enumerate_modes: p < 1");

        std::vector<Mode3D> group[4];
        for (int i = 0; i <= p; ++i)
            for (int j = 0; j <= p; ++j)
                for (int k = 0; k <= p; ++k) {
                    const int internal = (i >= 2) + (j >= 2) + (k >= 2);
                    group[internal].push_back(
                        Mode3D{i, j, k, static_cast<ModeKind>(internal)});
                }

        std::vector<Mode3D> all;
        for (int g = 0; g < 4; ++g)
            all.insert(all.end(), group[g].begin(), group[g].end());
        return all;
    }

    ShapeValues3D shape_3d(int p, const std::vector<Mode3D>& modes,
                           double xi, double eta, double zeta) {
        const std::vector<double> Nx = shape_functions(p, xi);
        const std::vector<double> Ny = shape_functions(p, eta);
        const std::vector<double> Nz = shape_functions(p, zeta);
        const std::vector<double> Dx = shape_function_derivs(p, xi);
        const std::vector<double> Dy = shape_function_derivs(p, eta);
        const std::vector<double> Dz = shape_function_derivs(p, zeta);

        ShapeValues3D s;
        const std::size_t n = modes.size();
        s.N.resize(n); s.dNdxi.resize(n); s.dNdeta.resize(n); s.dNdzeta.resize(n);

        for (std::size_t m = 0; m < n; ++m) {
            const std::size_t i = static_cast<std::size_t>(modes[m].i);
            const std::size_t j = static_cast<std::size_t>(modes[m].j);
            const std::size_t k = static_cast<std::size_t>(modes[m].k);
            s.N[m]       = Nx[i] * Ny[j] * Nz[k];
            s.dNdxi[m]   = Dx[i] * Ny[j] * Nz[k];
            s.dNdeta[m]  = Nx[i] * Dy[j] * Nz[k];
            s.dNdzeta[m] = Nx[i] * Ny[j] * Dz[k];
        }
        return s;
    }

}  // namespace fcm