#include "fcm/fcm3d/quadrature.hpp"

#include "fcm/core/gauss.hpp"

namespace fcm {

    CellQuadrature3D build_cell_quadrature(const Box& cell, const Config3D& cfg) {
        const GaussRule g = gauss_legendre(cfg.n_gauss());
        const std::size_t n = g.points.size();

        CellQuadrature3D q;
        const std::size_t total = n * n * n;
        q.xi.reserve(total); q.w.reserve(total);
        q.x.reserve(total);  q.mat.reserve(total);

        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = 0; b < n; ++b)
                for (std::size_t c = 0; c < n; ++c) {
                    const Vec3 xi{g.points[a], g.points[b], g.points[c]};
                    const Vec3 xg = map_to_cell(xi, cell);
                    q.xi.push_back(xi);
                    q.w.push_back(g.weights[a] * g.weights[b] * g.weights[c]);
                    q.x.push_back(xg);
                    q.mat.push_back(cfg.material_factor(xg));
                }
        return q;
    }

}  // namespace fcm