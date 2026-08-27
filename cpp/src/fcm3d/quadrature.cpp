#include "fcm/fcm3d/quadrature.hpp"

#include "fcm/core/gauss.hpp"
#include "fcm/fcm3d/partition.hpp"

namespace fcm {
    namespace {

        /// [lo,hi] araligini [-1,1]'e goturur.
        double to_ref(double x, double lo, double hi) {
            return 2.0 * (x - lo) / (hi - lo) - 1.0;
        }

    }  // namespace

    CellQuadrature3D build_cell_quadrature(const Box& cell, const Config3D& cfg) {
        const GaussRule g = gauss_legendre(cfg.n_gauss());
        const std::vector<Box> leaves = partition_cell(cell, cfg);
        const std::size_t ng = g.points.size();

        CellQuadrature3D q;
        const std::size_t total = leaves.size() * ng * ng * ng;
        q.xi.reserve(total); q.w.reserve(total);
        q.x.reserve(total);  q.mat.reserve(total);

        q.ng = static_cast<int>(ng);
        q.leaf_xi1d.reserve(leaves.size());

        for (const Box& leaf : leaves) {
            double jleaf = 1.0;
            for (int d = 0; d < 3; ++d) {
                const std::size_t s = static_cast<std::size_t>(d);
                jleaf *= (leaf.hi[s] - leaf.lo[s]) / (cell.hi[s] - cell.lo[s]);
            }

            std::array<std::vector<double>, 3> xg1d, xi1d;
            for (int d = 0; d < 3; ++d) {
                const std::size_t s = static_cast<std::size_t>(d);
                xg1d[s].resize(ng);
                xi1d[s].resize(ng);
                for (std::size_t a = 0; a < ng; ++a) {
                    const double xgv = 0.5 * ((1.0 - g.points[a]) * leaf.lo[s]
                                            + (1.0 + g.points[a]) * leaf.hi[s]);
                    xg1d[s][a] = xgv;
                    xi1d[s][a] = to_ref(xgv, cell.lo[s], cell.hi[s]);
                }
            }
            q.leaf_xi1d.push_back(xi1d);

            for (std::size_t a = 0; a < ng; ++a)
                for (std::size_t b = 0; b < ng; ++b)
                    for (std::size_t c = 0; c < ng; ++c) {
                        const Vec3 xg{xg1d[0][a], xg1d[1][b], xg1d[2][c]};
                        const Vec3 xi{xi1d[0][a], xi1d[1][b], xi1d[2][c]};
                        q.xi.push_back(xi);
                        q.w.push_back(g.weights[a] * g.weights[b] * g.weights[c] * jleaf);
                        q.x.push_back(xg);
                        q.mat.push_back(cfg.material_factor(xg));
                    }
        }
        return q;
    }

}  // namespace fcm