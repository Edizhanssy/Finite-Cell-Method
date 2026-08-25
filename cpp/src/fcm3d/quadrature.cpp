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

        for (const Box& leaf : leaves) {
            // yaprak -> hucre referans olcegi
            double jleaf = 1.0;
            for (int d = 0; d < 3; ++d) {
                const std::size_t s = static_cast<std::size_t>(d);
                jleaf *= (leaf.hi[s] - leaf.lo[s]) / (cell.hi[s] - cell.lo[s]);
            }

            for (std::size_t a = 0; a < ng; ++a)
                for (std::size_t b = 0; b < ng; ++b)
                    for (std::size_t c = 0; c < ng; ++c) {
                        const Vec3 xg{
                            0.5 * ((1.0 - g.points[a]) * leaf.lo[0] + (1.0 + g.points[a]) * leaf.hi[0]),
                            0.5 * ((1.0 - g.points[b]) * leaf.lo[1] + (1.0 + g.points[b]) * leaf.hi[1]),
                            0.5 * ((1.0 - g.points[c]) * leaf.lo[2] + (1.0 + g.points[c]) * leaf.hi[2])};

                        const Vec3 xi{to_ref(xg[0], cell.lo[0], cell.hi[0]),
                                      to_ref(xg[1], cell.lo[1], cell.hi[1]),
                                      to_ref(xg[2], cell.lo[2], cell.hi[2])};

                        q.xi.push_back(xi);
                        q.w.push_back(g.weights[a] * g.weights[b] * g.weights[c] * jleaf);
                        q.x.push_back(xg);
                        q.mat.push_back(cfg.material_factor(xg));
                    }
        }
        return q;
    }

}  // namespace fcm