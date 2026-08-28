#include "fcm/fcm3d/partition.hpp"

namespace fcm {
    namespace {

        constexpr int kSeed = 5;      // for each axis: 5^3 = 125 seeds

        bool is_cut(const Box& b, const Config3D& cfg) {
            bool first = false, init = false;
            for (int a = 0; a < kSeed; ++a)
                for (int c = 0; c < kSeed; ++c)
                    for (int e = 0; e < kSeed; ++e) {
                        const int t[3] = {a, c, e};
                        Vec3 x{};
                        for (int d = 0; d < 3; ++d) {
                            const std::size_t s = static_cast<std::size_t>(d);
                            const double u = static_cast<double>(t[d]) / (kSeed - 1);
                            x[s] = b.lo[s] + u * (b.hi[s] - b.lo[s]);
                        }
                        const bool in = cfg.inside_fictitious(x);
                        if (!init) { first = in; init = true; }
                        else if (in != first) return true;
                    }
            return false;
        }

        void recurse(const Box& b, const Config3D& cfg, int depth, std::vector<Box>& out) {
            if (depth < cfg.max_depth && is_cut(b, cfg)) {
                Vec3 mid{};
                for (int d = 0; d < 3; ++d) {
                    const std::size_t s = static_cast<std::size_t>(d);
                    mid[s] = 0.5 * (b.lo[s] + b.hi[s]);
                }
                for (int oct = 0; oct < 8; ++oct) {
                    Box child{};
                    for (int d = 0; d < 3; ++d) {
                        const std::size_t s = static_cast<std::size_t>(d);
                        const bool upper = (oct >> d) & 1;
                        child.lo[s] = upper ? mid[s] : b.lo[s];
                        child.hi[s] = upper ? b.hi[s] : mid[s];
                    }
                    recurse(child, cfg, depth + 1, out);
                }
            } else {
                out.push_back(b);
            }
        }

    }  // namespace

    std::vector<Box> partition_cell(const Box& cell, const Config3D& cfg) {
        std::vector<Box> out;
        recurse(cell, cfg, 0, out);
        return out;
    }

}  // namespace fcm