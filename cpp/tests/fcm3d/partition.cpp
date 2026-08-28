#include "fcm/fcm3d/partition.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

namespace {
int checked = 0, failed = 0;
void check(bool ok, const char* what, double a = 0, double b = 0) {
    ++checked;
    if (!ok && ++failed <= 10) std::printf("FAIL %-38s  %.10g vs %.10g\n", what, a, b);
}
double volume(const fcm::Box& b) {
    return (b.hi[0]-b.lo[0]) * (b.hi[1]-b.lo[1]) * (b.hi[2]-b.lo[2]);
}
}  // namespace

int main() {
    const fcm::Box cell{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};

    {
        fcm::Config3D cfg;                       // fictitious box
        cfg.max_depth = 3;
        const std::vector<fcm::Box> leaves = fcm::partition_cell(cell, cfg);
        check(leaves.size() == 1, "one leaf cell could not be partitioned!",
              static_cast<double>(leaves.size()), 1);
    }

    {
        fcm::Config3D cfg;
        cfg.max_depth = 3;
        cfg.fictitious = {{{-1.0, -1.0, -1.0}, {2.0, 2.0, 2.0}}};
        const std::vector<fcm::Box> leaves = fcm::partition_cell(cell, cfg);
        check(leaves.size() == 1, "one leaf is inside",
              static_cast<double>(leaves.size()), 1);
    }


    for (int d = 1; d <= 4; ++d) {
        fcm::Config3D cfg;
        cfg.max_depth = d;
        cfg.fictitious = {{{0.5, -1.0, -1.0}, {2.0, 2.0, 2.0}}};
        const std::vector<fcm::Box> leaves = fcm::partition_cell(cell, cfg);

        double vol = 0.0;
        for (const fcm::Box& b : leaves) vol += volume(b);
        check(std::fabs(vol - 1.0) < 1e-14, "yaprak hacimleri toplami = 1", vol, 1.0);

        const std::size_t expect = (d == 1) ? 8 : 8 + 28 * ((1u << (2*(d-1))) - 1) / 3;
        check(leaves.size() == expect, "longitudinal cut: total number of leaves",
              static_cast<double>(leaves.size()), static_cast<double>(expect));
    }

    for (int d = 0; d <= 4; ++d) {
        fcm::Config3D cfg;
        cfg.max_depth = d;
        cfg.fictitious = {{{1.0/3.0, -1.0, -1.0}, {2.0, 2.0, 2.0}}};
        const std::vector<fcm::Box> leaves = fcm::partition_cell(cell, cfg);

        double vol = 0.0;
        for (const fcm::Box& b : leaves) vol += volume(b);
        check(std::fabs(vol - 1.0) < 1e-13, "hacim korunumu (genel kesik)", vol, 1.0);

        const std::size_t upper = static_cast<std::size_t>(std::pow(8.0, d));
        check(leaves.size() >= 1 && leaves.size() <= upper + 1,
              "leave number is in acceptable range",
              static_cast<double>(leaves.size()), static_cast<double>(upper));
    }

    {
        fcm::Config3D cfg;
        cfg.max_depth = 3;
        cfg.fictitious = {{{1.0/3.0, -1.0, -1.0}, {2.0, 2.0, 2.0}}};
        const std::vector<fcm::Box> leaves = fcm::partition_cell(cell, cfg);

        unsigned seed = 12345;
        auto rnd = [&seed]() {
            seed = seed * 1103515245u + 12345u;
            return static_cast<double>((seed >> 16) & 0x7fff) / 32768.0;
        };
        int bad = 0;
        for (int t = 0; t < 2000; ++t) {
            const fcm::Vec3 x{rnd(), rnd(), rnd()};
            int hits = 0;
            for (const fcm::Box& b : leaves) {
                bool in = true;
                for (int dd = 0; dd < 3; ++dd) {
                    const std::size_t s = static_cast<std::size_t>(dd);
                    if (x[s] < b.lo[s] || x[s] >= b.hi[s]) { in = false; break; }
                }
                if (in) ++hits;
            }
            if (hits != 1) ++bad;
        }
        check(bad == 0, "every node is in one leaf", bad, 0);
    }

    if (checked < 15) { std::printf("FAILED  partition3d: only %d checks\n", checked); return 2; }
    std::printf("%s  partition3d: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}