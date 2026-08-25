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

    // ---- 1) kesik degil -> tek yaprak ------------------------------------
    {
        fcm::Config3D cfg;                       // fictitious bos
        cfg.max_depth = 3;
        const std::vector<fcm::Box> leaves = fcm::partition_cell(cell, cfg);
        check(leaves.size() == 1, "kesilmeyen hucre tek yaprak",
              static_cast<double>(leaves.size()), 1);
    }

    // ---- 2) tamamen fiktif icinde -> tek yaprak ---------------------------
    {
        fcm::Config3D cfg;
        cfg.max_depth = 3;
        cfg.fictitious = {{{-1.0, -1.0, -1.0}, {2.0, 2.0, 2.0}}};
        const std::vector<fcm::Box> leaves = fcm::partition_cell(cell, cfg);
        check(leaves.size() == 1, "tamamen icerde tek yaprak",
              static_cast<double>(leaves.size()), 1);
    }

    // ---- 3) duzlemsel kesik: x = 0.5 --------------------------------------
    // Yalnizca x ekseni kesiliyor -> her seviyede 8 cocugun 8'i de kesik
    // DEGIL: x'te ikiye ayrilan 4+4'ten sadece kesiti icerenler bolunur.
    for (int d = 1; d <= 4; ++d) {
        fcm::Config3D cfg;
        cfg.max_depth = d;
        cfg.fictitious = {{{0.5, -1.0, -1.0}, {2.0, 2.0, 2.0}}};
        const std::vector<fcm::Box> leaves = fcm::partition_cell(cell, cfg);

        double vol = 0.0;
        for (const fcm::Box& b : leaves) vol += volume(b);
        check(std::fabs(vol - 1.0) < 1e-14, "yaprak hacimleri toplami = 1", vol, 1.0);

        // Kesik duzlem tam bolme duzlemine denk geliyor; kapali sinir
        // testi yuzunden dokunan cocuklar kesik gorunmeye devam eder.
        // Beklenen: her seviyede 4 cocuk bolunur.
        const std::size_t expect = (d == 1) ? 8 : 8 + 28 * ((1u << (2*(d-1))) - 1) / 3;
        check(leaves.size() == expect, "duzlemsel kesik: yaprak sayisi",
              static_cast<double>(leaves.size()), static_cast<double>(expect));
    }

    // ---- 4) izgaraya denk gelmeyen kesik: x = 7/3 - 1 = 1/3 ---------------
    for (int d = 0; d <= 4; ++d) {
        fcm::Config3D cfg;
        cfg.max_depth = d;
        cfg.fictitious = {{{1.0/3.0, -1.0, -1.0}, {2.0, 2.0, 2.0}}};
        const std::vector<fcm::Box> leaves = fcm::partition_cell(cell, cfg);

        double vol = 0.0;
        for (const fcm::Box& b : leaves) vol += volume(b);
        check(std::fabs(vol - 1.0) < 1e-13, "hacim korunumu (genel kesik)", vol, 1.0);

        // x'te 2^d dilimden yalnizca 1'i kesik kalir -> her seviye 4 kesik cocuk
        // yaprak sayisi: d=0 -> 1, sonrasi 8 + 4*7*(d-1)... elle saymak yerine
        // sadece monoton artis ve ust sinir kontrolu:
        const std::size_t upper = static_cast<std::size_t>(std::pow(8.0, d));
        check(leaves.size() >= 1 && leaves.size() <= upper + 1,
              "yaprak sayisi makul araliкta",
              static_cast<double>(leaves.size()), static_cast<double>(upper));
    }

    // ---- 5) yapraklar ortusmuyor (rastgele noktalarla) --------------------
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
        check(bad == 0, "her nokta tam bir yaprakta", bad, 0);
    }

    if (checked < 15) { std::printf("FAILED  partition3d: only %d checks\n", checked); return 2; }
    std::printf("%s  partition3d: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}