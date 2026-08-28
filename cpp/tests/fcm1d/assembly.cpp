#include "fcm/fcm1d/assembly.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double kRelTol = 1e-11;
double g_max_dev = 0.0;

bool close_rel(double a, double b) {
    const double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
    const double dev = std::fabs(a - b) / scale;
    g_max_dev = std::max(g_max_dev, dev);
    return dev <= kRelTol;
}

}

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : std::string(FCM_REFERENCE_DIR);
    const fcm::Config cfg;
    const fcm::Mesh   m = fcm::build_mesh(cfg);

    std::vector<fcm::ElementQuadrature> quads;
    for (const fcm::Element1D& el : m.elements)
        quads.push_back(fcm::build_element_quadrature(el, cfg));

    const int nm = cfg.n_modes();
    std::vector<std::vector<double>> Ke;
    for (std::size_t e = 0; e < m.elements.size(); ++e)
        Ke.push_back(fcm::element_stiffness(m.elements[e], cfg, quads[e]));

    const std::vector<double> F = fcm::assemble_force(m, cfg, quads);

    int checked = 0, failed = 0;
    std::string line;

    {   // stiffness_1d.txt : element i j Ke
        std::ifstream in(dir + "/stiffness_1d.txt");
        if (!in) { std::fprintf(stderr, "cannot open stiffness_1d.txt\n"); return 2; }
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            int e, i, j; double ref;
            ss >> e >> i >> j >> ref;
            ++checked;
            const double got = Ke[static_cast<std::size_t>(e)]
                                 [static_cast<std::size_t>(i) * nm + j];
            if (!close_rel(got, ref) && ++failed <= 5)
                std::printf("FAIL Ke e=%d i=%2d j=%2d  %.17g vs %.17g\n", e, i, j, got, ref);
        }
    }

    {   // force_1d.txt : dof F
        std::ifstream in(dir + "/force_1d.txt");
        if (!in) { std::fprintf(stderr, "cannot open force_1d.txt\n"); return 2; }
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            int d; double ref;
            ss >> d >> ref;
            ++checked;
            if (!close_rel(F[static_cast<std::size_t>(d)], ref) && ++failed <= 5)
                std::printf("FAIL F dof=%2d  %.17g vs %.17g\n",
                            d, F[static_cast<std::size_t>(d)], ref);
        }
    }

    for (std::size_t e = 0; e < Ke.size(); ++e) {
        double scale = 0.0;
        for (double v : Ke[e]) scale = std::max(scale, std::fabs(v));

        for (int i = 0; i < nm; ++i) {
            for (int j = i + 1; j < nm; ++j) {
                ++checked;
                const double a = Ke[e][static_cast<std::size_t>(i) * nm + j];
                const double b = Ke[e][static_cast<std::size_t>(j) * nm + i];
                if (std::fabs(a - b) > 1e-12 * scale && ++failed <= 5)
                    std::printf("FAIL symmetry e=%zu (%d,%d): %.17g vs %.17g\n", e, i, j, a, b);
            }
            ++checked;
            const double rb = Ke[e][static_cast<std::size_t>(i) * nm + 0] +
                              Ke[e][static_cast<std::size_t>(i) * nm + 1];
            if (std::fabs(rb) > 1e-12 * scale && ++failed <= 5)
                std::printf("FAIL rigid body e=%zu row %d: %.17g\n", e, i, rb);
        }
    }

    ++checked;
    const double fsum = F[0] + F[1] + F[2];
    if (std::fabs(fsum) > 1e-8 && ++failed <= 5)
        std::printf("FAIL total nodal force %.17g\n", fsum);

    std::printf("max relative difference: %.3e\n", g_max_dev);
    std::printf("sum of nodel forces: %.6e\n", fsum);
    if (checked < 500) { std::printf("FAILED  assembly: only %d checks\n", checked); return 2; }
    std::printf("%s  assembly: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}