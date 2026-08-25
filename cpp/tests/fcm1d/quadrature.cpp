#include "fcm/fcm1d/quadrature.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
namespace {

constexpr double kRelTol = 1e-12;

bool close_rel(double a, double b) {
    const double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
    return std::fabs(a - b) <= kRelTol * scale;
}

}

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : std::string(FCM_REFERENCE_DIR);
    const fcm::Config cfg;
    const fcm::Mesh m = fcm::build_mesh(cfg);

    std::vector<fcm::ElementQuadrature> quads;
    for (const fcm::Element1D& el : m.elements)
        quads.push_back(fcm::build_element_quadrature(el, cfg));

    int checked = 0, failed = 0, mat_mismatch = 0;

    std::ifstream in(dir + "/quadrature_1d.txt");
    if (!in) { std::fprintf(stderr, "cannot open quadrature_1d.txt\n"); return 2; }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int e, i; double xi_r, w_r, x_r, mat_r;
        ss >> e >> i >> xi_r >> w_r >> x_r >> mat_r;
        const fcm::ElementQuadrature& q = quads[static_cast<std::size_t>(e)];
        const std::size_t k = static_cast<std::size_t>(i);
        ++checked;
        bool bad = !close_rel(q.xi[k], xi_r) || !close_rel(q.w[k], w_r) ||
                   !close_rel(q.x[k], x_r);
        if (q.mat[k] != mat_r) { bad = true; ++mat_mismatch; }
        if (bad && ++failed <= 5)
            std::printf("FAIL quad e=%d i=%3d\n  xi  %.17g vs %.17g\n  w   %.17g vs %.17g\n"
                        "  x   %.17g vs %.17g\n  mat %.17g vs %.17g\n",
                        e, i, q.xi[k], xi_r, q.w[k], w_r, q.x[k], x_r, q.mat[k], mat_r);
    }

    for (std::size_t e = 0; e < quads.size(); ++e) {
        const double sum = std::accumulate(quads[e].w.begin(), quads[e].w.end(), 0.0);
        const double len = std::fabs(m.elements[e].global_[1] - m.elements[e].global_[0]);
        ++checked;
        if (std::fabs(sum - len) > 1e-13 * len && ++failed <= 5)
            std::printf("FAIL element %zu: sum(w) = %.17g, length = %.17g\n", e, sum, len);

        ++checked;
        const std::size_t expect = static_cast<std::size_t>(quads[e].n_subdomains) *
                                   static_cast<std::size_t>(cfg.n_gauss());
        if (quads[e].xi.size() != expect && ++failed <= 5)
            std::printf("FAIL element %zu: %zu point, expected %zu\n",
                        e, quads[e].xi.size(), expect);
    }

    if (mat_mismatch)
        std::printf("NOT: %d mismatching materials - one quadrature point domain could"
                    "fall into near boundary\n", mat_mismatch);
    if (checked < 300) { std::printf("FAILED  quadrature: only %d checks\n", checked); return 2; }
    std::printf("%s  quadrature: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}