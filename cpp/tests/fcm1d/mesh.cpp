#include "fcm/fcm1d/mesh.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : std::string(FCM_REFERENCE_DIR);
    const fcm::Config cfg;
    const fcm::Mesh m = fcm::build_mesh(cfg);

    int checked = 0, failed = 0;

    ++checked;
    if (m.n_dof != 31) {
        ++failed;
        std::printf("FAIL n_dof = %d, beklenen 31\n", m.n_dof);
    }

    std::ifstream in(dir + "/ltog_1d.txt");
    if (!in) { std::fprintf(stderr, "cannot open ltog_1d.txt\n"); return 2; }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int e, i, dof_ref;
        ss >> e >> i >> dof_ref;
        ++checked;
        const int got = m.ltog[static_cast<std::size_t>(e)][static_cast<std::size_t>(i)];
        if (got != dof_ref && ++failed <= 5)
            std::printf("FAIL ltog e=%d i=%2d  %d vs %d\n", e, i, got, dof_ref);
    }

    const double expect[2][2] = {{0.0, 1.5}, {1.5, 3.0}};
    for (int e = 0; e < 2; ++e)
        for (int k = 0; k < 2; ++k) {
            ++checked;
            const double got = m.elements[static_cast<std::size_t>(e)].global_[static_cast<std::size_t>(k)];
            if (std::fabs(got - expect[e][k]) > 1e-15 && ++failed <= 5)
                std::printf("FAIL element %d coord %d: %.17g vs %.17g\n", e, k, got, expect[e][k]);
        }

    if (checked < 30) { std::printf("FAILED  mesh: only %d checks\n", checked); return 2; }
    std::printf("%s  mesh: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}