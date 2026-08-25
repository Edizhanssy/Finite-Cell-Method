#include "fcm/fcm1d/mesh.hpp"
#include "fcm/fcm1d/partition.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : std::string(FCM_REFERENCE_DIR);
    const fcm::Config cfg;
    const fcm::Mesh m = fcm::build_mesh(cfg);

    std::vector<std::vector<std::array<double, 2>>> subs;
    for (const fcm::Element1D& el : m.elements)
        subs.push_back(fcm::partition_element(el, cfg));

    int checked = 0, failed = 0;

    std::ifstream in(dir + "/subdomains_1d.txt");
    if (!in) { std::fprintf(stderr, "cannot open subdomains_1d.txt\n"); return 2; }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int e, i; double a_ref, b_ref;
        ss >> e >> i >> a_ref >> b_ref;
        ++checked;
        const auto& s = subs[static_cast<std::size_t>(e)][static_cast<std::size_t>(i)];
        if ((std::fabs(s[0] - a_ref) > 1e-15 || std::fabs(s[1] - b_ref) > 1e-15) && ++failed <= 5)
            std::printf("FAIL sub e=%d i=%2d  [%.17g, %.17g] vs [%.17g, %.17g]\n",
                        e, i, s[0], s[1], a_ref, b_ref);
    }

    for (std::size_t e = 0; e < subs.size(); ++e) {
        ++checked;
        if (static_cast<int>(subs[e].size()) != cfg.max_depth + 1 && ++failed <= 5)
            std::printf("FAIL element %zu: %zu sub-domain, expected %d\n",
                        e, subs[e].size(), cfg.max_depth + 1);
        ++checked;
        if (std::fabs(subs[e].front()[0] + 1.0) > 1e-15 ||
            std::fabs(subs[e].back()[1] - 1.0) > 1e-15) {
            if (++failed <= 5) std::printf("FAIL element %zu: boundaries are not [-1,1]\n", e);
        }
        for (std::size_t i = 1; i < subs[e].size(); ++i) {
            ++checked;
            if (std::fabs(subs[e][i][0] - subs[e][i - 1][1]) > 1e-15 && ++failed <= 5)
                std::printf("FAIL element %zu: %zu ile %zu the gap between \n", e, i - 1, i);
        }
    }

    if (checked < 20) { std::printf("FAILED  partition: only %d checks\n", checked); return 2; }
    std::printf("%s  partition: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}