//
// Created by Edizhan Yildiz on 06.08.26.
//

#include "fcm/core/gauss.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : std::string(FCM_REFERENCE_DIR);
    int checked = 0, failed = 0;

    {
        std::ifstream in(dir + "/gauss.txt");
        if (!in) { std::fprintf(stderr, "cannot open gauss.txt\n"); return 2; }

        std::string line;
        fcm::GaussRule rule;
        int cached_n = -1;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            int n, i; double p_ref, w_ref;
            ss >> n >> i >> p_ref >> w_ref;
            if (n != cached_n) { rule = fcm::gauss_legendre(n); cached_n = n; }
            ++checked;
            const std::size_t k = static_cast<std::size_t>(i);
            if (std::fabs(rule.points[k]  - p_ref) > 1e-12 ||
                std::fabs(rule.weights[k] - w_ref) > 1e-12) {
                if (++failed <= 5)
                    std::printf("FAIL table n=%2d i=%2d  x %.17g vs %.17g   w %.17g vs %.17g\n",
                                n, i, rule.points[k], p_ref, rule.weights[k], w_ref);
            }
        }
    }

    for (int n = 1; n <= 40; ++n) {
        const fcm::GaussRule rule = fcm::gauss_legendre(n);
        for (int k = 0; k <= 2 * n - 1; ++k) {
            double sum = 0.0;
            for (std::size_t j = 0; j < rule.points.size(); ++j)
                sum += rule.weights[j] * std::pow(rule.points[j], k);
            const double exact = (k % 2 == 0) ? 2.0 / (k + 1) : 0.0;
            ++checked;
            if (std::fabs(sum - exact) > 1e-11 * std::max(1.0, std::fabs(exact))) {
                if (++failed <= 5)
                    std::printf("FAIL moment n=%2d k=%2d  %.17g vs %.17g\n",
                                n, k, sum, exact);
            }
        }
    }

    if (checked < 1000) {
        std::printf("FAILED  gauss: only %d checks\n", checked);
        return 2;
    }
    std::printf("%s  gauss: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}