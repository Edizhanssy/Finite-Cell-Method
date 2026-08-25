#include "fcm/fcm1d/solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : std::string(FCM_REFERENCE_DIR);
    const fcm::Config cfg;
    const fcm::SolveResult r = fcm::solve(cfg);

    int checked = 0, failed = 0;
    double max_dev = 0.0;

    std::ifstream in(dir + "/solution_1d.txt");
    if (!in) { std::fprintf(stderr, "cannot open solution_1d.txt\n"); return 2; }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int d; double ref;
        ss >> d >> ref;
        ++checked;
        const double got = r.u[static_cast<std::size_t>(d)];
        const double dev = std::fabs(got - ref) / std::max(1.0, std::fabs(ref));
        max_dev = std::max(max_dev, dev);
        if (dev > 1e-7 && ++failed <= 5)
            std::printf("FAIL u[%2d]  %.17g vs %.17g   dev %.3e\n", d, got, ref, dev);
    }

    double res_inf = 0.0, f_inf = 0.0;
    for (int i = 0; i < r.K.n; ++i) {
        double s = 0.0;
        for (int j = 0; j < r.K.n; ++j) s += r.K(i, j) * r.u[static_cast<std::size_t>(j)];
        res_inf = std::max(res_inf, std::fabs(s - r.F[static_cast<std::size_t>(i)]));
        f_inf   = std::max(f_inf,   std::fabs(r.F[static_cast<std::size_t>(i)]));
    }
    ++checked;
    if (res_inf / f_inf > 1e-12 && ++failed <= 5)
        std::printf("FAIL artik %.3e\n", res_inf / f_inf);

    std::printf("max realative difference: %.3e   "
                "(kappa*eps limit ~1.7e-03)\n", max_dev);
    std::printf("relative fraction ||Ku-F||/||F||:       %.3e\n", res_inf / f_inf);
    std::printf("total nodal force:             %.6e\n", r.nodal_force_sum);
    std::printf("%s  solve: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}