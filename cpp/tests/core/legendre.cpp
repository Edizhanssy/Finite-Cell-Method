//
// Created by Edizhan Yildiz on 06.08.26.
//

#include "fcm/core/legendre.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double kTol = 1e-13;

bool close_enough(double a, double b) {
    const double scale = std::max(1.0, std::max(std::fabs(a), std::fabs(b)));
    return std::fabs(a - b) <= kTol * scale;
}

std::ifstream open_or_die(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        std::fprintf(stderr, "cannot open %s\n", path.c_str());
        std::exit(2);
    }
    return in;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string dir = (argc > 1) ? argv[1] : std::string(FCM_REFERENCE_DIR);
    int checked = 0, failed = 0;
    std::string line;

    {   // legendre.txt : n xi P dP
        std::ifstream in = open_or_die(dir + "/legendre.txt");
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            int n; double xi, P_ref, dP_ref;
            ss >> n >> xi >> P_ref >> dP_ref;
            const fcm::LegendreValues v = fcm::legendre_and_derivs(n, xi);
            ++checked;
            if (!close_enough(v.P[n], P_ref) || !close_enough(v.dP[n], dP_ref)) {
                if (++failed <= 5)
                    std::printf("FAIL legendre n=%2d xi=%6.3f  P %.17g vs %.17g   "
                                "dP %.17g vs %.17g\n",
                                n, xi, v.P[n], P_ref, v.dP[n], dP_ref);
            }
        }
    }

    {   // shapefunc.txt : p xi i N dN
        std::ifstream in = open_or_die(dir + "/shapefunc.txt");
        std::vector<double> N, dN;
        int cached_p = -1;
        double cached_xi = 1e300;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            int p, i; double xi, N_ref, dN_ref;
            ss >> p >> xi >> i >> N_ref >> dN_ref;
            if (p != cached_p || xi != cached_xi) {
                N  = fcm::shape_functions(p, xi);
                dN = fcm::shape_function_derivs(p, xi);
                cached_p = p;
                cached_xi = xi;
            }
            ++checked;
            const std::size_t k = static_cast<std::size_t>(i);
            if (!close_enough(N[k], N_ref) || !close_enough(dN[k], dN_ref)) {
                if (++failed <= 5)
                    std::printf("FAIL shapefunc p=%2d xi=%6.3f i=%2d  N %.17g vs %.17g   "
                                "dN %.17g vs %.17g\n",
                                p, xi, i, N[k], N_ref, dN[k], dN_ref);
            }
        }
    }
    if (checked < 300) {
        std::printf("FAILED  legendre: only %d checks parsed (expected 320)\n", checked);
        return 2;
    }
    std::printf("%s  legendre: %d checks, %d failures\n",
                failed ? "FAILED" : "PASSED", checked, failed);
    return failed ? 1 : 0;
}