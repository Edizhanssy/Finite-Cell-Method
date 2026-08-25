#pragma once
#include <cstddef>
#include <vector>

namespace fcm {

struct DenseMatrix {
    int n = 0;
    std::vector<double> a;

    DenseMatrix() = default;
    explicit DenseMatrix(int n_)
        : n(n_), a(static_cast<std::size_t>(n_) * n_, 0.0) {}

    double& operator()(int i, int j)       { return a[static_cast<std::size_t>(i) * n + j]; }
    double  operator()(int i, int j) const { return a[static_cast<std::size_t>(i) * n + j]; }
};

}  // namespace fcm