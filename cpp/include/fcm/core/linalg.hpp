#pragma once
#include <vector>

#include "fcm/core/dense_matrix.hpp"

namespace fcm {

std::vector<double> solve_dense(DenseMatrix A, std::vector<double> b);

}  // namespace fcm