#pragma once

#include <cstddef>
#include <span>
#include <vector>

namespace zfit::detail {

struct ResidualWhitener {
    std::size_t n { };
    std::vector<double> lower;

    ResidualWhitener(std::size_t n, std::span<const double> covariance);

    void whiten(std::span<const double> residual, std::span<double> output) const;
};
}
