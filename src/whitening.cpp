#include "zfit/whitening.hpp"

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace zfit::detail {

ResidualWhitener::ResidualWhitener(std::size_t count, std::span<const double> covariance)
    : n(count)
{
    if (n == 0 || n > std::numeric_limits<std::size_t>::max() / n || covariance.size() != n * n
        || n > static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max())) {
        throw std::invalid_argument("invalid covariance dimensions");
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(covariance[i * n + i]) || covariance[i * n + i] <= 0) {
            throw std::invalid_argument("invalid covariance diagonal");
        }

        for (std::size_t j = 0; j < i; ++j) {
            const double a = covariance[i * n + j];
            const double b = covariance[j * n + i];
            if (!std::isfinite(a) || !std::isfinite(b)) {
                throw std::invalid_argument("nonfinite covariance");
            }

            const double scale = std::max({ std::abs(a), std::abs(b),
                std::sqrt(covariance[i * n + i]) * std::sqrt(covariance[j * n + j]) });
            if (std::abs(a - b) > 64.0 * std::numeric_limits<double>::epsilon() * scale) {
                throw std::invalid_argument("nonsymmetric covariance");
            }
        }
    }

    using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    const auto dim = static_cast<Eigen::Index>(n);
    const Eigen::Map<const RowMajorMatrix> matrix(covariance.data(), dim, dim);
    const Eigen::LLT<Eigen::MatrixXd> factor(matrix);

    if (factor.info() != Eigen::Success) {
        throw std::invalid_argument("covariance is not positive definite");
    }

    lower.assign(n * n, 0.0);
    const auto L = factor.matrixL();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            const double value = L(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
            if (!std::isfinite(value) || (i == j && value <= 0.0)) {
                throw std::invalid_argument("invalid covariance factor");
            }
            lower[i * n + j] = value;
        }
    }
}

void ResidualWhitener::whiten(std::span<const double> residual, std::span<double> output) const
{
    if (n == 0 || n > std::numeric_limits<std::size_t>::max() / n || lower.size() != n * n
        || residual.size() != n || output.size() != n) {
        throw std::invalid_argument("invalid whitener dimensions");
    }

    for (std::size_t i = 0; i < n; ++i) {
        double value = residual[i];
        if (!std::isfinite(value)) {
            throw std::invalid_argument("nonfinite residual");
        }

        for (std::size_t j = 0; j < i; ++j) {
            const double entry = lower[i * n + j];
            if (!std::isfinite(entry)) {
                throw std::invalid_argument("invalid covariance factor");
            }
            value -= entry * output[j];
        }

        const double diagonal = lower[i * n + i];
        if (!std::isfinite(diagonal) || diagonal <= 0.0) {
            throw std::invalid_argument("invalid covariance factor");
        }

        output[i] = value / diagonal;
        if (!std::isfinite(output[i])) {
            throw std::invalid_argument("nonfinite whitened residual");
        }
    }
}

}
