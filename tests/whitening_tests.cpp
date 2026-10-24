#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/whitening.hpp"

TEST_CASE("whitening divides independent residuals by their standard deviations")
{
    const std::vector<double> covariance{4.0, 0.0, 0.0, 9.0};
    const zfit::detail::ResidualWhitener whitener{2, covariance};
    const std::vector<double> residual{4.0, 6.0};
    std::vector<double> whitened(2);

    whitener.whiten(residual, whitened);

    REQUIRE(whitened[0] == Catch::Approx(2.0));
    REQUIRE(whitened[1] == Catch::Approx(2.0));
    REQUIRE(whitener.lower.size() == 4);
    REQUIRE(whitener.lower[0] == Catch::Approx(2.0));
    REQUIRE(whitener.lower[2] == Catch::Approx(0.0));
    REQUIRE(whitener.lower[3] == Catch::Approx(3.0));
}

TEST_CASE("whitening accounts for correlated measurement residuals")
{
    // Sigma = [[4, 2], [2, 3]], with L = [[2, 0], [1, sqrt(2)]].
    const std::vector<double> covariance{4.0, 2.0, 2.0, 3.0};
    const zfit::detail::ResidualWhitener whitener{2, covariance};
    const std::vector<double> residual{2.0, 3.0};
    std::vector<double> whitened(2);

    whitener.whiten(residual, whitened);

    REQUIRE(whitened[0] == Catch::Approx(1.0));
    REQUIRE(whitened[1] == Catch::Approx(std::sqrt(2.0)));
    const double chi_squared = whitened[0] * whitened[0] + whitened[1] * whitened[1];
    REQUIRE(chi_squared == Catch::Approx(3.0));
    REQUIRE(whitener.lower[2] * whitened[0] + whitener.lower[3] * whitened[1] ==
            Catch::Approx(residual[1]));
}

TEST_CASE("whitening supports the same vector as residual and output")
{
    const std::vector<double> covariance{4.0, 2.0, 2.0, 3.0};
    const zfit::detail::ResidualWhitener whitener{2, covariance};
    std::vector<double> residual{2.0, 3.0};

    whitener.whiten(residual, residual);

    REQUIRE(residual[0] == Catch::Approx(1.0));
    REQUIRE(residual[1] == Catch::Approx(std::sqrt(2.0)));
}

TEST_CASE("whitening rejects singular, indefinite, or asymmetric covariance")
{
    const std::vector<double> singular{1.0, 1.0, 1.0, 1.0};
    REQUIRE_THROWS_AS(zfit::detail::ResidualWhitener(2, singular), std::invalid_argument);

    const std::vector<double> indefinite{1.0, 2.0, 2.0, 1.0};
    REQUIRE_THROWS_AS(zfit::detail::ResidualWhitener(2, indefinite), std::invalid_argument);

    const std::vector<double> asymmetric{1.0, 0.1, 0.2, 1.0};
    REQUIRE_THROWS_AS(zfit::detail::ResidualWhitener(2, asymmetric), std::invalid_argument);

    const std::vector<double> nonfinite{1.0, 0.0, 0.0, std::numeric_limits<double>::infinity()};
    REQUIRE_THROWS_AS(zfit::detail::ResidualWhitener(2, nonfinite), std::invalid_argument);
}

TEST_CASE("whitening rejects invalid dimensions and a mutated factor")
{
    const std::vector<double> covariance{4.0, 0.0, 0.0, 9.0};
    REQUIRE_THROWS_AS(zfit::detail::ResidualWhitener(0, covariance), std::invalid_argument);
    REQUIRE_THROWS_AS(zfit::detail::ResidualWhitener(3, covariance), std::invalid_argument);

    zfit::detail::ResidualWhitener whitener{2, covariance};
    const std::vector<double> residual{4.0, 6.0};
    std::vector<double> short_output(1);
    REQUIRE_THROWS_AS(whitener.whiten(residual, short_output), std::invalid_argument);

    std::vector<double> output(2);
    const std::vector<double> nonfinite_residual{4.0, std::numeric_limits<double>::quiet_NaN()};
    REQUIRE_THROWS_AS(whitener.whiten(nonfinite_residual, output), std::invalid_argument);

    whitener.lower[3] = 0.0;
    REQUIRE_THROWS_AS(whitener.whiten(residual, output), std::invalid_argument);
}
