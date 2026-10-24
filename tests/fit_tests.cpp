#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/cvodes_integrator.hpp"
#include "zfit/fit.hpp"

namespace {

constexpr zfit::ModelParameters true_parameters {
    .sigma_g = 4.94e-18,
    .sigma_s = 1.60e-17,
    .sigma_t = 1.50e-17,
    .tau10 = 1.00e-12,
    .tau21 = 1.00e-12,
    .tau13 = 1.00e-12,
    .tau30 = 1.29e-7,
    .tau43 = 1.00e-12,
};

constexpr zfit::SampleParameters sample {
    .thickness_m = 1.0e-3,
    .concentration_mol_per_l = 1.75e18 * 1000.0 / 6.02214076e23,
};

constexpr zfit::LaserParameters laser {
    .wavelength_m = 532.0e-9,
    .waist_m = 14.5e-6,
    .pulse_energy_j = 50.0e-9,
    .pulse_fwhm_s = 9.0e-9,
    .m2 = 1.2,
};

constexpr zfit::SolverOptions options {
    .nz = 11,
    .waveform_knots = 129,
    .pulse_start_fwhm = -0.5,
    .pulse_end_fwhm = 0.5,
};

std::vector<double> scan_positions()
{
    std::vector<double> positions;
    for (int i = -10; i <= 10; ++i) {
        positions.push_back(static_cast<double>(i) * 1.0e-3);
    }
    return positions;
}

std::vector<double> simulate(const zfit::ModelParameters& parameters,
    const std::vector<double>& positions, zfit::SolverOptions solver_options = options)
{
    zfit::MolModel model { parameters, sample, zfit::Laser { laser }, solver_options, positions };
    zfit::CvodesIntegrator integrator { model };
    std::vector<double> transmission;
    transmission.reserve(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i) {
        transmission.push_back(integrator.simulate(i));
    }
    return transmission;
}

double residual_norm(const std::vector<double>& predicted, const std::vector<double>& observed)
{
    double sum_squares { };
    for (std::size_t i = 0; i < observed.size(); ++i) {
        const double difference = predicted[i] - observed[i];
        sum_squares += difference * difference;
    }
    return std::sqrt(sum_squares);
}

zfit::FitRequest make_request()
{
    zfit::FitRequest request;
    request.initial_parameters = true_parameters;
    request.sample = sample;
    request.laser = laser;
    request.solver_options = options;
    request.scan_positions_m = scan_positions();
    request.observed_transmission = simulate(true_parameters, request.scan_positions_m);
    request.fitted_parameters = { zfit::FitParameter::Sigma_G };
    return request;
}

void set_diagonal_covariance(zfit::FitRequest& request, double sigma)
{
    const std::size_t m = request.scan_positions_m.size();
    request.measurement_covariance_row_major.assign(m * m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        request.measurement_covariance_row_major[i * m + i] = sigma * sigma;
    }
}

} // namespace

TEST_CASE("fit recovers two selected parameters from a synthetic curve")
{
    auto request = make_request();
    request.fitted_parameters = { zfit::FitParameter::Sigma_G, zfit::FitParameter::Sigma_T };
    request.initial_parameters.sigma_g *= 1.3;
    request.initial_parameters.sigma_t *= 0.7;

    const double initial_norm
        = residual_norm(simulate(request.initial_parameters, request.scan_positions_m),
            request.observed_transmission);
    REQUIRE(initial_norm > 0.0);

    const auto result = zfit::fit(request);
    INFO("status=" << result.status << ", evaluations=" << result.function_evaluations
                   << ", residual=" << result.residual_norm);

    REQUIRE(result.status >= 1);
    REQUIRE(result.status <= 4);
    REQUIRE_FALSE(result.chi_squared.has_value());
    REQUIRE_FALSE(result.uncertainty.has_value());
    REQUIRE(result.function_evaluations > 0);
    REQUIRE(result.residual_norm < initial_norm);
    REQUIRE(result.parameters.sigma_g == Catch::Approx(true_parameters.sigma_g).epsilon(0.05));
    REQUIRE(result.parameters.sigma_t == Catch::Approx(true_parameters.sigma_t).epsilon(0.05));
    REQUIRE(result.parameters.sigma_s == true_parameters.sigma_s);
    REQUIRE(result.parameters.tau13 == true_parameters.tau13);
    REQUIRE(result.predicted_transmission.size() == request.observed_transmission.size());
    REQUIRE(residual_norm(result.predicted_transmission, request.observed_transmission)
        == Catch::Approx(result.residual_norm));

    auto verification_options = options;
    verification_options.waveform_knots = 513;
    const auto verified
        = simulate(result.parameters, request.scan_positions_m, verification_options);
    const auto reference
        = simulate(true_parameters, request.scan_positions_m, verification_options);
    REQUIRE(residual_norm(verified, reference) < 1.0e-5);
}

TEST_CASE("fit accepts a single dynamically selected parameter")
{
    auto request = make_request();
    request.initial_parameters.sigma_g *= 1.3;

    const auto result = zfit::fit(request);

    REQUIRE(result.status >= 1);
    REQUIRE(result.status <= 4);
    REQUIRE(result.parameters.sigma_g == Catch::Approx(true_parameters.sigma_g).epsilon(0.05));
    REQUIRE(result.parameters.sigma_t == true_parameters.sigma_t);
}

TEST_CASE("fit uses diagonal measurement covariance as inverse-variance weights")
{
    auto request = make_request();
    request.initial_parameters.sigma_g *= 1.2;
    for (std::size_t i = 0; i < request.observed_transmission.size(); ++i) {
        request.observed_transmission[i] += 2.0e-4 * std::sin(0.8 * static_cast<double>(i));
    }

    constexpr double sigma = 5.0e-3;
    const std::size_t m = request.scan_positions_m.size();
    request.measurement_covariance_row_major.assign(m * m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        request.measurement_covariance_row_major[i * m + i] = sigma * sigma;
    }

    const auto result = zfit::fit(request);
    REQUIRE(result.status >= 1);
    REQUIRE(result.status <= 4);
    REQUIRE(result.chi_squared.has_value());
    REQUIRE(*result.chi_squared ==
            Catch::Approx(result.residual_norm * result.residual_norm / (sigma * sigma)));
    REQUIRE(result.parameters.sigma_g ==
            Catch::Approx(true_parameters.sigma_g).epsilon(0.05));
    REQUIRE(result.uncertainty.has_value());
    REQUIRE(result.uncertainty->parameter_order ==
            std::vector<zfit::FitParameter>{zfit::FitParameter::Sigma_G});
    REQUIRE(result.uncertainty->jacobian_rank == 1);
    REQUIRE(result.uncertainty->covariance_row_major.size() == 1);
    REQUIRE(result.uncertainty->standard_errors.size() == 1);
    REQUIRE(result.uncertainty->standard_errors[0] > 0.0);
    REQUIRE(result.uncertainty->standard_errors[0] ==
            Catch::Approx(std::sqrt(result.uncertainty->covariance_row_major[0])));

    // Independently check the local SE with a central log-parameter derivative.
    constexpr double step = 1.0e-3;
    auto plus = result.parameters;
    auto minus = result.parameters;
    plus.sigma_g *= std::exp(step);
    minus.sigma_g *= std::exp(-step);
    const auto plus_curve = simulate(plus, request.scan_positions_m);
    const auto minus_curve = simulate(minus, request.scan_positions_m);
    double whitened_jacobian_norm_squared{};
    for (std::size_t i = 0; i < m; ++i) {
        const double derivative = (plus_curve[i] - minus_curve[i]) / (2.0 * step * sigma);
        whitened_jacobian_norm_squared += derivative * derivative;
    }
    const double expected_se = result.parameters.sigma_g /
                               std::sqrt(whitened_jacobian_norm_squared);
    REQUIRE(result.uncertainty->standard_errors[0] ==
            Catch::Approx(expected_se).epsilon(0.1));
}

TEST_CASE("fit reports symmetric physical covariance in selected parameter order")
{
    auto request = make_request();
    request.fitted_parameters = {zfit::FitParameter::Sigma_T,
                                 zfit::FitParameter::Sigma_G};
    request.initial_parameters.sigma_t *= 0.7;
    request.initial_parameters.sigma_g *= 1.3;
    set_diagonal_covariance(request, 5.0e-3);

    const auto result = zfit::fit(request);
    REQUIRE(result.status >= 1);
    REQUIRE(result.status <= 4);
    REQUIRE(result.uncertainty.has_value());
    const auto& uncertainty = *result.uncertainty;
    REQUIRE(uncertainty.parameter_order == request.fitted_parameters);
    REQUIRE(uncertainty.jacobian_rank == 2);
    REQUIRE(uncertainty.covariance_row_major.size() == 4);
    REQUIRE(uncertainty.standard_errors.size() == 2);
    REQUIRE(uncertainty.covariance_row_major[1] == uncertainty.covariance_row_major[2]);
    for (std::size_t i = 0; i < 2; ++i) {
        REQUIRE(std::isfinite(uncertainty.standard_errors[i]));
        REQUIRE(uncertainty.standard_errors[i] > 0.0);
        REQUIRE(uncertainty.standard_errors[i] ==
                Catch::Approx(std::sqrt(uncertainty.covariance_row_major[i * 2 + i])));
    }
}

TEST_CASE("fit does not report SEs for a rank-deficient weighted fit")
{
    auto request = make_request();
    request.initial_parameters.sigma_g = 0.0;
    request.initial_parameters.sigma_s = 0.0;
    request.initial_parameters.sigma_t = 0.0;
    request.observed_transmission = simulate(request.initial_parameters,
                                             request.scan_positions_m);
    request.fitted_parameters = {zfit::FitParameter::Tau10,
                                 zfit::FitParameter::Tau30};
    set_diagonal_covariance(request, 5.0e-3);

    const auto result = zfit::fit(request);
    REQUIRE(result.status >= 1);
    REQUIRE(result.status <= 4);
    REQUIRE(result.uncertainty.has_value());
    REQUIRE(result.uncertainty->parameter_order == request.fitted_parameters);
    REQUIRE(result.uncertainty->jacobian_rank < request.fitted_parameters.size());
    REQUIRE(result.uncertainty->covariance_row_major.empty());
    REQUIRE(result.uncertainty->standard_errors.empty());
}

TEST_CASE("fit does not compute uncertainty after its evaluation limit")
{
    auto request = make_request();
    request.initial_parameters.sigma_g *= 1.3;
    request.max_evaluations = 1;
    set_diagonal_covariance(request, 5.0e-3);

    const auto result = zfit::fit(request);
    REQUIRE(result.status == 5);
    REQUIRE(result.chi_squared.has_value());
    REQUIRE_FALSE(result.uncertainty.has_value());
}

TEST_CASE("fit reports chi-squared for correlated measurement covariance")
{
    auto request = make_request();
    request.initial_parameters.sigma_g *= 1.2;
    for (std::size_t i = 0; i < request.observed_transmission.size(); ++i) {
        request.observed_transmission[i] += 2.0e-4 * std::sin(0.8 * static_cast<double>(i));
    }

    constexpr double sigma = 5.0e-3;
    constexpr double rho = 0.35;
    const std::size_t m = request.scan_positions_m.size();
    request.measurement_covariance_row_major.resize(m * m);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            request.measurement_covariance_row_major[i * m + j] =
                sigma * sigma * (i == j ? 1.0 : rho);
        }
    }

    const auto result = zfit::fit(request);
    REQUIRE(result.status >= 1);
    REQUIRE(result.status <= 4);
    REQUIRE(result.chi_squared.has_value());

    // Inverse of Sigma = sigma² [(1-rho) I + rho 11ᵀ].
    double sum_residuals{};
    double sum_squares{};
    for (std::size_t i = 0; i < m; ++i) {
        const double residual =
            result.predicted_transmission[i] - request.observed_transmission[i];
        sum_residuals += residual;
        sum_squares += residual * residual;
    }
    const double expected_chi_squared =
        (sum_squares - rho * sum_residuals * sum_residuals /
                           (1.0 - rho + static_cast<double>(m) * rho)) /
        (sigma * sigma * (1.0 - rho));
    REQUIRE(*result.chi_squared == Catch::Approx(expected_chi_squared).epsilon(1e-8));
    REQUIRE(result.residual_norm == Catch::Approx(std::sqrt(sum_squares)));
}

TEST_CASE("fit downweights a noisy outlier in the optimized parameters")
{
    auto request = make_request();
    request.scan_positions_m = {request.scan_positions_m.front(), request.scan_positions_m[10]};
    request.observed_transmission = simulate(true_parameters, request.scan_positions_m);
    request.observed_transmission[0] += 0.05;
    request.initial_parameters.sigma_g *= 1.2;

    const auto unweighted = zfit::fit(request);
    request.measurement_covariance_row_major = {
        0.1 * 0.1, 0.0,
        0.0, 0.001 * 0.001,
    };
    const auto weighted = zfit::fit(request);

    INFO("unweighted sigma_g=" << unweighted.parameters.sigma_g
                               << ", weighted sigma_g=" << weighted.parameters.sigma_g);
    REQUIRE_FALSE(unweighted.chi_squared.has_value());
    REQUIRE(weighted.chi_squared.has_value());
    REQUIRE(std::abs(weighted.parameters.sigma_g - true_parameters.sigma_g) <
            std::abs(unweighted.parameters.sigma_g - true_parameters.sigma_g));
}

TEST_CASE("fit rejects malformed or singular measurement covariance")
{
    auto request = make_request();
    const std::size_t m = request.scan_positions_m.size();
    request.measurement_covariance_row_major.assign(m * m - 1, 0.0);
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);

    request.measurement_covariance_row_major.assign(m * m, 1.0e-6);
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);

    request.measurement_covariance_row_major.assign(m * m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        request.measurement_covariance_row_major[i * m + i] = 1.0e-6;
    }
    request.measurement_covariance_row_major[0] =
        std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);
}

TEST_CASE("fit rejects duplicate selected parameters")
{
    auto request = make_request();
    request.fitted_parameters = { zfit::FitParameter::Sigma_G, zfit::FitParameter::Sigma_G };
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);
}

TEST_CASE("fit rejects nonpositive starting values for selected parameters")
{
    auto request = make_request();
    request.initial_parameters.sigma_g = 0.0;
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);
}

TEST_CASE("fit rejects mismatched scan and transmission sizes")
{
    auto request = make_request();
    request.observed_transmission.pop_back();
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);
}
