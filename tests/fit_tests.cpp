#include <algorithm>
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
    const std::vector<double>& positions, zfit::SolverOptions solver_options = options,
    zfit::LaserParameters curve_laser = laser)
{
    zfit::MolModel model { parameters, sample, zfit::Laser { curve_laser }, solver_options,
        positions };
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

double max_abs_residual(const std::vector<double>& predicted, const std::vector<double>& observed)
{
    double maximum { };
    for (std::size_t i = 0; i < observed.size(); ++i) {
        maximum = std::max(maximum, std::abs(predicted[i] - observed[i]));
    }
    return maximum;
}

zfit::LaserParameters laser_at(double pulse_energy_j)
{
    auto adjusted = laser;
    adjusted.pulse_energy_j = pulse_energy_j;
    return adjusted;
}

zfit::EnergyCurve make_curve(const std::vector<double>& positions,
    const std::vector<double>& transmission, zfit::LaserParameters curve_laser = laser)
{
    zfit::EnergyCurve curve;
    curve.laser = curve_laser;
    curve.scan_positions_m = positions;
    curve.observed_transmission = transmission;
    return curve;
}

zfit::FitRequest make_request()
{
    zfit::FitRequest request;
    request.initial_parameters = true_parameters;
    request.sample = sample;
    request.solver_options = options;
    const auto positions = scan_positions();
    request.curves.push_back(make_curve(positions, simulate(true_parameters, positions)));
    request.fitted_parameters = { zfit::FitParameter::Sigma_G };
    return request;
}

zfit::FitRequest make_optical_request()
{
    zfit::FitRequest request;
    request.initial_parameters = true_parameters;
    request.sample = sample;
    request.solver_options = options;

    std::vector<double> positions;
    for (int i = -10; i <= 10; ++i) {
        positions.push_back(static_cast<double>(i) * 0.5e-3);
    }
    for (double energy_j : { 50.0e-9, 150.0e-9 }) {
        const auto curve_laser = laser_at(energy_j);
        request.curves.push_back(make_curve(positions,
            simulate(true_parameters, positions, options, curve_laser), curve_laser));
    }
    return request;
}

void set_diagonal_covariance(zfit::FitRequest& request, double sigma)
{
    for (auto& curve : request.curves) {
        const std::size_t m = curve.scan_positions_m.size();
        curve.measurement_covariance_row_major.assign(m * m, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            curve.measurement_covariance_row_major[i * m + i] = sigma * sigma;
        }
    }
}

} // namespace

TEST_CASE("fit recovers two selected parameters from a synthetic curve")
{
    auto request = make_request();
    request.fitted_parameters = { zfit::FitParameter::Sigma_G, zfit::FitParameter::Sigma_T };
    request.initial_parameters.sigma_g *= 1.3;
    request.initial_parameters.sigma_t *= 0.7;

    const auto& positions = request.curves.front().scan_positions_m;
    const double initial_norm = residual_norm(simulate(request.initial_parameters, positions),
        request.curves.front().observed_transmission);
    REQUIRE(initial_norm > 0.0);

    const auto result = zfit::fit(request);
    INFO("status=" << result.status << ", evaluations=" << result.function_evaluations
                   << ", residual=" << result.residual_norm);

    REQUIRE(result.status >= 1);
    REQUIRE(result.status <= 4);
    REQUIRE(result.converged);
    REQUIRE_FALSE(result.chi_squared.has_value());
    REQUIRE_FALSE(result.uncertainty.has_value());
    REQUIRE(result.curves.size() == 1);
    REQUIRE(result.function_evaluations > 0);
    REQUIRE(result.residual_norm < initial_norm);
    REQUIRE(result.parameters.sigma_g == Catch::Approx(true_parameters.sigma_g).epsilon(0.05));
    REQUIRE(result.parameters.sigma_t == Catch::Approx(true_parameters.sigma_t).epsilon(0.05));
    REQUIRE(result.parameters.sigma_s == true_parameters.sigma_s);
    REQUIRE(result.parameters.tau13 == true_parameters.tau13);
    REQUIRE(result.curves.front().predicted_transmission.size()
        == request.curves.front().observed_transmission.size());
    REQUIRE(residual_norm(result.curves.front().predicted_transmission,
                request.curves.front().observed_transmission)
        == Catch::Approx(result.residual_norm));

    auto verification_options = options;
    verification_options.waveform_knots = 513;
    const auto verified = simulate(result.parameters, positions, verification_options);
    const auto reference = simulate(true_parameters, positions, verification_options);
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
    auto& transmission = request.curves.front().observed_transmission;
    for (std::size_t i = 0; i < transmission.size(); ++i) {
        transmission[i] += 2.0e-4 * std::sin(0.8 * static_cast<double>(i));
    }

    constexpr double sigma = 5.0e-3;
    set_diagonal_covariance(request, sigma);
    const std::size_t m = request.curves.front().scan_positions_m.size();

    const auto result = zfit::fit(request);
    REQUIRE(result.status >= 1);
    REQUIRE(result.status <= 4);
    REQUIRE(result.converged);
    REQUIRE(result.chi_squared.has_value());
    REQUIRE(*result.chi_squared
        == Catch::Approx(result.residual_norm * result.residual_norm / (sigma * sigma)));
    REQUIRE(result.parameters.sigma_g == Catch::Approx(true_parameters.sigma_g).epsilon(0.05));
    REQUIRE(result.uncertainty.has_value());
    REQUIRE(result.uncertainty->parameter_order
        == std::vector<zfit::FitParameter> { zfit::FitParameter::Sigma_G });
    REQUIRE(result.uncertainty->jacobian_rank == 1);
    REQUIRE(result.uncertainty->covariance_row_major.size() == 1);
    REQUIRE(result.uncertainty->standard_errors.size() == 1);
    REQUIRE(result.uncertainty->standard_errors[0] > 0.0);
    REQUIRE(result.uncertainty->standard_errors[0]
        == Catch::Approx(std::sqrt(result.uncertainty->covariance_row_major[0])));

    // Independently check the local SE with a central log-parameter derivative.
    constexpr double step = 1.0e-3;
    auto plus = result.parameters;
    auto minus = result.parameters;
    plus.sigma_g *= std::exp(step);
    minus.sigma_g *= std::exp(-step);
    const auto plus_curve = simulate(plus, request.curves.front().scan_positions_m);
    const auto minus_curve = simulate(minus, request.curves.front().scan_positions_m);
    double whitened_jacobian_norm_squared { };
    for (std::size_t i = 0; i < m; ++i) {
        const double derivative = (plus_curve[i] - minus_curve[i]) / (2.0 * step * sigma);
        whitened_jacobian_norm_squared += derivative * derivative;
    }
    const double expected_se
        = result.parameters.sigma_g / std::sqrt(whitened_jacobian_norm_squared);
    REQUIRE(result.uncertainty->standard_errors[0] == Catch::Approx(expected_se).epsilon(0.1));
}

TEST_CASE("fit reports symmetric physical covariance in selected parameter order")
{
    auto request = make_request();
    request.fitted_parameters = { zfit::FitParameter::Sigma_T, zfit::FitParameter::Sigma_G };
    request.initial_parameters.sigma_t *= 0.7;
    request.initial_parameters.sigma_g *= 1.3;
    set_diagonal_covariance(request, 5.0e-3);

    const auto result = zfit::fit(request);
    REQUIRE(result.status >= 1);
    REQUIRE(result.status <= 4);
    REQUIRE(result.converged);
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
        REQUIRE(uncertainty.standard_errors[i]
            == Catch::Approx(std::sqrt(uncertainty.covariance_row_major[i * 2 + i])));
    }
}

TEST_CASE("fit does not report SEs for a rank-deficient weighted fit")
{
    auto request = make_request();
    request.initial_parameters.sigma_g = 0.0;
    request.initial_parameters.sigma_s = 0.0;
    request.initial_parameters.sigma_t = 0.0;
    request.curves.front().observed_transmission
        = simulate(request.initial_parameters, request.curves.front().scan_positions_m);
    request.fitted_parameters = { zfit::FitParameter::Tau10, zfit::FitParameter::Tau30 };
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
    REQUIRE_FALSE(result.converged);
    REQUIRE(result.initial_residual_norm > 0.0);
    REQUIRE(result.chi_squared.has_value());
    REQUIRE_FALSE(result.uncertainty.has_value());
}

TEST_CASE("fit reports raw diagnostics and a converged flag")
{
    auto request = make_request();
    request.initial_parameters.sigma_g *= 1.3;
    const auto& positions = request.curves.front().scan_positions_m;
    const double initial_reference = residual_norm(simulate(request.initial_parameters, positions),
        request.curves.front().observed_transmission);

    const auto result = zfit::fit(request);
    REQUIRE(result.converged);
    REQUIRE(result.initial_residual_norm == Catch::Approx(initial_reference).epsilon(1e-9));
    REQUIRE(result.initial_residual_norm > result.residual_norm);
    REQUIRE(result.max_abs_residual
        == Catch::Approx(max_abs_residual(result.curves.front().predicted_transmission,
            request.curves.front().observed_transmission)));
    REQUIRE(result.max_abs_residual <= result.residual_norm);

    // The diagnostic stays raw when a measurement covariance is present.
    set_diagonal_covariance(request, 5.0e-3);
    const auto weighted = zfit::fit(request);
    REQUIRE(weighted.converged);
    REQUIRE(weighted.initial_residual_norm == Catch::Approx(initial_reference).epsilon(1e-9));
    REQUIRE(weighted.initial_residual_norm > weighted.residual_norm);
    REQUIRE(weighted.max_abs_residual <= weighted.residual_norm);
}

TEST_CASE("fit uncertainties are stable against the finite-difference step")
{
    auto request = make_request();
    request.initial_parameters.sigma_g *= 1.2;
    set_diagonal_covariance(request, 5.0e-3);

    std::vector<double> standard_errors;
    for (double epsfcn : { 1e-3, 1e-6, 1e-9 }) {
        // epsfcn sets both the optimizer's forward-difference step and the
        // final Jacobian's step (root of epsfcn in log coordinates), so the
        // fits take slightly different paths.
        request.epsfcn = epsfcn;
        const auto result = zfit::fit(request);
        REQUIRE(result.converged);
        REQUIRE(result.uncertainty.has_value());
        standard_errors.push_back(result.uncertainty->standard_errors[0]);
    }
    for (double se : standard_errors) {
        REQUIRE(se == Catch::Approx(standard_errors[1]).epsilon(0.10));
    }
}

// Offline validation note: epsfcn below ~1e-9 is invalid, not just noisy.
// lmdif's internal step is sqrt(epsfcn)*|x| in log coordinates; at 1e-12
// (a 1e-6 relative parameter change) the step falls below CVODES solver
// noise, the fit itself converges to wrong parameters, and any SE it
// reports is meaningless. Keep epsfcn at or above the default 1e-6.

TEST_CASE("fit uncertainties are stable at 513 waveform knots")
{
    auto request = make_request();
    request.initial_parameters.sigma_g *= 1.2;
    set_diagonal_covariance(request, 5.0e-3);

    auto fine_request = request;
    fine_request.solver_options.waveform_knots = 513;

    const auto coarse = zfit::fit(request);
    const auto fine = zfit::fit(fine_request);

    REQUIRE(coarse.converged);
    REQUIRE(fine.converged);
    REQUIRE(fine.uncertainty.has_value());
    REQUIRE(coarse.uncertainty->standard_errors[0]
        == Catch::Approx(fine.uncertainty->standard_errors[0]).epsilon(0.10));
}

TEST_CASE("fit reports chi-squared for correlated measurement covariance")
{
    auto request = make_request();
    request.initial_parameters.sigma_g *= 1.2;
    auto& transmission = request.curves.front().observed_transmission;
    for (std::size_t i = 0; i < transmission.size(); ++i) {
        transmission[i] += 2.0e-4 * std::sin(0.8 * static_cast<double>(i));
    }

    constexpr double sigma = 5.0e-3;
    constexpr double rho = 0.35;
    auto& curve = request.curves.front();
    const std::size_t m = curve.scan_positions_m.size();
    curve.measurement_covariance_row_major.resize(m * m);
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            curve.measurement_covariance_row_major[i * m + j]
                = sigma * sigma * (i == j ? 1.0 : rho);
        }
    }

    const auto result = zfit::fit(request);
    REQUIRE(result.status >= 1);
    REQUIRE(result.status <= 4);
    REQUIRE(result.chi_squared.has_value());

    // Inverse of Sigma = sigma² [(1-rho) I + rho 11ᵀ].
    double sum_residuals { };
    double sum_squares { };
    for (std::size_t i = 0; i < m; ++i) {
        const double residual
            = result.curves.front().predicted_transmission[i] - curve.observed_transmission[i];
        sum_residuals += residual;
        sum_squares += residual * residual;
    }
    const double expected_chi_squared
        = (sum_squares
              - rho * sum_residuals * sum_residuals / (1.0 - rho + static_cast<double>(m) * rho))
        / (sigma * sigma * (1.0 - rho));
    REQUIRE(*result.chi_squared == Catch::Approx(expected_chi_squared).epsilon(1e-8));
    REQUIRE(result.residual_norm == Catch::Approx(std::sqrt(sum_squares)));
}

TEST_CASE("fit downweights a noisy outlier in the optimized parameters")
{
    auto request = make_request();
    auto& curve = request.curves.front();
    curve.scan_positions_m = { curve.scan_positions_m.front(), curve.scan_positions_m[10] };
    curve.observed_transmission = simulate(true_parameters, curve.scan_positions_m);
    curve.observed_transmission[0] += 0.05;
    request.initial_parameters.sigma_g *= 1.2;

    const auto unweighted = zfit::fit(request);
    curve.measurement_covariance_row_major = {
        0.1 * 0.1,
        0.0,
        0.0,
        0.001 * 0.001,
    };
    const auto weighted = zfit::fit(request);

    INFO("unweighted sigma_g=" << unweighted.parameters.sigma_g
                               << ", weighted sigma_g=" << weighted.parameters.sigma_g);
    REQUIRE_FALSE(unweighted.chi_squared.has_value());
    REQUIRE(weighted.chi_squared.has_value());
    REQUIRE(std::abs(weighted.parameters.sigma_g - true_parameters.sigma_g)
        < std::abs(unweighted.parameters.sigma_g - true_parameters.sigma_g));
}

TEST_CASE("fit rejects malformed or singular measurement covariance")
{
    auto request = make_request();
    auto& curve = request.curves.front();
    const std::size_t m = curve.scan_positions_m.size();
    curve.measurement_covariance_row_major.assign(m * m - 1, 0.0);
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);

    curve.measurement_covariance_row_major.assign(m * m, 1.0e-6);
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);

    curve.measurement_covariance_row_major.assign(m * m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        curve.measurement_covariance_row_major[i * m + i] = 1.0e-6;
    }
    curve.measurement_covariance_row_major[0] = std::numeric_limits<double>::quiet_NaN();
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
    request.curves.front().observed_transmission.pop_back();
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);
}

TEST_CASE("fit rejects an empty curve list")
{
    auto request = make_request();
    request.curves.clear();
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);
}

TEST_CASE("joint fit requires all or no curve covariances")
{
    auto request = make_request();
    set_diagonal_covariance(request, 5.0e-3);
    const auto positions = scan_positions();
    request.curves.push_back(make_curve(positions, simulate(true_parameters, positions)));
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);
}

TEST_CASE("joint fit recovers shared parameters across pulse energies")
{
    auto request = make_request();
    request.initial_parameters.sigma_g *= 1.3;
    request.initial_parameters.tau13 *= 0.7;
    request.fitted_parameters = { zfit::FitParameter::Sigma_G, zfit::FitParameter::Tau13 };
    request.curves.clear();

    for (double energy_j : { 30.0e-9, 50.0e-9, 80.0e-9 }) {
        const auto positions = scan_positions();
        const auto curve_laser = laser_at(energy_j);
        request.curves.push_back(make_curve(
            positions, simulate(true_parameters, positions, options, curve_laser), curve_laser));
    }

    const auto result = zfit::fit(request);
    INFO("status=" << result.status << ", evaluations=" << result.function_evaluations
                   << ", residual=" << result.residual_norm);

    REQUIRE(result.status >= 1);
    REQUIRE(result.status <= 4);
    REQUIRE(result.converged);
    REQUIRE_FALSE(result.chi_squared.has_value());
    REQUIRE_FALSE(result.uncertainty.has_value());
    REQUIRE(result.curves.size() == request.curves.size());

    for (std::size_t k = 0; k < result.curves.size(); ++k) {
        REQUIRE(result.curves[k].predicted_transmission.size()
            == request.curves[k].scan_positions_m.size());
        // Loose bound: the optimize stops on xtol, not machine precision.
        REQUIRE(residual_norm(result.curves[k].predicted_transmission,
                    request.curves[k].observed_transmission)
            < 1.0e-5);
    }

    REQUIRE(result.parameters.sigma_g == Catch::Approx(true_parameters.sigma_g).epsilon(0.05));
    REQUIRE(result.parameters.tau13 == Catch::Approx(true_parameters.tau13).epsilon(0.05));
    REQUIRE(result.parameters.sigma_t == true_parameters.sigma_t);
}

TEST_CASE("joint fit narrows standard errors relative to a single curve")
{
    // The same truth and covariance: a joint fit of three pulse energies
    // contains strictly more information than one curve, so its standard
    // errors must be smaller.
    auto single = make_request();
    single.initial_parameters.sigma_g *= 1.2;
    set_diagonal_covariance(single, 5.0e-3);
    auto joint = single;

    joint.curves.clear();
    for (double energy_j : { 30.0e-9, 50.0e-9, 80.0e-9 }) {
        const auto positions = scan_positions();
        const auto curve_laser = laser_at(energy_j);
        joint.curves.push_back(make_curve(
            positions, simulate(true_parameters, positions, options, curve_laser), curve_laser));
    }
    set_diagonal_covariance(joint, 5.0e-3);

    const auto single_result = zfit::fit(single);
    const auto joint_result = zfit::fit(joint);

    INFO("single se=" << single_result.uncertainty->standard_errors[0]
                      << ", joint se=" << joint_result.uncertainty->standard_errors[0]);
    REQUIRE(single_result.converged);
    REQUIRE(joint_result.converged);
    REQUIRE(single_result.uncertainty.has_value());
    REQUIRE(joint_result.uncertainty.has_value());
    REQUIRE(
        single_result.uncertainty->parameter_order == joint_result.uncertainty->parameter_order);
    REQUIRE(single_result.uncertainty->standard_errors[0]
        > joint_result.uncertainty->standard_errors[0]);
}

TEST_CASE("joint fit honors distinct per-curve covariance weights")
{
    const auto positions = scan_positions();
    // At this pulse energy sigma_t is strongly identified, so a biased
    // curve visibly drags it when fit alone.
    const auto curve_laser = laser_at(80.0e-9);

    zfit::FitRequest request;
    request.initial_parameters = true_parameters;
    request.initial_parameters.sigma_t *= 1.3;
    request.sample = sample;
    request.solver_options = options;
    request.fitted_parameters = { zfit::FitParameter::Sigma_T };

    auto set_curve_diagonal_covariance = [](zfit::EnergyCurve& curve, double sigma) {
        const std::size_t m = curve.scan_positions_m.size();
        curve.measurement_covariance_row_major.assign(m * m, 0.0);
        for (std::size_t i = 0; i < m; ++i) {
            curve.measurement_covariance_row_major[i * m + i] = sigma * sigma;
        }
    };

    // A loose curve biased by twenty of its standard deviations: on its own
    // it drags sigma_t far from the truth.
    auto biased_transmission = simulate(true_parameters, positions, options, curve_laser);
    for (auto& value : biased_transmission) {
        value += 2.0e-2;
    }

    zfit::EnergyCurve biased_curve = make_curve(positions, biased_transmission, curve_laser);
    zfit::EnergyCurve truth_curve = make_curve(
        positions, simulate(true_parameters, positions, options, curve_laser), curve_laser);
    set_curve_diagonal_covariance(biased_curve, 1.0e-3);
    set_curve_diagonal_covariance(truth_curve, 1.0e-6);

    request.curves.push_back(biased_curve);
    request.curves.push_back(truth_curve);
    const auto joint_result = zfit::fit(request);

    auto loose_only = request;
    loose_only.curves.pop_back();
    const auto loose_result = zfit::fit(loose_only);

    const double joint_error = std::abs(joint_result.parameters.sigma_t - true_parameters.sigma_t);
    const double loose_error = std::abs(loose_result.parameters.sigma_t - true_parameters.sigma_t);
    INFO("joint_error=" << joint_error << ", loose_error=" << loose_error);
    REQUIRE(joint_result.converged);
    REQUIRE(loose_result.converged);
    // The tight curve dominates the joint pull: sigma_t stays within 1%
    // of the truth, while the biased curve alone drags it far away.
    REQUIRE(joint_error < 0.01 * true_parameters.sigma_t);
    REQUIRE(loose_error > 0.2 * true_parameters.sigma_t);
    REQUIRE(joint_error < 0.01 * loose_error);
}

TEST_CASE("joint fit recovers each shared optical parameter independently")
{
    for (auto id : { zfit::FitParameter::Waist, zfit::FitParameter::M2 }) {
        auto request = make_optical_request();
        request.fitted_parameters = { id };
        for (auto& curve : request.curves) {
            if (id == zfit::FitParameter::Waist) {
                curve.laser.waist_m *= 1.15;
            } else {
                curve.laser.m2 = 1.4;
            }
        }

        const auto result = zfit::fit(request);
        INFO("parameter=" << static_cast<int>(id) << ", status=" << result.status
                          << ", residual=" << result.residual_norm);
        REQUIRE(result.converged);
        REQUIRE(result.residual_norm < 1.0e-5);
        REQUIRE(result.fitted_waist_m.has_value() == (id == zfit::FitParameter::Waist));
        REQUIRE(result.fitted_m2.has_value() == (id == zfit::FitParameter::M2));
        for (std::size_t k = 0; k < result.curves.size(); ++k) {
            REQUIRE(result.curves[k].laser.pulse_energy_j == request.curves[k].laser.pulse_energy_j);
            REQUIRE(result.curves[k].laser.waist_m
                == Catch::Approx(laser.waist_m).epsilon(0.05));
            REQUIRE(result.curves[k].laser.m2 == Catch::Approx(laser.m2).epsilon(0.05));
        }
    }
}

TEST_CASE("joint fit of kinetics and optics reports ordered covariance without double counting")
{
    auto request = make_optical_request();
    request.fitted_parameters
        = { zfit::FitParameter::M2, zfit::FitParameter::Sigma_T, zfit::FitParameter::Waist };
    request.initial_parameters.sigma_t *= 1.02;
    for (auto& curve : request.curves) {
        curve.laser.waist_m *= 1.05;
        curve.laser.m2 = 1.25;
    }
    set_diagonal_covariance(request, 5.0e-3);
    request.calibration.waist_std_m = 0.5e-6;
    request.calibration.m2_std = 0.04;
    request.calibration.energy_meter_std_j = 1.0e-9;

    const auto result = zfit::fit(request);
    INFO("status=" << result.status << ", residual=" << result.residual_norm);
    REQUIRE(result.converged);
    REQUIRE(result.fitted_waist_m.has_value());
    REQUIRE(result.fitted_m2.has_value());
    REQUIRE(*result.fitted_waist_m == Catch::Approx(laser.waist_m).epsilon(0.05));
    REQUIRE(*result.fitted_m2 == Catch::Approx(laser.m2).epsilon(0.05));
    REQUIRE(result.parameters.sigma_t == Catch::Approx(true_parameters.sigma_t).epsilon(0.05));
    REQUIRE(result.uncertainty.has_value());
    const auto& uncertainty = *result.uncertainty;
    REQUIRE(uncertainty.parameter_order == request.fitted_parameters);
    REQUIRE(uncertainty.jacobian_rank == request.fitted_parameters.size());
    REQUIRE(uncertainty.covariance_row_major.size() == 9);
    REQUIRE(uncertainty.total_covariance_row_major.size() == 9);
    REQUIRE(uncertainty.calibration_variance_contributions.size() == 3);
    for (std::size_t i = 0; i < 3; ++i) {
        REQUIRE(uncertainty.standard_errors[i] > 0.0);
        REQUIRE(uncertainty.total_standard_errors[i] >= uncertainty.standard_errors[i]);
        REQUIRE(uncertainty.calibration_variance_contributions[i][3] == 0.0);
        REQUIRE(uncertainty.calibration_variance_contributions[i][4] == 0.0);
        double total = uncertainty.covariance_row_major[i * 3 + i];
        for (double contribution : uncertainty.calibration_variance_contributions[i]) {
            total += contribution;
        }
        REQUIRE(uncertainty.total_covariance_row_major[i * 3 + i]
            == Catch::Approx(total).epsilon(1e-8));
        for (std::size_t j = 0; j < 3; ++j) {
            REQUIRE(uncertainty.covariance_row_major[i * 3 + j]
                == Catch::Approx(uncertainty.covariance_row_major[j * 3 + i]));
        }
    }
    for (const auto& curve : result.curves) {
        REQUIRE(curve.laser.waist_m == *result.fitted_waist_m);
        REQUIRE(curve.laser.m2 == *result.fitted_m2);
    }
}

TEST_CASE("optical fit validates initial values and preserves fixed ideal M2")
{
    auto request = make_optical_request();
    request.fitted_parameters = { zfit::FitParameter::M2 };
    request.curves.front().laser.m2 = 1.0;
    request.curves.back().laser.m2 = 1.0;
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);
    request.curves.front().laser.m2 = 1.4;
    request.curves.back().laser.m2 = 1.3;
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);

    request.fitted_parameters = { zfit::FitParameter::Waist };
    request.curves.front().laser.waist_m = 0.0;
    request.curves.back().laser.waist_m = 0.0;
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);

    request = make_optical_request();
    request.fitted_parameters = { zfit::FitParameter::Waist, zfit::FitParameter::Waist };
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);
    request.fitted_parameters = { zfit::FitParameter::Sigma_G };
    request.calibration.fixed_parameters = { { zfit::FitParameter::Waist, 1.0e-6 } };
    REQUIRE_THROWS_AS(zfit::fit(request), std::invalid_argument);

    auto ideal = make_request();
    ideal.curves.front().laser.m2 = 1.0;
    ideal.curves.front().observed_transmission
        = simulate(true_parameters, ideal.curves.front().scan_positions_m, options,
            ideal.curves.front().laser);
    ideal.initial_parameters.sigma_g *= 1.1;
    const auto fitted = zfit::fit(ideal);
    REQUIRE(fitted.converged);
    REQUIRE_FALSE(fitted.fitted_m2.has_value());
    REQUIRE(fitted.curves.front().laser.m2 == 1.0);
}

TEST_CASE("M2 SE uses the shifted-log derivative and stays stable at 513 knots")
{
    auto request = make_optical_request();
    request.fitted_parameters = { zfit::FitParameter::M2 };
    for (auto& curve : request.curves) {
        curve.laser.m2 = 1.4;
    }
    set_diagonal_covariance(request, 5.0e-3);

    const auto result = zfit::fit(request);
    REQUIRE(result.converged);
    REQUIRE(result.fitted_m2.has_value());
    REQUIRE(result.uncertainty.has_value());
    REQUIRE(result.uncertainty->jacobian_rank == 1);
    REQUIRE(result.uncertainty->standard_errors.size() == 1);

    // An independent central derivative in y = log(M2 - 1) checks that the
    // covariance is converted with dM2/dy = M2 - 1, not M2.
    constexpr double step = 1.0e-3;
    constexpr double noise = 5.0e-3;
    double norm_squared = 0.0;
    for (std::size_t k = 0; k < result.curves.size(); ++k) {
        const auto& curve = result.curves[k];
        auto plus = curve.laser;
        auto minus = curve.laser;
        plus.m2 = 1.0 + (*result.fitted_m2 - 1.0) * std::exp(step);
        minus.m2 = 1.0 + (*result.fitted_m2 - 1.0) * std::exp(-step);
        const auto& positions = request.curves[k].scan_positions_m;
        const auto plus_curve = simulate(result.parameters, positions, options, plus);
        const auto minus_curve = simulate(result.parameters, positions, options, minus);
        for (std::size_t i = 0; i < positions.size(); ++i) {
            const double derivative = (plus_curve[i] - minus_curve[i]) / (2.0 * step * noise);
            norm_squared += derivative * derivative;
        }
    }
    REQUIRE(result.uncertainty->standard_errors[0]
        == Catch::Approx((*result.fitted_m2 - 1.0) / std::sqrt(norm_squared)).epsilon(0.1));

    auto fine = request;
    fine.solver_options.waveform_knots = 513;
    const auto verified = zfit::fit(fine);
    REQUIRE(verified.converged);
    REQUIRE(verified.uncertainty.has_value());
    REQUIRE(verified.uncertainty->standard_errors[0]
        == Catch::Approx(result.uncertainty->standard_errors[0]).epsilon(0.2));
}

TEST_CASE("fitting only waist still propagates fixed M2 uncertainty")
{
    auto request = make_optical_request();
    request.fitted_parameters = { zfit::FitParameter::Waist };
    for (auto& curve : request.curves) {
        curve.laser.waist_m *= 1.15;
    }
    set_diagonal_covariance(request, 5.0e-3);
    request.calibration.waist_std_m = 0.5e-6;
    request.calibration.m2_std = 0.04;

    const auto result = zfit::fit(request);
    REQUIRE(result.converged);
    REQUIRE(result.uncertainty.has_value());
    REQUIRE(result.uncertainty->calibration_variance_contributions.size() == 1);
    REQUIRE(result.uncertainty->calibration_variance_contributions[0][3] == 0.0);
    REQUIRE(result.uncertainty->calibration_variance_contributions[0][4] > 0.0);
    REQUIRE(result.uncertainty->total_standard_errors[0]
        > result.uncertainty->standard_errors[0]);
}

TEST_CASE("M2 unidentifiable at focus reports rank without covariance")
{
    auto request = make_optical_request();
    for (auto& curve : request.curves) {
        // At z=0, the irradiance has no dependence on the Rayleigh range,
        // hence no dependence on M2; waist is still identifiable.
        curve.scan_positions_m = { 0.0, 0.0, 0.0 };
        curve.observed_transmission = simulate(true_parameters, curve.scan_positions_m,
            options, curve.laser);
        curve.laser.waist_m *= 1.15;
        curve.laser.m2 = 1.4;
    }
    request.fitted_parameters = { zfit::FitParameter::Waist, zfit::FitParameter::M2 };
    set_diagonal_covariance(request, 5.0e-3);

    const auto result = zfit::fit(request);
    REQUIRE(result.converged);
    REQUIRE(result.uncertainty.has_value());
    REQUIRE(result.uncertainty->jacobian_rank < 2);
    REQUIRE(result.uncertainty->covariance_row_major.empty());
    REQUIRE(result.uncertainty->standard_errors.empty());
}
