#include <numeric>
#include <random>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/cvodes_integrator.hpp"
#include "zfit/fit.hpp"
#include "zfit/model.hpp"

namespace {

constexpr zfit::ModelParameters truth_parameters {
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

constexpr zfit::LaserParameters base_laser {
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
    const std::vector<double>& positions, const zfit::LaserParameters& laser)
{
    zfit::MolModel model { parameters, sample, zfit::Laser { laser }, options, positions };
    zfit::CvodesIntegrator integrator { model };
    std::vector<double> transmission;
    transmission.reserve(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i) {
        transmission.push_back(integrator.simulate(i));
    }
    return transmission;
}

zfit::FitRequest make_two_curve_request(std::vector<double> pulse_energies_j)
{
    zfit::FitRequest request;
    request.initial_parameters = truth_parameters;
    request.sample = sample;
    request.solver_options = options;

    for (double energy : pulse_energies_j) {
        auto laser = base_laser;
        laser.pulse_energy_j = energy;
        const auto positions = scan_positions();
        zfit::EnergyCurve curve;
        curve.laser = laser;
        curve.scan_positions_m = positions;
        curve.observed_transmission = simulate(truth_parameters, positions, laser);
        curve.measurement_covariance_row_major.assign(positions.size() * positions.size(), 0.0);
        for (std::size_t i = 0; i < positions.size(); ++i) {
            curve.measurement_covariance_row_major[i * positions.size() + i] = 0.01 * 0.01;
        }
        request.curves.push_back(std::move(curve));
    }
    request.fitted_parameters = { zfit::FitParameter::Sigma_G };
    return request;
}

constexpr double waist = base_laser.waist_m;
constexpr double pulse_energies[2] = { 50.0e-9, 80.0e-9 };

} // namespace

TEST_CASE("fit leaves total covariance absent without calibration uncertainties")
{
    auto request = make_two_curve_request({ pulse_energies[0], pulse_energies[1] });
    const auto result = zfit::fit(request);

    REQUIRE(result.uncertainty.has_value());
    REQUIRE(result.uncertainty->total_covariance_row_major.empty());
    REQUIRE(result.uncertainty->total_standard_errors.empty());
    REQUIRE(result.uncertainty->calibration_variance_contributions.empty());
}

TEST_CASE("fit adds calibration propagation to conditional SEs")
{
    auto request = make_two_curve_request({ pulse_energies[0], pulse_energies[1] });
    request.calibration.pulse_energy_stds_j = { 2.0e-9, 2.0e-9 };
    request.calibration.entrance_chain_std = 0.006;
    request.calibration.waist_std_m = 0.725e-6;
    const auto result = zfit::fit(request);

    REQUIRE(result.uncertainty.has_value());
    const auto& uncertainty = *result.uncertainty;
    REQUIRE(uncertainty.total_standard_errors.size() == 1);
    REQUIRE(uncertainty.total_standard_errors[0] > uncertainty.standard_errors[0]);
    REQUIRE(uncertainty.total_covariance_row_major.size() == 1);

    // Contributions decompose the total variance exactly: conditional plus
    // the supplied groups 0 (energy measurements), 2 (entrance chain), and
    // 3 (waist).
    const double conditional = uncertainty.standard_errors[0] * uncertainty.standard_errors[0];
    const double energy_group = uncertainty.calibration_variance_contributions[0][0];
    const double entrance_group = uncertainty.calibration_variance_contributions[0][2];
    const double waist_group = uncertainty.calibration_variance_contributions[0][3];
    const double total
        = uncertainty.total_standard_errors[0] * uncertainty.total_standard_errors[0];
    REQUIRE(conditional + energy_group + entrance_group + waist_group
        == Catch::Approx(total).epsilon(1e-9));
    REQUIRE(energy_group > 0.0);
    REQUIRE(entrance_group > 0.0);
    REQUIRE(waist_group > 0.0);

    // Empty meter and M2 columns: those sources were not supplied.
    REQUIRE(uncertainty.calibration_variance_contributions[0][1] == 0.0);
    REQUIRE(uncertainty.calibration_variance_contributions[0][4] == 0.0);
}

TEST_CASE("entrance-chain-only calibration reaches the propagation")
{
    // Regression: the propagation gate used to ignore the entrance chain
    // when it was the only supplied source, silently returning
    // conditional-only SEs.
    auto request = make_two_curve_request({ pulse_energies[0], pulse_energies[1] });
    request.calibration.entrance_chain_std = 0.006;
    const auto result = zfit::fit(request);

    REQUIRE(result.uncertainty.has_value());
    const auto& uncertainty = *result.uncertainty;
    REQUIRE(uncertainty.total_standard_errors.size() == 1);
    REQUIRE(uncertainty.total_standard_errors[0] > uncertainty.standard_errors[0]);
    REQUIRE(uncertainty.calibration_variance_contributions[0][2] > 0.0);
}

TEST_CASE("fit propagates fixed-parameter priors into total SEs")
{
    // The held-out parameter tau30 (fixed at truth) carries a prior; the
    // group-7 column behaves like any shared calibration direction.
    auto request = make_two_curve_request({ pulse_energies[0], pulse_energies[1] });
    request.calibration.fixed_parameters
        = { { zfit::FitParameter::Tau30, 2.0e-9 } };
    const auto with_prior = zfit::fit(request);

    REQUIRE(with_prior.uncertainty.has_value());
    const auto& uncertainty = *with_prior.uncertainty;
    REQUIRE(uncertainty.total_standard_errors.size() == 1);
    REQUIRE(uncertainty.total_standard_errors[0] > uncertainty.standard_errors[0]);
    REQUIRE(uncertainty.total_covariance_row_major.size() == 1);

    // The identity extends to eight groups: conditional plus every
    // contribution, with only the fixed-parameter group supplied here.
    double contributions = 0.0;
    for (std::size_t g = 0; g < uncertainty.calibration_variance_contributions[0].size();
        ++g) {
        contributions += uncertainty.calibration_variance_contributions[0][g];
    }
    REQUIRE(contributions > 0.0);
    for (std::size_t g = 0; g + 1
        < uncertainty.calibration_variance_contributions[0].size(); ++g) {
        REQUIRE(uncertainty.calibration_variance_contributions[0][g] == 0.0);
    }
    const double conditional = uncertainty.standard_errors[0] * uncertainty.standard_errors[0];
    const double total
        = uncertainty.total_standard_errors[0] * uncertainty.total_standard_errors[0];
    REQUIRE(conditional + contributions == Catch::Approx(total).epsilon(1e-9));

    // A wider prior moves the total SE monotonically.
    auto wider = request;
    wider.calibration.fixed_parameters = { { zfit::FitParameter::Tau30, 4.0e-9 } };
    const auto wide = zfit::fit(wider);
    REQUIRE(wide.uncertainty.has_value());
    REQUIRE(wide.uncertainty->total_standard_errors[0]
        > uncertainty.total_standard_errors[0]);
}

TEST_CASE("a zero fixed-parameter prior claims the value exact")
{
    auto request = make_two_curve_request({ pulse_energies[0], pulse_energies[1] });
    request.calibration.fixed_parameters = { { zfit::FitParameter::Tau30, 0.0 } };
    const auto exact = zfit::fit(request);

    auto uncalibrated = make_two_curve_request({ pulse_energies[0], pulse_energies[1] });
    const auto reference = zfit::fit(uncalibrated);

    REQUIRE(exact.uncertainty.has_value());
    REQUIRE(reference.uncertainty.has_value());
    REQUIRE(exact.uncertainty->total_covariance_row_major.empty());
    REQUIRE(exact.uncertainty->total_standard_errors.empty());
    REQUIRE(exact.uncertainty->standard_errors[0]
        == Catch::Approx(reference.uncertainty->standard_errors[0]).epsilon(1e-12));
}

TEST_CASE("fixed-parameter priors are validated")
{
    auto request = make_two_curve_request({ pulse_energies[0], pulse_energies[1] });

    auto negative = request;
    negative.calibration.fixed_parameters
        = { { zfit::FitParameter::Tau30, -1.0e-9 } };
    REQUIRE_THROWS_AS(zfit::fit(negative), std::invalid_argument);

    auto fitted_collision = request;
    fitted_collision.calibration.fixed_parameters
        = { { zfit::FitParameter::Sigma_G, 1.0e-18 } };
    REQUIRE_THROWS_AS(zfit::fit(fitted_collision), std::invalid_argument);

    auto duplicate = request;
    duplicate.calibration.fixed_parameters
        = { { zfit::FitParameter::Tau30, 2.0e-9 }, { zfit::FitParameter::Tau30, 1.0e-9 } };
    REQUIRE_THROWS_AS(zfit::fit(duplicate), std::invalid_argument);
}

TEST_CASE("fit propagation matches seeded coherent Monte Carlo refits")
{
    // Instrument-level and per-curve calibration errors, shared ones fully
    // correlated across both curves.
    constexpr double energy_stds[2] = { 2.0e-9, 2.0e-9 };
    constexpr double meter_std = 1.5e-9;
    constexpr double entrance_rel_std = 0.01;
    constexpr double waist_std = 0.725e-6;
    constexpr double m2_std = 0.03;
    constexpr double concentration_std = 2.9e-6;
    constexpr double concentration = sample.concentration_mol_per_l;
    constexpr double path_length_std = 1.0e-5;

    auto request = make_two_curve_request({ pulse_energies[0], pulse_energies[1] });
    request.calibration.pulse_energy_stds_j = { energy_stds[0], energy_stds[1] };
    request.calibration.energy_meter_std_j = meter_std;
    request.calibration.entrance_chain_std = entrance_rel_std;
    request.calibration.waist_std_m = waist_std;
    request.calibration.m2_std = m2_std;
    request.calibration.concentration_std = concentration_std;
    request.calibration.path_length_std_m = path_length_std;

    const auto linearized = zfit::fit(request);
    REQUIRE(linearized.uncertainty.has_value());
    REQUIRE(linearized.uncertainty->total_standard_errors.size() == 1);

    // Monte Carlo target: the calibration-induced part of the sigma_g
    // variance, which the observations here lack entirely (they are the
    // noiseless truth curve; only the model inputs are drawn).
    const double conditional_variance
        = linearized.uncertainty->standard_errors[0] * linearized.uncertainty->standard_errors[0];
    const double total_variance = linearized.uncertainty->total_standard_errors[0]
        * linearized.uncertainty->total_standard_errors[0];
    const double calibration_variance = total_variance - conditional_variance;
    REQUIRE(calibration_variance > 0.0);

    // 48 draws: the seeded variance estimate over 24 draws carries ~29%
    // relative sampling error on its own, which is uncomfortably close to
    // the acceptance window; the doubled count halves it.
    constexpr int draw_count = 48;
    std::mt19937_64 engine(20241024);
    std::normal_distribution<double> normal(0.0, 1.0);

    std::vector<double> draws;
    draws.reserve(draw_count);
    for (int d = 0; d < draw_count; ++d) {
        const double meter_draw = meter_std * normal(engine);
        // The entrance chain is a shared multiplicative factor on every
        // curve's energy (relative standard deviation).
        const double entrance_draw = 1.0 + entrance_rel_std * normal(engine);
        const double waist_draw = waist + waist_std * normal(engine);
        const double m2_draw = base_laser.m2 + m2_std * normal(engine);
        const double concentration_draw = concentration + concentration_std * normal(engine);
        const double thickness_draw = sample.thickness_m + path_length_std * normal(engine);

        auto drawn_request = request;
        drawn_request.calibration = { }; // only sigma_g draws are compared
        for (std::size_t k = 0; k < drawn_request.curves.size(); ++k) {
            drawn_request.curves[k].laser.pulse_energy_j
                = pulse_energies[k] * entrance_draw + energy_stds[k] * normal(engine) + meter_draw;
            drawn_request.curves[k].laser.waist_m = waist_draw;
            drawn_request.curves[k].laser.m2 = m2_draw;
        }
        drawn_request.sample.concentration_mol_per_l = concentration_draw;
        drawn_request.sample.thickness_m = thickness_draw;

        const auto drawn = zfit::fit(drawn_request);
        REQUIRE(drawn.converged);
        draws.push_back(drawn.parameters.sigma_g);
    }

    const double mean
        = std::accumulate(draws.begin(), draws.end(), 0.0) / static_cast<double>(draw_count);
    double mc_variance = 0.0;
    for (double draw : draws) {
        mc_variance += (draw - mean) * (draw - mean);
    }
    mc_variance /= static_cast<double>(draw_count - 1);

    const double ratio = mc_variance / calibration_variance;

    CHECK(ratio > 0.7);
    CHECK(ratio < 1.45);
}
