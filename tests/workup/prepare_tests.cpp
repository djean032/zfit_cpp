#include <cmath>
#include <limits>
#include <stdexcept>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/workup/prepare.hpp"

namespace {

struct Fixture {
    zfit::workup::RawScan sample;
    zfit::workup::RawScan blank;
    zfit::workup::CuvettePairCalibration cuvettes;
    zfit::workup::DarkSelection dark;
};

Fixture make_fixture()
{
    Fixture fixture;
    auto& sample = fixture.sample;
    sample.kind = zfit::workup::ScanKind::Sample;
    sample.session_id = "session-1";
    sample.sample_id = "sample-1";
    sample.solvent_id = "solvent-1";
    sample.cuvette_id = "sample-cuvette";
    sample.acquisition_id = "sample-scan";
    sample.instrument_configuration_id = "energy-1-filters-1";
    sample.measured_pulse_energy_j = 50.0e-9;
    sample.num_readings_per_point = 2;
    // With darks (ai0, ai1) = (0.1, 0.2), these sample pulse ratios are
    // (1.6, 2.0) and (1.0, 1.2).
    sample.points = {
        { .z_m = -0.001, .readings = { { .ai0 = 1.1, .ai1 = 1.8 }, { .ai0 = 1.1, .ai1 = 2.2 } } },
        { .z_m = 0.001, .readings = { { .ai0 = 1.1, .ai1 = 1.2 }, { .ai0 = 1.1, .ai1 = 1.4 } } },
    };

    fixture.blank = sample;
    auto& blank = fixture.blank;
    blank.kind = zfit::workup::ScanKind::Blank;
    blank.sample_id.clear();
    blank.cuvette_id = "blank-cuvette";
    blank.acquisition_id = "blank-scan";
    // Blank pulse ratios are (1.8, 2.2) and (2.0, 2.0).
    blank.points = {
        { .z_m = -0.001, .readings = { { .ai0 = 1.1, .ai1 = 2.0 }, { .ai0 = 1.1, .ai1 = 2.4 } } },
        { .z_m = 0.001, .readings = { { .ai0 = 1.1, .ai1 = 2.2 }, { .ai0 = 1.1, .ai1 = 2.2 } } },
    };

    fixture.cuvettes.calibration_id = "cuvette-calibration-1";
    fixture.cuvettes.sample_cuvette_id = sample.cuvette_id;
    fixture.cuvettes.blank_cuvette_id = blank.cuvette_id;
    fixture.cuvettes.instrument_configuration_id = sample.instrument_configuration_id;
    fixture.cuvettes.correction_factor = 1.1;
    fixture.cuvettes.standard_error = 0.05;

    auto& dark = fixture.dark;
    dark.selected.record_id = "dark-1";
    dark.selected.session_id = sample.session_id;
    dark.selected.instrument_configuration_id = sample.instrument_configuration_id;
    dark.selected.ai0_mean = 0.1;
    dark.selected.ai1_mean = 0.2;
    dark.selected.ai0_mean_variance = 0.01;
    dark.selected.ai1_mean_variance = 0.04;
    dark.selected.mean_covariance = 0.02;
    dark.drift_checked = true;
    return fixture;
}

} // namespace

TEST_CASE("workup prepares absolute transmission and reading-noise SEs")
{
    const auto fixture = make_fixture();
    const auto curve = zfit::workup::prepare_curve(
        fixture.sample, fixture.blank, fixture.cuvettes, fixture.dark);

    REQUIRE(curve.scan_positions_m.size() == 2);
    REQUIRE(curve.dark_sensitivities.size() == curve.scan_positions_m.size());
    REQUIRE(curve.scan_positions_m[0] == -0.001);
    REQUIRE(curve.scan_positions_m[1] == 0.001);
    // Raw ai1/ai0 readings exceed one; matched blank and cuvette correction
    // still recover physical transmission below one.
    REQUIRE(
        fixture.sample.points[0].readings[0].ai1 / fixture.sample.points[0].readings[0].ai0 > 1.0);
    REQUIRE(curve.transmission[0] == Catch::Approx(1.8 / (2.0 * 1.1)));
    REQUIRE(curve.transmission[1] == Catch::Approx(1.1 / (2.0 * 1.1)));

    // The first sample and blank ratio means each have variance 0.04;
    // only the sample ratio varies at the second point (variance 0.01).
    const double first_T = 1.8 / (2.0 * 1.1);
    const double first_variance
        = 0.04 / (2.0 * 1.1 * 2.0 * 1.1) + 0.04 * std::pow(first_T / 2.0, 2);
    REQUIRE(curve.transmission_standard_errors[0] == Catch::Approx(std::sqrt(first_variance)));
    REQUIRE(curve.transmission_standard_errors[1] == Catch::Approx(0.1 / (2.0 * 1.1)));

    REQUIRE(curve.session_id == fixture.sample.session_id);
    REQUIRE(curve.sample_id == fixture.sample.sample_id);
    REQUIRE(curve.sample_acquisition_id == fixture.sample.acquisition_id);
    REQUIRE(curve.blank_acquisition_id == fixture.blank.acquisition_id);
    REQUIRE(curve.dark_record_id == fixture.dark.selected.record_id);
    REQUIRE(curve.cuvette_calibration_id == fixture.cuvettes.calibration_id);
    REQUIRE(curve.dark_drift_checked);
    REQUIRE(curve.measured_pulse_energy_j == fixture.sample.measured_pulse_energy_j);
}

TEST_CASE("workup prepared dark sensitivities match finite differences")
{
    constexpr double step = 1.0e-5;
    auto fixture = make_fixture();
    // Give one sample pulse a different reference level so the ai0 effect
    // does not cancel exactly against the blank at the first scan position.
    fixture.sample.points[0].readings[1] = { .ai0 = 2.1, .ai1 = 4.2 };
    const auto baseline = zfit::workup::prepare_curve(
        fixture.sample, fixture.blank, fixture.cuvettes, fixture.dark);
    REQUIRE(baseline.dark_sensitivities.size() == baseline.transmission.size());
    REQUIRE(std::abs(baseline.dark_sensitivities[0][0]) > 0.1);
    REQUIRE(std::abs(baseline.dark_sensitivities[0][1]) > 0.01);

    for (std::size_t channel = 0; channel < 2; ++channel) {
        auto plus = fixture.dark;
        auto minus = fixture.dark;
        if (channel == 0) {
            plus.selected.ai0_mean += step;
            minus.selected.ai0_mean -= step;
        } else {
            plus.selected.ai1_mean += step;
            minus.selected.ai1_mean -= step;
        }

        const auto curve_plus
            = zfit::workup::prepare_curve(fixture.sample, fixture.blank, fixture.cuvettes, plus);
        const auto curve_minus
            = zfit::workup::prepare_curve(fixture.sample, fixture.blank, fixture.cuvettes, minus);
        for (std::size_t i = 0; i < baseline.transmission.size(); ++i) {
            const double finite_difference
                = (curve_plus.transmission[i] - curve_minus.transmission[i]) / (2.0 * step);
            REQUIRE(baseline.dark_sensitivities[i][channel]
                == Catch::Approx(finite_difference).margin(1.0e-8));
        }
    }
}

TEST_CASE("workup keeps shared dark and cuvette uncertainty out of point SEs")
{
    auto fixture = make_fixture();
    const auto baseline = zfit::workup::prepare_curve(
        fixture.sample, fixture.blank, fixture.cuvettes, fixture.dark);
    fixture.dark.selected.ai0_mean_variance = 100.0;
    fixture.dark.selected.ai1_mean_variance = 100.0;
    fixture.cuvettes.standard_error = 0.5;
    const auto revised = zfit::workup::prepare_curve(
        fixture.sample, fixture.blank, fixture.cuvettes, fixture.dark);
    REQUIRE(revised.transmission_standard_errors == baseline.transmission_standard_errors);
}

TEST_CASE("workup retains measured transmission slightly above one")
{
    auto fixture = make_fixture();
    fixture.sample.points[0].readings = {
        { .ai0 = 1.1, .ai1 = 2.6 },
        { .ai0 = 1.1, .ai1 = 2.8 },
    };
    const auto curve = zfit::workup::prepare_curve(
        fixture.sample, fixture.blank, fixture.cuvettes, fixture.dark);
    REQUIRE(curve.transmission[0] > 1.0);
}

TEST_CASE("workup rejects an invalid dark selection or pulse energy")
{
    auto fixture = make_fixture();
    fixture.dark.selected.session_id = "another-session";
    REQUIRE_THROWS_AS(
        zfit::workup::prepare_curve(fixture.sample, fixture.blank, fixture.cuvettes, fixture.dark),
        std::invalid_argument);

    fixture = make_fixture();
    fixture.sample.measured_pulse_energy_j = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS(
        zfit::workup::prepare_curve(fixture.sample, fixture.blank, fixture.cuvettes, fixture.dark),
        std::invalid_argument);
}

TEST_CASE("workup rejects a nonpositive sample or blank mean ratio")
{
    auto fixture = make_fixture();
    fixture.blank.points[0].readings = {
        { .ai0 = 1.1, .ai1 = 0.2 },
        { .ai0 = 1.1, .ai1 = 0.2 },
    };
    REQUIRE_THROWS_AS(
        zfit::workup::prepare_curve(fixture.sample, fixture.blank, fixture.cuvettes, fixture.dark),
        std::invalid_argument);

    fixture = make_fixture();
    fixture.sample.points[0].readings = {
        { .ai0 = 1.1, .ai1 = 0.2 },
        { .ai0 = 1.1, .ai1 = 0.2 },
    };
    REQUIRE_THROWS_AS(
        zfit::workup::prepare_curve(fixture.sample, fixture.blank, fixture.cuvettes, fixture.dark),
        std::invalid_argument);
}
