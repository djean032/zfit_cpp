#include <cmath>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/cvodes_integrator.hpp"
#include "zfit/fit_bridge.hpp"
#include "zfit/model.hpp"
#include "zfit/workup/dark.hpp"
#include "zfit/workup/entrance_chain.hpp"
#include "zfit/workup/prepare.hpp"

namespace {

using zfit::workup::CuvettePairCalibration;
using zfit::workup::DarkRecord;
using zfit::workup::DarkSelection;
using zfit::workup::DetectorReading;
using zfit::workup::EnergyCalibration;
using zfit::workup::RawScan;
using zfit::workup::RawScanPoint;
using zfit::workup::SampleCalibration;
using zfit::workup::ScanKind;
using zfit::workup::summarize_dark;

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

constexpr zfit::SampleParameters truth_sample {
    .thickness_m = 1.0e-3,
    .concentration_mol_per_l = 1.75e18 * 1000.0 / 6.02214076e23,
};

constexpr zfit::SolverOptions options {
    .nz = 11,
    .waveform_knots = 129,
    .pulse_start_fwhm = -0.5,
    .pulse_end_fwhm = 0.5,
};

constexpr double cuvette_correction = 1.05;
constexpr double window_index = 1.4610; // fused silica at 532 nm
constexpr double window_index_std = 0.0050; // catalogue tolerance
constexpr double solvent_index = 1.4969; // toluene at 532 nm
constexpr double solvent_index_std = 0.0010;
constexpr double reading_scatter = 0.03;

CuvettePairCalibration make_cuvettes()
{
    return CuvettePairCalibration {
        .calibration_id = "cuvette-campaign",
        .sample_cuvette_id = "cuvette-s",
        .blank_cuvette_id = "cuvette-b",
        .instrument_configuration_id = "config-1",
        .correction_factor = cuvette_correction,
        .standard_error = 0.03,
    };
}

// The entrance chain, computed exactly the way production code will get it:
// raw material indices with catalogue tolerances, Fresnel product computed
// in the workup layer. No hand-typed kappa anywhere.
zfit::workup::EntranceChainValue entrance_chain()
{
    zfit::workup::EntranceChainInputs inputs;
    inputs.window_refractive_index = window_index;
    inputs.window_index_std = window_index_std;
    inputs.solvent_refractive_index = solvent_index;
    inputs.solvent_index_std = solvent_index_std;
    inputs.index_source_id = "index-catalogue";
    return zfit::workup::compute_entrance_chain(inputs);
}

SampleCalibration make_sample_calibration()
{
    SampleCalibration calibration;
    calibration.laser_characterization_id = "laser-campaign";
    calibration.wavelength_m = 532.0e-9;
    calibration.waist_m = 14.5e-6;
    calibration.waist_standard_error_m = 0.725e-6;
    calibration.pulse_fwhm_s = 9.0e-9;
    calibration.m2 = 1.2;
    calibration.m2_standard_error = 0.03;
    calibration.energy_meter_calibration_id = "meter-campaign";
    calibration.energy_meter_standard_error_j = 1.5e-9;
    const auto chain = entrance_chain();
    calibration.entrance_chain_calibration_id = chain.index_source_id;
    calibration.entrance_chain_transmission = chain.transmission;
    calibration.entrance_chain_standard_error = chain.standard_error;
    calibration.concentration_mol_per_l = truth_sample.concentration_mol_per_l;
    calibration.concentration_standard_error_mol_per_l = 2.9e-6;
    calibration.path_length_m = truth_sample.thickness_m;
    calibration.path_length_standard_error_m = 1.0e-5;
    return calibration;
}

zfit::LaserParameters laser_at(double pulse_energy_j)
{
    return zfit::LaserParameters {
        .wavelength_m = 532.0e-9,
        .waist_m = 14.5e-6,
        .pulse_energy_j = pulse_energy_j,
        .pulse_fwhm_s = 9.0e-9,
        .m2 = 1.2,
    };
}

std::vector<double> scan_positions()
{
    std::vector<double> positions;
    for (int i = -10; i <= 10; ++i) {
        positions.push_back(static_cast<double>(i) * 1.0e-3);
    }
    return positions;
}

std::vector<double> simulate_curve(double pulse_energy_j, const std::vector<double>& positions)
{
    zfit::MolModel model { truth_parameters, truth_sample, zfit::Laser { laser_at(pulse_energy_j) },
        options, positions };
    zfit::CvodesIntegrator integrator { model };
    std::vector<double> transmission;
    transmission.reserve(positions.size());
    for (std::size_t i = 0; i < positions.size(); ++i) {
        transmission.push_back(integrator.simulate(i));
    }
    return transmission;
}

// Detector readings whose dark-corrected per-reading ratios bracket the
// target exactly, so the mean ratio is the target with a nonzero scatter.
std::vector<DetectorReading> readings_for_ratio(
    double target_ratio, double reference_counts, const DetectorReading& dark)
{
    const double reference = reference_counts + dark.ai0;
    return {
        DetectorReading { .ai0 = reference,
            .ai1 = target_ratio * (1.0 - reading_scatter) * reference_counts + dark.ai1 },
        DetectorReading { .ai0 = reference,
            .ai1 = target_ratio * (1.0 + reading_scatter) * reference_counts + dark.ai1 },
    };
}

RawScan make_scan(ScanKind kind, double pulse_energy_j, const std::vector<double>& positions,
    const std::vector<double>& target_ratios, const DetectorReading& dark,
    const std::string& acquisition_id)
{
    RawScan scan;
    scan.kind = kind;
    scan.session_id = "session-1";
    scan.sample_id = "sample-1";
    scan.solvent_id = "solvent-1";
    scan.cuvette_id = (kind == ScanKind::Sample) ? "cuvette-s" : "cuvette-b";
    scan.measured_pulse_energy_j = pulse_energy_j;
    scan.acquisition_id = acquisition_id;
    scan.instrument_configuration_id = "config-1";
    scan.num_readings_per_point = 2;
    scan.points.reserve(positions.size());

    // Blank and sample channels read at different detector scales.
    const double reference_counts = (kind == ScanKind::Sample) ? 50000.0 : 60000.0;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        scan.points.push_back(RawScanPoint {
            .z_m = positions[i],
            .readings = readings_for_ratio(target_ratios[i], reference_counts, dark),
        });
    }
    return scan;
}

struct CurveInputs {
    zfit::workup::PreparedCurve curve;
    zfit::workup::DarkSummary dark;
    zfit::workup::CuvettePairCalibration cuvettes;
    EnergyCalibration energy;
};

CurveInputs make_curve_inputs(
    double pulse_energy_j, const EnergyCalibration& energy, const std::string& suffix)
{
    const auto positions = scan_positions();
    // The molecule physically sees the entrance-chain-scaled energy: the
    // synthetic readings must come from that local fluence so the workup
    // output and the locally-scaled fit inputs stay consistent.
    const auto molecular
        = simulate_curve(entrance_chain().transmission * pulse_energy_j, positions);

    std::vector<double> blank_ratios;
    std::vector<double> sample_ratios;
    blank_ratios.reserve(positions.size());
    sample_ratios.reserve(positions.size());
    // Blank ratio with mild structure: a wrong blank or cuvette correction
    // cannot masquerade as a molecular transmission dip.
    const double blank_cycle[3] = { 0.88, 0.92, 0.95 };
    for (std::size_t i = 0; i < positions.size(); ++i) {
        const double blank = blank_cycle[i % 3];
        blank_ratios.push_back(blank);
        sample_ratios.push_back(molecular[i] * blank * cuvette_correction);
    }

    const DetectorReading dark_reading { 20.0, 30.0 };
    DarkRecord dark_record {
        .record_id = "dark-" + suffix,
        .session_id = "session-1",
        .instrument_configuration_id = "config-1",
        .acquired_at = { },
        .readings = { dark_reading, dark_reading },
    };

    CurveInputs inputs;
    inputs.energy = energy;
    inputs.cuvettes = make_cuvettes();
    const DarkSelection dark_selection { summarize_dark(dark_record), true };
    inputs.curve
        = zfit::workup::prepare_curve(make_scan(ScanKind::Sample, pulse_energy_j, positions,
                                          sample_ratios, dark_reading, "sample-" + suffix),
            make_scan(ScanKind::Blank, pulse_energy_j, positions, blank_ratios, dark_reading,
                "blank-" + suffix),
            inputs.cuvettes, dark_selection);
    inputs.dark = dark_selection.selected;
    return inputs;
}

EnergyCalibration burst(std::vector<double> energies, const std::string& id)
{
    return zfit::workup::summarize_energy_measurement(energies, id);
}

struct Fixture {
    CurveInputs first;
    CurveInputs second;
    std::vector<zfit::workup::PreparedCurve> curves;
    std::vector<EnergyCalibration> energies;
    std::vector<zfit::workup::DarkSummary> darks;
    std::vector<CuvettePairCalibration> cuvettes;
    SampleCalibration calibration;
};

Fixture make_two_curve_fixture()
{
    Fixture fixture;
    fixture.first = make_curve_inputs(50.0e-9, burst({ 48.5e-9, 51.5e-9 }, "meter-burst-50"), "50");
    fixture.second
        = make_curve_inputs(80.0e-9, burst({ 78.0e-9, 82.0e-9 }, "meter-burst-80"), "80");

    fixture.curves = { fixture.first.curve, fixture.second.curve };
    fixture.energies = { fixture.first.energy, fixture.second.energy };
    fixture.darks = { fixture.first.dark, fixture.second.dark };
    fixture.cuvettes = { fixture.first.cuvettes, fixture.second.cuvettes };
    fixture.calibration = make_sample_calibration();
    return fixture;
}

} // namespace

TEST_CASE("bridge couples workup curves into a joint fit request")
{
    const Fixture fixture = make_two_curve_fixture();
    const auto request = zfit::make_fit_request(
        fixture.curves, fixture.energies, fixture.darks, fixture.cuvettes, fixture.calibration);

    REQUIRE(request.curves.size() == 2);
    REQUIRE(request.sample.thickness_m == Catch::Approx(fixture.calibration.path_length_m));
    REQUIRE(request.sample.concentration_mol_per_l
        == Catch::Approx(fixture.calibration.concentration_mol_per_l));

    for (std::size_t k = 0; k < request.curves.size(); ++k) {
        const double expected_energy = (k == 0) ? 50.0e-9 : 80.0e-9;
        const auto chain = entrance_chain();
        const auto& laser = request.curves[k].laser;
        // Fit energies are physically local: the energymeter burst mean
        // scaled by the entrance chain.
        REQUIRE(laser.pulse_energy_j == Catch::Approx(chain.transmission * expected_energy));
        REQUIRE(laser.wavelength_m == Catch::Approx(532.0e-9));
        REQUIRE(laser.waist_m == Catch::Approx(14.5e-6));
        REQUIRE(laser.pulse_fwhm_s == Catch::Approx(9.0e-9));
        REQUIRE(laser.m2 == Catch::Approx(1.2));

        // Workup output reproduces the local-fluence molecular transmission
        // of the burst mean, so these observations must match the physics.
        const auto molecular
            = simulate_curve(chain.transmission * expected_energy, scan_positions());
        for (std::size_t i = 0; i < molecular.size(); ++i) {
            REQUIRE(request.curves[k].observed_transmission[i]
                == Catch::Approx(molecular[i]).epsilon(1e-9));
        }
        REQUIRE(request.curves[k].measurement_covariance_row_major.empty());
    }

    // Calibration uncertainty block: per-curve burst SEM (scaled through
    // the local energies) plus the shared instrument-level errors, all
    // forwarded.
    REQUIRE(request.calibration.pulse_energy_stds_j.size() == 2);
    REQUIRE(request.calibration.pulse_energy_stds_j[0]
        == Catch::Approx(entrance_chain().transmission * 1.5e-9));
    REQUIRE(request.calibration.pulse_energy_stds_j[1]
        == Catch::Approx(entrance_chain().transmission * 2.0e-9));
    REQUIRE(request.calibration.energy_meter_std_j == Catch::Approx(1.5e-9));
    REQUIRE(request.calibration.entrance_chain_std
        == Catch::Approx(entrance_chain().standard_error / entrance_chain().transmission));
    REQUIRE(request.calibration.waist_std_m == Catch::Approx(0.725e-6));
    REQUIRE(request.calibration.m2_std == Catch::Approx(0.03));
    REQUIRE(request.calibration.concentration_std == Catch::Approx(2.9e-6));
    REQUIRE(request.calibration.path_length_std_m == Catch::Approx(1.0e-5));

    // Joint covariance: 42 x 42 with shared-cuvette cross blocks.
    const std::size_t points = scan_positions().size();
    const std::size_t total = 2 * points;
    REQUIRE(request.joint_measurement_covariance_row_major.size() == total * total);
    const double relative = 0.03 / cuvette_correction;
    REQUIRE(request.joint_measurement_covariance_row_major[0 * total + points]
        == Catch::Approx(relative * relative * fixture.curves[0].transmission[0]
            * fixture.curves[1].transmission[0]));
    // Curve 1's last diagonal entry: its own reading noise plus the cuvette
    // term of its own transmission.
    REQUIRE(request.joint_measurement_covariance_row_major[(total - 1) * total + (total - 1)]
        == Catch::Approx(relative * relative * fixture.curves[1].transmission[points - 1]
                * fixture.curves[1].transmission[points - 1]
            + fixture.curves[1].transmission_standard_errors[points - 1]
                * fixture.curves[1].transmission_standard_errors[points - 1])
            .epsilon(1e-9));
}

TEST_CASE("bridge refuses uncalibrated or inconsistent inputs")
{
    const Fixture fixture = make_two_curve_fixture();
    const auto build = [&fixture](const SampleCalibration& calibration) {
        return zfit::make_fit_request(
            fixture.curves, fixture.energies, fixture.darks, fixture.cuvettes, calibration);
    };

    auto calibration = make_sample_calibration();
    calibration.waist_standard_error_m = 0.0;
    REQUIRE_THROWS_AS(build(calibration), std::invalid_argument);

    calibration = make_sample_calibration();
    calibration.m2_standard_error = 0.0;
    REQUIRE_THROWS_AS(build(calibration), std::invalid_argument);

    calibration = make_sample_calibration();
    calibration.energy_meter_standard_error_j = 0.0;
    REQUIRE_THROWS_AS(build(calibration), std::invalid_argument);

    calibration = make_sample_calibration();
    calibration.concentration_standard_error_mol_per_l = 0.0;
    REQUIRE_THROWS_AS(build(calibration), std::invalid_argument);

    calibration = make_sample_calibration();
    calibration.path_length_standard_error_m = 0.0;
    REQUIRE_THROWS_AS(build(calibration), std::invalid_argument);

    calibration = make_sample_calibration();
    calibration.entrance_chain_standard_error = 0.0;
    REQUIRE_THROWS_AS(build(calibration), std::invalid_argument);

    calibration = make_sample_calibration();
    calibration.entrance_chain_calibration_id.clear();
    REQUIRE_THROWS_AS(build(calibration), std::invalid_argument);

    calibration = make_sample_calibration();
    calibration.entrance_chain_transmission = 1.2;
    REQUIRE_THROWS_AS(build(calibration), std::invalid_argument);

    calibration = make_sample_calibration();
    calibration.entrance_chain_transmission = 0.0;
    REQUIRE_THROWS_AS(build(calibration), std::invalid_argument);

    calibration = make_sample_calibration();
    calibration.laser_characterization_id.clear();
    REQUIRE_THROWS_AS(build(calibration), std::invalid_argument);

    calibration = make_sample_calibration();
    calibration.energy_meter_calibration_id.clear();
    REQUIRE_THROWS_AS(build(calibration), std::invalid_argument);

    {
        std::vector<zfit::workup::PreparedCurve> one_curve { fixture.curves[0] };
        std::vector<EnergyCalibration> one_energy { fixture.energies[0] };
        REQUIRE_THROWS_AS(zfit::make_fit_request(one_curve, fixture.energies, fixture.darks,
                              fixture.cuvettes, fixture.calibration),
            std::invalid_argument);

        REQUIRE_THROWS_AS(zfit::make_fit_request(fixture.curves, one_energy, fixture.darks,
                              fixture.cuvettes, fixture.calibration),
            std::invalid_argument);
    }

    {
        auto malformed_energy = fixture.energies[0];
        malformed_energy.mean_j = 0.0;
        std::vector<EnergyCalibration> energies { malformed_energy, fixture.energies[1] };
        REQUIRE_THROWS_AS(zfit::make_fit_request(fixture.curves, energies, fixture.darks,
                              fixture.cuvettes, fixture.calibration),
            std::invalid_argument);
    }

    {
        auto missing_provenance = fixture.energies[0];
        missing_provenance.provenance_id.clear();
        std::vector<EnergyCalibration> energies { missing_provenance, fixture.energies[1] };
        REQUIRE_THROWS_AS(zfit::make_fit_request(fixture.curves, energies, fixture.darks,
                              fixture.cuvettes, fixture.calibration),
            std::invalid_argument);
    }
}

TEST_CASE("bridge fits detector-readings-derived curves end to end")
{
    const Fixture fixture = make_two_curve_fixture();
    auto request = zfit::make_fit_request(
        fixture.curves, fixture.energies, fixture.darks, fixture.cuvettes, fixture.calibration);

    request.initial_parameters = truth_parameters;
    request.initial_parameters.sigma_g *= 2.0;
    request.solver_options = options;
    request.fitted_parameters = { zfit::FitParameter::Sigma_G };

    const auto result = zfit::fit(request);

    REQUIRE(result.converged);
    REQUIRE((result.status >= 1 && result.status <= 4));
    REQUIRE(
        result.parameters.sigma_g / truth_parameters.sigma_g == Catch::Approx(1.0).margin(0.02));
    REQUIRE(result.chi_squared.has_value());
    REQUIRE(result.uncertainty.has_value());
    REQUIRE(result.uncertainty->standard_errors.size() == 1);
    // The fit passes through the shared-cuvette joint covariance, so any
    // squeezed uncertainty must reflect that correlation: nothing here may
    // pretend the cuvette correction was exactly known.
    REQUIRE(result.uncertainty->covariance_row_major[0] > 0.0);
}
