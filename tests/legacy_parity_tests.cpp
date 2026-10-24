#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>

#include "zfit/cvodes_integrator.hpp"

namespace {

struct ReferencePoint {
    double z_m { };
    double transmission { };
};

std::vector<ReferencePoint> load_legacy_reference()
{
    std::ifstream input { std::string { ZFIT_SOURCE_DIR }
        + "/tests/data/legacy_zscan_reference.csv" };
    if (!input.good()) {
        throw std::runtime_error("failed to open legacy reference fixture");
    }

    std::string line;
    std::getline(input, line);

    std::vector<ReferencePoint> reference;
    while (std::getline(input, line)) {
        std::stringstream row { line };
        std::string z_m;
        std::string transmission;

        std::getline(row, z_m, ',');
        std::getline(row, transmission);
        reference.push_back({ std::stod(z_m), std::stod(transmission) });
    }

    return reference;
}

void apply_legacy_wing_shift(std::vector<double>& transmission)
{
    constexpr std::size_t wing_points = 5;
    REQUIRE(transmission.size() >= 2 * wing_points);

    double wing_sum { };
    for (std::size_t index = 0; index < wing_points; ++index) {
        wing_sum += transmission[index];
        wing_sum += transmission[transmission.size() - wing_points + index];
    }

    const double wing_average = wing_sum / (2.0 * wing_points);
    for (double& value : transmission) {
        value += 1.0 - wing_average;
    }
}

} // namespace

namespace {

struct ParityFixture {
    zfit::ModelParameters model_parameters {
        .sigma_g = 4.94e-18,
        .sigma_s = 1.60e-17,
        .sigma_t = 1.50e-17,
        .tau10 = 1.00e-12,
        .tau21 = 1.00e-12,
        .tau13 = 1.00e-12,
        .tau30 = 1.29e-7,
        .tau43 = 1.00e-12,
    };

    zfit::SampleParameters sample_parameters {
        // The legacy solver treats its sample_width input as centimeters.
        .thickness_m = 1.0e-3,
        .concentration_mol_per_l = 1.75e18 * 1000.0 / 6.02214076e23,
    };

    zfit::LaserParameters laser_parameters {
        .wavelength_m = 532.0e-9,
        .waist_m = 14.5e-6,
        .pulse_energy_j = 50.0e-9,
        .pulse_fwhm_s = 9.0e-9,
        .m2 = 1.2,
    };
};

struct ParityCurve {
    std::vector<double> transmission;
};

ParityCurve simulate_parity_curve(const ParityFixture& fixture,
    const std::vector<double>& scan_positions_m, zfit::SolverOptions solver_options)
{
    zfit::MolModel model { fixture.model_parameters, fixture.sample_parameters,
        zfit::Laser { fixture.laser_parameters }, solver_options, scan_positions_m };
    zfit::CvodesIntegrator integrator { model };

    ParityCurve curve;
    curve.transmission.reserve(scan_positions_m.size());
    for (std::size_t index = 0; index < scan_positions_m.size(); ++index) {
        curve.transmission.push_back(integrator.simulate(index));
    }

    apply_legacy_wing_shift(curve.transmission);
    return curve;
}

} // namespace

TEST_CASE("CVODES matches the legacy fixed-parameter Z-scan shape")
{
    const auto reference = load_legacy_reference();

    std::vector<double> scan_positions_m;
    scan_positions_m.reserve(reference.size());
    for (const auto& point : reference) {
        scan_positions_m.push_back(point.z_m);
    }

    const ParityFixture fixture;
    const zfit::SolverOptions solver_options {
        .nz = 11,
        .population_rtol = 1.0e-6,
        .quadrature_rtol = 1.0e-7,
        .atol = 1.0e-10,
        .max_steps = 100000,
        .pulse_start_fwhm = -0.5,
        .pulse_end_fwhm = 0.5,
    };

    const auto curve = simulate_parity_curve(fixture, scan_positions_m, solver_options);
    const auto& transmission = curve.transmission;

    REQUIRE(transmission.size() == reference.size());

    double maximum_difference { };
    for (std::size_t index = 0; index < reference.size(); ++index) {
        maximum_difference = std::max(
            maximum_difference, std::abs(transmission[index] - reference[index].transmission));
    }

    INFO("maximum legacy-curve difference: " << maximum_difference);
    REQUIRE(maximum_difference <= 0.02);
}

TEST_CASE("CVODES decoupled tolerances match the tight-tolerance curve")
{
    const auto reference = load_legacy_reference();

    std::vector<double> scan_positions_m;
    scan_positions_m.reserve(reference.size());
    for (const auto& point : reference) {
        scan_positions_m.push_back(point.z_m);
    }

    const ParityFixture fixture;

    const zfit::SolverOptions decoupled_options {
        .nz = 11,
        .population_rtol = 1.0e-6,
        .quadrature_rtol = 1.0e-7,
        .atol = 1.0e-10,
        .max_steps = 100000,
        .pulse_start_fwhm = -0.5,
        .pulse_end_fwhm = 0.5,
    };

    zfit::SolverOptions tight_options = decoupled_options;
    tight_options.population_rtol = 1.0e-8;
    tight_options.quadrature_rtol = 1.0e-8;
    tight_options.atol = 1.0e-12;

    const auto decoupled = simulate_parity_curve(fixture, scan_positions_m, decoupled_options);
    const auto tight = simulate_parity_curve(fixture, scan_positions_m, tight_options);

    REQUIRE(decoupled.transmission.size() == tight.transmission.size());

    double maximum_difference { };
    for (std::size_t index = 0; index < tight.transmission.size(); ++index) {
        maximum_difference = std::max(maximum_difference,
            std::abs(decoupled.transmission[index] - tight.transmission[index]));
    }

    INFO("maximum decoupled-vs-tight difference: " << maximum_difference);
    REQUIRE(maximum_difference <= 1.0e-6);
}
