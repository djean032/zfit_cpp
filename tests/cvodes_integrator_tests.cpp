#include <array>
#include <cmath>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/cvodes_integrator.hpp"

namespace {

constexpr zfit::SampleParameters sample_parameters {
    .thickness_m = 1.0e-2,
    .concentration_mol_per_l = 1.0e-3,
};

constexpr zfit::LaserParameters laser_parameters {
    .wavelength_m = 800.0e-9,
    .waist_m = 20.0e-6,
    .pulse_energy_j = 1.0e-6,
    .pulse_fwhm_s = 100.0e-15,
};

constexpr zfit::SolverOptions solver_options { .nz = 2 };

constexpr zfit::ModelParameters absorbing_parameters {
    .sigma_g = 1.0e-18,
    .sigma_s = 2.0e-18,
    .sigma_t = 3.0e-18,
    .tau10 = 2.0e-9,
    .tau21 = 4.0e-9,
    .tau13 = 5.0e-9,
    .tau30 = 10.0e-9,
    .tau43 = 20.0e-9,
};

zfit::MolModel make_model(
    zfit::ModelParameters model_parameters, zfit::SolverOptions options = solver_options)
{
    const std::array<double, 1> scan_positions_m { 0.0 };
    return zfit::MolModel { model_parameters, sample_parameters, zfit::Laser { laser_parameters },
        options, scan_positions_m };
}

double simulate_focus(zfit::SolverOptions options)
{
    auto model = make_model(absorbing_parameters, options);
    zfit::CvodesIntegrator integrator { model };
    return integrator.simulate(0);
}

} // namespace

TEST_CASE("CVODES transmission is one for a transparent sample")
{
    auto transparent_parameters = absorbing_parameters;
    transparent_parameters.sigma_g = 0.0;
    transparent_parameters.sigma_s = 0.0;
    transparent_parameters.sigma_t = 0.0;

    auto model = make_model(transparent_parameters);
    zfit::CvodesIntegrator integrator { model };

    REQUIRE(integrator.simulate(0) == Catch::Approx(1.0).epsilon(1.0e-5));
}

TEST_CASE("CVODES resets populations before every independent scan")
{
    auto model = make_model(absorbing_parameters);
    zfit::CvodesIntegrator integrator { model };

    const double first = integrator.simulate(0);
    const double second = integrator.simulate(0);

    REQUIRE(first > 0.0);
    REQUIRE(first < 1.0);
    REQUIRE(second == Catch::Approx(first).epsilon(1.0e-8));
}

TEST_CASE("CVODES transmission converges under spatial refinement")
{
    auto options = solver_options;
    options.nz = 16;
    const double transmission_16 = simulate_focus(options);

    options.nz = 32;
    const double transmission_32 = simulate_focus(options);

    options.nz = 64;
    const double transmission_64 = simulate_focus(options);

    options.nz = 128;
    const double transmission_128 = simulate_focus(options);

    const double difference_16_32 = std::abs(transmission_16 - transmission_32);
    const double difference_32_64 = std::abs(transmission_32 - transmission_64);
    const double difference_64_128 = std::abs(transmission_64 - transmission_128);

    REQUIRE(difference_32_64 < difference_16_32);
    REQUIRE(difference_64_128 < difference_32_64);
}

TEST_CASE("CVODES transmission converges as tolerances tighten")
{
    auto options = solver_options;
    const double default_transmission = simulate_focus(options);

    options.population_rtol = 1.0e-8;
    options.quadrature_rtol = 1.0e-8;
    options.atol = 1.0e-12;
    const double tight_transmission = simulate_focus(options);

    options.population_rtol = 1.0e-10;
    options.quadrature_rtol = 1.0e-10;
    options.atol = 1.0e-14;
    const double reference_transmission = simulate_focus(options);

    const double default_difference = std::abs(default_transmission - reference_transmission);
    const double tight_difference = std::abs(tight_transmission - reference_transmission);

    REQUIRE(tight_difference < default_difference);
}

TEST_CASE("CVODES transmission converges as waveform knots refine")
{
    auto options = solver_options;
    options.waveform_knots = 129;
    const double transmission_129 = simulate_focus(options);

    options.waveform_knots = 257;
    const double transmission_257 = simulate_focus(options);

    options.waveform_knots = 513;
    const double transmission_513 = simulate_focus(options);

    REQUIRE(std::abs(transmission_257 - transmission_513)
        < std::abs(transmission_129 - transmission_257));
}
