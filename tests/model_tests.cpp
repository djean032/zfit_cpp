#include <array>
#include <cmath>
#include <limits>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/model.hpp"

namespace {

constexpr double avogadro_per_mol = 6.02214076e23;

constexpr zfit::ModelParameters model_parameters {
    .sigma_g = 1.0e-18,
    .sigma_s = 2.0e-18,
    .sigma_t = 3.0e-18,
    .tau10 = 2.0e-9,
    .tau21 = 4.0e-9,
    .tau13 = 5.0e-9,
    .tau30 = 10.0e-9,
    .tau43 = 20.0e-9,
};

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

zfit::MolModel make_model(zfit::SolverOptions options = solver_options)
{
    const std::array<double, 2> scan_positions_m { 0.0, 1.0e-3 };
    return zfit::MolModel { model_parameters, sample_parameters, zfit::Laser { laser_parameters },
        options, scan_positions_m };
}

zfit::PulseWaveform make_linear_waveform()
{
    zfit::PulseWaveform waveform;
    waveform.times = { -2.0, -0.5, 1.0, 3.0 };
    waveform.log_sensitivities.assign(waveform.times.size() * zfit::fit_parameter_count, 0.0);

    for (std::size_t knot = 0; knot < waveform.times.size(); ++knot) {
        const double time = waveform.times[knot];
        waveform.log_intensity.push_back(1.0 + 0.5 * time);
        waveform.log_sensitivities[knot * zfit::fit_parameter_count + 2] = -0.25 + 0.75 * time;
    }

    return waveform;
}

} // namespace

TEST_CASE("model derives spatial and density quantities")
{
    const auto model = make_model();

    REQUIRE(model.incident_peak_w_per_cm2.size() == 2);
    REQUIRE(model.dz_m == Catch::Approx(5.0e-3));
    REQUIRE(model.dz_cm == Catch::Approx(0.5));
    REQUIRE(model.number_density_per_cm3 == Catch::Approx(avogadro_per_mol * 1.0e-6));
    REQUIRE(model.incident_peak_w_per_cm2[0] > 0.0);
    REQUIRE(model.incident_peak_w_per_cm2[1] < model.incident_peak_w_per_cm2[0]);
}

TEST_CASE("model initializes a local cell in the ground state")
{
    const auto model = make_model();
    std::array<double, zfit::n_levels> state { };

    model.initialize_cell_state(state);

    REQUIRE(state[0] == Catch::Approx(model.number_density_per_cm3));
    for (std::size_t level = 1; level < zfit::n_levels; ++level) {
        REQUIRE(state[level] == Catch::Approx(0.0));
    }
}

TEST_CASE("pulse waveform validates storage")
{
    const auto waveform = make_linear_waveform();
    REQUIRE(waveform.valid());

    auto bad_shape = waveform;
    bad_shape.log_intensity.pop_back();
    REQUIRE_FALSE(bad_shape.valid());

    auto unordered = waveform;
    unordered.times[2] = unordered.times[1];
    REQUIRE_FALSE(unordered.valid());

    auto nonfinite = waveform;
    nonfinite.log_sensitivities[0] = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_FALSE(nonfinite.valid());
}

TEST_CASE("pulse waveform interpolates log intensity and sensitivity")
{
    const auto waveform = make_linear_waveform();
    constexpr double time = 0.25;
    const double expected_log_intensity = 1.0 + 0.5 * time;
    const double expected_intensity = std::exp(expected_log_intensity);
    const double expected_log_sensitivity = -0.25 + 0.75 * time;

    REQUIRE(waveform.intensity(time) == Catch::Approx(expected_intensity));
    REQUIRE(waveform.intensity_sensitivity(time, 2)
        == Catch::Approx(expected_intensity * expected_log_sensitivity));
    REQUIRE(waveform.intensity_sensitivity(time, 1) == Catch::Approx(0.0));
}

TEST_CASE("pulse waveform reports invalid evaluations with NaN")
{
    const auto waveform = make_linear_waveform();

    REQUIRE(std::isnan(waveform.intensity(-3.0)));
    REQUIRE(std::isnan(waveform.intensity_sensitivity(0.0, zfit::fit_parameter_count)));

    auto bad_shape = waveform;
    bad_shape.log_intensity.pop_back();
    REQUIRE(std::isnan(bad_shape.intensity(0.0)));
}

TEST_CASE("model constructs a uniform incident waveform")
{
    const auto model = make_model();
    const auto waveform = model.incident_waveform(0);
    const double t_start = model.solver_options.pulse_start_fwhm * model.laser.params.pulse_fwhm_s;
    const double t_end = model.solver_options.pulse_end_fwhm * model.laser.params.pulse_fwhm_s;

    REQUIRE(waveform.valid());
    REQUIRE(waveform.times.size() == model.solver_options.waveform_knots);
    REQUIRE(waveform.times.front() == Catch::Approx(t_start));
    REQUIRE(waveform.times.back() == Catch::Approx(t_end));

    for (std::size_t knot = 0; knot < waveform.times.size(); ++knot) {
        REQUIRE(waveform.intensity(waveform.times[knot])
            == Catch::Approx(model.incident_irradiance(waveform.times[knot], 0)));
        for (std::size_t parameter = 0; parameter < zfit::fit_parameter_count; ++parameter) {
            REQUIRE(waveform.log_sensitivity(knot, parameter) == Catch::Approx(0.0));
        }
    }
}
