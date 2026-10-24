#include <cmath>
#include <numbers>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/laser.hpp"

namespace {

constexpr double speed_of_light_m_per_s = 299'792'458.0;
constexpr double planck_constant_j_s = 6.62607015e-34;

constexpr zfit::LaserParameters parameters {
    .wavelength_m = 800.0e-9,
    .waist_m = 20.0e-6,
    .pulse_energy_j = 1.0e-6,
    .pulse_fwhm_s = 100.0e-15,
};

} // namespace

TEST_CASE("laser photon energy and inverse are consistent")
{
    const zfit::Laser laser { parameters };
    const double expected_energy
        = planck_constant_j_s * speed_of_light_m_per_s / parameters.wavelength_m;

    REQUIRE(laser.photon_energy_j == Catch::Approx(expected_energy));
    REQUIRE(laser.photon_energy_j * laser.inverse_photon_energy_per_j == Catch::Approx(1.0));
}

TEST_CASE("laser temporal Gaussian follows the FWHM convention")
{
    const zfit::Laser laser { parameters };

    REQUIRE(laser.temporal_factor(0.0) == Catch::Approx(1.0));
    REQUIRE(laser.temporal_factor(parameters.pulse_fwhm_s / 2.0) == Catch::Approx(0.5));
}

TEST_CASE("laser beam radius grows by square root two at Rayleigh range")
{
    const zfit::Laser laser { parameters };

    REQUIRE(laser.beam_radius_m(0.0) == Catch::Approx(parameters.waist_m));
    REQUIRE(laser.beam_radius_m(laser.rayleigh_range_m)
        == Catch::Approx(parameters.waist_m * std::sqrt(2.0)));
}

TEST_CASE("laser irradiance applies temporal and radial Gaussian factors")
{
    const zfit::Laser laser { parameters };
    const double z_m = 0.5 * laser.rayleigh_range_m;
    const double peak = laser.peak_irradiance_w_per_cm2(z_m);
    const double radius_m = laser.beam_radius_m(z_m);

    REQUIRE(laser.irradiance_w_per_cm2(0.0, 0.0, z_m) == Catch::Approx(peak));
    REQUIRE(laser.irradiance_w_per_cm2(radius_m, parameters.pulse_fwhm_s / 2.0, z_m)
        == Catch::Approx(peak * 0.5 * std::exp(-2.0)));
}
