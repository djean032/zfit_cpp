#include "../include/zfit/laser.hpp"

namespace zfit {

namespace {
    constexpr double speed_of_light_m_per_s = 299'792'458.0;
    constexpr double planck_constant_j_s = 6.62607015e-34;

    constexpr double m2_to_cm2 = 1.0e4;
    constexpr double per_m2_to_per_cm2 = 1.0e-4;
} // namespace

Laser::Laser(LaserParameters params) noexcept
    : params(params)
{
    rayleigh_range_m
        = std::numbers::pi * params.waist_m * params.waist_m / (params.m2 * params.wavelength_m);

    frequency_hz = speed_of_light_m_per_s / params.wavelength_m;

    photon_energy_j = planck_constant_j_s * frequency_hz;

    inverse_photon_energy_per_j = 1.0 / photon_energy_j;

    const double ln2 = std::log(2.0);

    const double gaussian_prefactor = 2.0 * std::sqrt(ln2 / std::numbers::pi);

    peak_power_w = gaussian_prefactor * params.pulse_energy_j / params.pulse_fwhm_s;

    temporal_decay_per_s2 = 4.0 * ln2 / (params.pulse_fwhm_s * params.pulse_fwhm_s);
}

double Laser::temporal_factor(double t_s) const noexcept
{
    return std::exp(-temporal_decay_per_s2 * t_s * t_s);
}

double Laser::beam_radius_m(double z_m) const noexcept
{
    const double ratio = z_m / rayleigh_range_m;

    return params.waist_m * std::sqrt(1.0 + ratio * ratio);
}

double Laser::pulse_power_w(double t_s) const noexcept
{
    return peak_power_w * temporal_factor(t_s);
}

double Laser::irradiance_w_per_cm2(double radius_m, double t_s, double z_m) const noexcept
{
    const double wz = beam_radius_m(z_m);

    const double radial_factor = std::exp(-2.0 * radius_m * radius_m / (wz * wz));

    const double irradiance_w_per_cm2
        = peak_irradiance_w_per_cm2(z_m) * temporal_factor(t_s) * radial_factor;

    return irradiance_w_per_cm2;
}

double Laser::photon_flux_per_cm2_s(double radius_m, double t_s, double z_m) const noexcept
{
    return irradiance_w_per_cm2(radius_m, t_s, z_m) * inverse_photon_energy_per_j;
}

double Laser::peak_irradiance_w_per_cm2(double z_m) const noexcept
{
    const double wz = beam_radius_m(z_m);

    const double irradiance_w_per_m2 = (2.0 * peak_power_w) / (std::numbers::pi * wz * wz);

    return irradiance_w_per_m2 * per_m2_to_per_cm2;
}

} // namespace zfit
