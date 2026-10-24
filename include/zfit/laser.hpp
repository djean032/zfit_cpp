#pragma once

#include <cmath>
#include <numbers>

namespace zfit {

struct LaserParameters {
    double wavelength_m { };
    double waist_m { };
    double pulse_energy_j { };
    double pulse_fwhm_s { };
    double m2 { 1.0 };
};

class Laser {
public:
    explicit Laser(LaserParameters params) noexcept;

    [[nodiscard]]
    double temporal_factor(double t_s) const noexcept;

    [[nodiscard]]
    double beam_radius_m(double z_m) const noexcept;

    [[nodiscard]]
    double peak_irradiance_w_per_cm2(double z_m) const noexcept;

    [[nodiscard]]
    double pulse_power_w(double t_s) const noexcept;

    [[nodiscard]]
    double irradiance_w_per_cm2(double radius_m, double t_s, double z_m) const noexcept;

    [[nodiscard]]
    double photon_flux_per_cm2_s(double radius_m, double t_s, double z_m) const noexcept;

    LaserParameters params;

    double rayleigh_range_m;
    double frequency_hz;
    double photon_energy_j;
    double inverse_photon_energy_per_j;

    double peak_power_w { };
    double temporal_decay_per_s2 { };
};

} // namespace zfit
