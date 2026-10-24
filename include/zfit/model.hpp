#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "kinetics.hpp"
#include "laser.hpp"
#include "parameters.hpp"

namespace zfit {

struct SampleParameters {
    double thickness_m { };
    double concentration_mol_per_l { };
};

struct SolverOptions {
    std::size_t nz { 64 };

    double population_rtol { 1e-6 };
    double quadrature_rtol { 1e-7 };
    double atol { 1e-10 };
    std::size_t max_steps { 100000 };
    std::size_t waveform_knots { 513 };

    double pulse_start_fwhm { -4.0 };
    double pulse_end_fwhm { 4.0 };
};

struct PulseWaveform {
    std::vector<double> times;
    std::vector<double> log_intensity;
    std::vector<double> log_sensitivities;

    double log_sensitivity(std::size_t knot, std::size_t parameter) const noexcept;
    double intensity(double time) const noexcept;
    double intensity_sensitivity(double time, std::size_t parameter) const noexcept;
    [[nodiscard]] bool valid() const noexcept;
};

class MolModel {
public:
    MolModel(ModelParameters model_params, SampleParameters sample_params, Laser laser,
        SolverOptions solver_options, std::span<const double> scan_positions_m);

    [[nodiscard]]
    double transmitted_irradiance(
        double t_s, std::size_t scan_index, std::span<const double> state) const;

    void initialize_cell_state(std::span<double, n_levels> state) const noexcept;

    /*  void rhs(double t_s, double scan_z_m, std::span<const double> state,
               std::span<double> derivative) const;*/

    void explicit_rhs(double t_s, std::size_t scan_index, std::span<const double> state,
        std::span<double> derivative) const;

    [[nodiscard]]
    double explicit_rhs_with_transmission(double t_s, std::size_t scan_index,
        std::span<const double> state, std::span<double> derivative) const;

    void implicit_rhs(std::span<const double> state, std::span<double> derivative) const;

    void implicit_cell_jacobian(std::span<double, n_levels * n_levels> jacobian) const noexcept;

    [[nodiscard]]
    double incident_irradiance(double t_s, std::size_t scan_index) const noexcept;

    [[nodiscard]]
    PulseWaveform incident_waveform(std::size_t scan_index) const;

    ModelParameters model_params;
    SampleParameters sample_params;
    Laser laser;
    SolverOptions solver_options;

    double dz_m { };
    double dz_cm { };
    double number_density_per_cm3 { };
    double inverse_photon_energy_per_j { };

    std::vector<double> incident_peak_w_per_cm2;
};

} // namespace zfit
