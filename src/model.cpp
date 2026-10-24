#include "../include/zfit/model.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace zfit {

namespace {

    constexpr double avogadro_per_mol = 6.02214076e23;
    constexpr double liters_to_cm3 = 1000.0;
    constexpr double meters_to_cm = 100.0;

    template <typename ValueAt>
    double cubic_hermite(std::span<const double> times, double time, ValueAt value_at) noexcept
    {
        if (times.size() < 2 || !std::isfinite(time) || time < times.front()
            || time > times.back()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        if (time == times.back()) {
            return value_at(times.size() - 1);
        }

        const auto upper = std::upper_bound(times.begin(), times.end(), time);
        const std::size_t right = static_cast<std::size_t>(upper - times.begin());
        if (right == 0 || right >= times.size()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        const std::size_t left = right - 1;
        const auto slope_at = [&](std::size_t knot) noexcept {
            if (knot == 0) {
                return (value_at(1) - value_at(0)) / (times[1] - times[0]);
            }

            if (knot + 1 == times.size()) {
                return (value_at(knot) - value_at(knot - 1)) / (times[knot] - times[knot - 1]);
            }

            return (value_at(knot + 1) - value_at(knot - 1)) / (times[knot + 1] - times[knot - 1]);
        };

        const double interval = times[right] - times[left];
        const double u = (time - times[left]) / interval;
        const double u2 = u * u;
        const double u3 = u2 * u;

        const double y0 = value_at(left);
        const double y1 = value_at(right);
        const double m0 = slope_at(left);
        const double m1 = slope_at(right);

        return (2.0 * u3 - 3.0 * u2 + 1.0) * y0 + (u3 - 2.0 * u2 + u) * interval * m0
            + (-2.0 * u3 + 3.0 * u2) * y1 + (u3 - u2) * interval * m1;
    }

} // namespace

bool PulseWaveform::valid() const noexcept
{
    if (times.size() < 2 || times.size() != log_intensity.size()
        || times.size() > std::numeric_limits<std::size_t>::max() / fit_parameter_count
        || log_sensitivities.size() != times.size() * fit_parameter_count) {
        return false;
    }

    for (std::size_t knot = 0; knot < times.size(); ++knot) {
        if (!std::isfinite(times[knot]) || !std::isfinite(log_intensity[knot])
            || (knot > 0 && !(times[knot] > times[knot - 1]))) {
            return false;
        }
    }

    for (double sensitivity : log_sensitivities) {
        if (!std::isfinite(sensitivity)) {
            return false;
        }
    }

    return true;
}

double PulseWaveform::log_sensitivity(std::size_t knot, std::size_t parameter) const noexcept
{
    if (parameter >= fit_parameter_count || knot >= times.size()
        || log_intensity.size() != times.size()
        || times.size() > std::numeric_limits<std::size_t>::max() / fit_parameter_count
        || log_sensitivities.size() != times.size() * fit_parameter_count) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return log_sensitivities[knot * fit_parameter_count + parameter];
}

MolModel::MolModel(ModelParameters model_params, SampleParameters sample_params, Laser laser,
    SolverOptions solver_options, std::span<const double> scan_positions_m)
    : model_params(model_params)
    , sample_params(sample_params)
    , laser(std::move(laser))
    , solver_options(solver_options)
    , dz_m(sample_params.thickness_m / static_cast<double>(solver_options.nz))
    , dz_cm(dz_m * meters_to_cm)
    , number_density_per_cm3(
          sample_params.concentration_mol_per_l * avogadro_per_mol / liters_to_cm3)
    , inverse_photon_energy_per_j(laser.inverse_photon_energy_per_j)
{
    incident_peak_w_per_cm2.reserve(scan_positions_m.size());

    for (double z_m : scan_positions_m) {
        incident_peak_w_per_cm2.push_back(laser.peak_irradiance_w_per_cm2(z_m));
    }
}

double PulseWaveform::intensity(double time) const noexcept
{
    if (times.size() < 2 || log_intensity.size() != times.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double log_i = cubic_hermite(
        times, time, [this](std::size_t knot) noexcept { return log_intensity[knot]; });

    if (!std::isfinite(log_i)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double intensity = std::exp(log_i);
    return std::isfinite(intensity) ? intensity : std::numeric_limits<double>::quiet_NaN();
}

double PulseWaveform::intensity_sensitivity(double time, std::size_t parameter) const noexcept
{
    if (parameter >= fit_parameter_count || times.size() < 2 || log_intensity.size() != times.size()
        || times.size() > std::numeric_limits<std::size_t>::max() / fit_parameter_count
        || log_sensitivities.size() != times.size() * fit_parameter_count) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double incoming_intensity = intensity(time);
    if (!std::isfinite(incoming_intensity)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double log_intensity_sensitivity = cubic_hermite(times, time,
        [this, parameter](std::size_t knot) noexcept { return log_sensitivity(knot, parameter); });

    if (!std::isfinite(log_intensity_sensitivity)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double sensitivity = incoming_intensity * log_intensity_sensitivity;
    return std::isfinite(sensitivity) ? sensitivity : std::numeric_limits<double>::quiet_NaN();
}

void MolModel::initialize_cell_state(std::span<double, n_levels> state) const noexcept
{
    state[0] = number_density_per_cm3;
    state[1] = state[2] = state[3] = state[4] = 0.0;
}

double MolModel::incident_irradiance(double t_s, std::size_t scan_index) const noexcept
{
    return incident_peak_w_per_cm2[scan_index] * laser.temporal_factor(t_s);
}

PulseWaveform MolModel::incident_waveform(std::size_t scan_index) const
{
    if (scan_index >= incident_peak_w_per_cm2.size()) {
        throw std::out_of_range("scan index is out of range");
    }

    const std::size_t knot_count = solver_options.waveform_knots;
    if (knot_count < 2
        || knot_count > std::numeric_limits<std::size_t>::max() / fit_parameter_count) {
        throw std::invalid_argument("invalid waveform knot count");
    }

    const double t_start = solver_options.pulse_start_fwhm * laser.params.pulse_fwhm_s;
    const double t_end = solver_options.pulse_end_fwhm * laser.params.pulse_fwhm_s;

    if (!(t_end > t_start) || !std::isfinite(t_start) || !std::isfinite(t_end)) {
        throw std::invalid_argument("invalid waveform time interval");
    }

    PulseWaveform waveform;
    waveform.times.resize(knot_count);
    waveform.log_intensity.resize(knot_count);
    waveform.log_sensitivities.assign(knot_count * fit_parameter_count, 0.0);

    const double dt = (t_end - t_start) / static_cast<double>(knot_count - 1);
    if (!(dt > 0.0) || !std::isfinite(dt)) {
        throw std::invalid_argument("invalid waveform time spacing");
    }
    for (std::size_t knot = 0; knot < knot_count; ++knot) {
        const double time
            = knot + 1 == knot_count ? t_end : t_start + static_cast<double>(knot) * dt;
        const double intensity = incident_irradiance(time, scan_index);

        if (!(intensity > 0.0) || !std::isfinite(intensity)) {
            throw std::invalid_argument("incident waveform intensity is invalid");
        }

        waveform.times[knot] = time;
        waveform.log_intensity[knot] = std::log(intensity);
    }

    return waveform;
};

void MolModel::explicit_rhs(double t_s, std::size_t scan_index, std::span<const double> state,
    std::span<double> derivative) const
{
    assert(state.size() == solver_options.nz * n_levels);
    assert(derivative.size() == solver_options.nz * n_levels);

    double intensity = incident_irradiance(t_s, scan_index);

    for (std::size_t cell = 0; cell < solver_options.nz; ++cell) {
        const std::size_t offset = cell * n_levels;
        const std::span<const double, n_levels> population { state.data() + offset, n_levels };
        const std::span<double, n_levels> rate { derivative.data() + offset, n_levels };

        const double photon_flux = intensity * inverse_photon_energy_per_j;

        optical_rhs(photon_flux, population, rate, model_params);

        const double alpha_per_cm = absorption_coefficient(population, model_params);

        intensity *= std::exp(-alpha_per_cm * dz_cm);
    }
}

double MolModel::explicit_rhs_with_transmission(double t_s, std::size_t scan_index,
    std::span<const double> state, std::span<double> derivative) const
{
    assert(state.size() == solver_options.nz * n_levels);
    assert(derivative.size() == solver_options.nz * n_levels);

    double intensity = incident_irradiance(t_s, scan_index);

    for (std::size_t cell = 0; cell < solver_options.nz; ++cell) {
        const std::size_t offset = cell * n_levels;
        const std::span<const double, n_levels> population { state.data() + offset, n_levels };
        const std::span<double, n_levels> rate { derivative.data() + offset, n_levels };

        const double photon_flux = intensity * inverse_photon_energy_per_j;

        optical_rhs(photon_flux, population, rate, model_params);

        const double alpha_per_cm = absorption_coefficient(population, model_params);

        intensity *= std::exp(-alpha_per_cm * dz_cm);
    }
    return intensity;
}

void MolModel::implicit_rhs(std::span<const double> state, std::span<double> derivative) const
{
    assert(state.size() == solver_options.nz * n_levels);
    assert(derivative.size() == solver_options.nz * n_levels);

    for (std::size_t cell = 0; cell < solver_options.nz; ++cell) {
        const std::size_t offset = cell * n_levels;
        const std::span<const double, n_levels> population { state.data() + offset, n_levels };
        const std::span<double, n_levels> rate { derivative.data() + offset, n_levels };

        relaxation_rhs(population, rate, model_params);
    }
}

void MolModel::implicit_cell_jacobian(
    std::span<double, n_levels * n_levels> jacobian) const noexcept
{
    relaxation_jacobian(jacobian, model_params);
}

double MolModel::transmitted_irradiance(
    double t_s, std::size_t scan_index, std::span<const double> state) const
{
    assert(state.size() == solver_options.nz * n_levels);

    double intensity = incident_irradiance(t_s, scan_index);
    for (std::size_t cell = 0; cell < solver_options.nz; ++cell) {
        const std::size_t offset = cell * n_levels;
        const std::span<const double, n_levels> population { state.data() + offset, n_levels };

        const double alpha_per_cm = absorption_coefficient(population, model_params);

        intensity *= std::exp(-alpha_per_cm * dz_cm);
    }

    return intensity;
}

} // namespace zfit
