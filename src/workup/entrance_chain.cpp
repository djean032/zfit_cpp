#include "zfit/workup/entrance_chain.hpp"

#include <cmath>
#include <stdexcept>

namespace zfit::workup {
namespace {

    double fresnel_amplitude_reflection(double from_index, double to_index)
    {
        // Normal-incidence amplitude reflection: (n2 - n1)/(n1 + n2).
        const double amplitude = (to_index - from_index) / (from_index + to_index);
        return amplitude;
    }

    double transmission_from_lower(
        double lower_index, double upper_index, double amplitude_reflection)
    {
        // dT/da with the upper index fixed. Since T = 1 - r^2 with r the
        // amplitude reflection, dT/da = 4 * n_upper * r / (n_lower + n_upper)^2:
        // raising the lower index toward the upper one reduces the reflection.
        const double sum = lower_index + upper_index;
        return 4.0 * upper_index * amplitude_reflection / (sum * sum);
    }

    double transmission_from_upper(
        double lower_index, double upper_index, double amplitude_reflection)
    {
        // dT/db with the lower index fixed: dT/db = -4 * n_lower * r / (sum)^2.
        // Raising the upper index increases the contrast and lowers T.
        const double sum = lower_index + upper_index;
        return -4.0 * lower_index * amplitude_reflection / (sum * sum);
    }

    void throw_invalid() { throw std::invalid_argument("invalid entrance chain inputs"); }

} // namespace

EntranceChainValue compute_entrance_chain(const EntranceChainInputs& inputs)
{
    for (double value : { inputs.window_refractive_index, inputs.solvent_refractive_index,
             inputs.window_index_std, inputs.solvent_index_std }) {
        if (!std::isfinite(value)) {
            throw_invalid();
        }
    }
    if (inputs.window_refractive_index <= 0.0 || inputs.solvent_refractive_index <= 0.0
        || inputs.window_index_std <= 0.0 || inputs.solvent_index_std <= 0.0
        || inputs.index_source_id.empty()) {
        throw_invalid();
    }

    const double air_window_amplitude
        = fresnel_amplitude_reflection(1.0, inputs.window_refractive_index);
    const double window_solvent_amplitude = fresnel_amplitude_reflection(
        inputs.window_refractive_index, inputs.solvent_refractive_index);
    const double reflection_air_window = air_window_amplitude * air_window_amplitude;
    const double reflection_window_solvent = window_solvent_amplitude * window_solvent_amplitude;
    const double transmission_air_window = 1.0 - reflection_air_window;
    const double transmission_window_solvent = 1.0 - reflection_window_solvent;
    const double transmission = transmission_air_window * transmission_window_solvent;
    if (!std::isfinite(transmission) || transmission <= 0.0 || transmission > 1.0) {
        throw std::invalid_argument("entrance chain transmission out of range");
    }

    // d(kappa)/d(window index): the window index is the upper index of the
    // air->window surface and the lower index of the window->solvent surface.
    const double kappa_dw
        = transmission_from_upper(1.0, inputs.window_refractive_index, air_window_amplitude)
            * transmission_window_solvent
        + transmission_air_window
            * transmission_from_lower(inputs.window_refractive_index,
                inputs.solvent_refractive_index, window_solvent_amplitude);
    // d(kappa)/d(solvent index): the solvent index is only the upper index
    // of the window->solvent surface.
    const double kappa_ds = transmission_air_window
        * transmission_from_upper(inputs.window_refractive_index, inputs.solvent_refractive_index,
            window_solvent_amplitude);

    const double variance = kappa_dw * kappa_dw * inputs.window_index_std * inputs.window_index_std
        + kappa_ds * kappa_ds * inputs.solvent_index_std * inputs.solvent_index_std;
    const double standard_error = std::sqrt(variance);
    if (!std::isfinite(standard_error) || standard_error <= 0.0) {
        throw std::invalid_argument("entrance chain standard error out of range");
    }

    return EntranceChainValue {
        .transmission = transmission,
        .standard_error = standard_error,
        .index_source_id = inputs.index_source_id,
    };
}

} // namespace zfit::workup
