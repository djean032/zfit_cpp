#include "zfit/workup/energy.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace zfit::workup {

EnergyCalibration summarize_energy_measurement(
    std::span<const double> pulse_energies_j, std::string provenance_id)
{
    if (pulse_energies_j.size() < 2 || provenance_id.empty()) {
        throw std::invalid_argument("invalid energy measurement inputs");
    }

    double mean = 0.0;
    double m2 = 0.0;
    double count = 0.0;

    for (const double energy : pulse_energies_j) {
        if (!std::isfinite(energy) || energy <= 0.0) {
            throw std::invalid_argument("invalid pulse energy reading");
        }

        count += 1.0;
        const double delta = energy - mean;
        mean += delta / count;
        m2 += delta * (energy - mean);
    }

    const double mean_variance = m2 / (count * (count - 1.0));
    if (!std::isfinite(mean) || !std::isfinite(mean_variance) || mean_variance < 0.0) {
        throw std::invalid_argument("invalid energy statistics");
    }

    return {
        .mean_j = mean,
        .standard_error_j = std::sqrt(mean_variance),
        .reading_count = pulse_energies_j.size(),
        .provenance_id = std::move(provenance_id),
    };
}

} // namespace zfit::workup
