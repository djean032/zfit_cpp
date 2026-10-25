#pragma once

#include <cstddef>
#include <span>
#include <string>

namespace zfit::workup {

// Per-curve multi-shot energymeter statistics: the incident pulse energy the
// fit consumes, with this curve's own standard error of the mean (an
// independent per-curve uncertainty). The energymeter's instrument-level
// calibration uncertainty is a shared error recorded on SampleCalibration,
// never folded in here.
struct EnergyCalibration {
    double mean_j { };
    double standard_error_j { };
    std::size_t reading_count { };
    std::string provenance_id;
};

// Welford statistics over a multi-shot energy burst: the mean is what the
// fit consumes; the standard error of the mean is this curve's energy
// uncertainty.
EnergyCalibration summarize_energy_measurement(
    std::span<const double> pulse_energies_j, std::string provenance_id);

} // namespace zfit::workup
