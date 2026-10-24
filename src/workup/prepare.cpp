#include "zfit/workup/prepare.hpp"

#include <cmath>
#include <ranges>
#include <stdexcept>

#include "zfit/workup/matching.hpp"
#include "zfit/workup/ratio.hpp"
#include "zfit/workup/scan.hpp"

namespace zfit::workup {

PreparedCurve prepare_curve(const RawScan& sample, const RawScan& blank,
    const CuvettePairCalibration& cuvettes, const DarkSelection& dark)
{
    validate_pair(sample, blank, cuvettes);

    if (sample.sample_id.empty() || !std::isfinite(sample.measured_pulse_energy_j)
        || sample.measured_pulse_energy_j <= 0.0 || dark.selected.record_id.empty()
        || dark.selected.session_id != sample.session_id
        || dark.selected.instrument_configuration_id != sample.instrument_configuration_id) {
        throw std::invalid_argument("invalid prepared-curve inputs");
    }

    PreparedCurve curve;
    curve.session_id = sample.session_id;
    curve.sample_id = sample.sample_id;
    curve.sample_acquisition_id = sample.acquisition_id;
    curve.blank_acquisition_id = blank.acquisition_id;
    curve.dark_record_id = dark.selected.record_id;
    curve.cuvette_calibration_id = cuvettes.calibration_id;
    curve.dark_drift_checked = dark.drift_checked;
    curve.measured_pulse_energy_j = sample.measured_pulse_energy_j;

    curve.scan_positions_m.reserve(sample.points.size());
    curve.transmission.reserve(sample.points.size());
    curve.transmission_standard_errors.reserve(sample.points.size());
    curve.dark_sensitivities.reserve(sample.points.size());

    for (const auto& [sample_point, blank_point] : std::views::zip(sample.points, blank.points)) {
        const auto S = summarize_ratios(sample_point, dark.selected);
        const auto B = summarize_ratios(blank_point, dark.selected);

        if (S.mean <= 0.0 || B.mean <= 0.0) {
            throw std::invalid_argument("nonpositive sample or blank ratio");
        }

        const double denominator = B.mean * cuvettes.correction_factor;
        const double T = S.mean / denominator;
        if (!std::isfinite(denominator) || denominator <= 0.0 || !std::isfinite(T) || T <= 0.0) {
            throw std::invalid_argument("invalid calibrated transmission");
        }

        const double dT_dS = 1.0 / denominator;
        const double dT_dB = -T / B.mean;
        const double variance = dT_dS * dT_dS * S.mean_variance + dT_dB * dT_dB * B.mean_variance;
        const double dT_d_dark_ai0 = dT_dS * S.dmean_d_dark_ai0 + dT_dB * B.dmean_d_dark_ai0;
        const double dT_d_dark_ai1 = dT_dS * S.dmean_d_dark_ai1 + dT_dB * B.dmean_d_dark_ai1;

        if (!std::isfinite(dT_d_dark_ai0) || !std::isfinite(dT_d_dark_ai1)) {
            throw std::invalid_argument("invalid dark count sensitivities");
        }

        if (!std::isfinite(variance) || variance < 0.0) {
            throw std::invalid_argument("invalid transmission variance");
        }

        curve.scan_positions_m.push_back(sample_point.z_m);
        curve.transmission.push_back(T);
        curve.transmission_standard_errors.push_back(std::sqrt(variance));
        curve.dark_sensitivities.push_back({ dT_d_dark_ai0, dT_d_dark_ai1 });
    }
    return curve;
}
}
