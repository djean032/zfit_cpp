#include "zfit/workup/matching.hpp"
#include "zfit/workup/scan.hpp"

#include <cmath>
#include <ranges>
#include <stdexcept>

namespace zfit::workup {

void validate_pair(
    const RawScan& sample, const RawScan& blank, const CuvettePairCalibration& cuvettes)
{
    if (sample.kind != ScanKind::Sample || blank.kind != ScanKind::Blank) {
        throw std::invalid_argument("expected a sample scan and a blank scan");
    }

    if (sample.acquisition_id.empty() || blank.acquisition_id.empty()
        || sample.acquisition_id == blank.acquisition_id || sample.session_id.empty()
        || sample.session_id != blank.session_id || sample.solvent_id.empty()
        || sample.solvent_id != blank.solvent_id || sample.instrument_configuration_id.empty()
        || sample.instrument_configuration_id != blank.instrument_configuration_id) {
        throw std::invalid_argument("sample and blank do not match");
    }
    if (cuvettes.calibration_id.empty() || sample.cuvette_id != cuvettes.sample_cuvette_id
        || blank.cuvette_id != cuvettes.blank_cuvette_id
        || sample.instrument_configuration_id != cuvettes.instrument_configuration_id
        || !std::isfinite(cuvettes.correction_factor) || cuvettes.correction_factor <= 0.0
        || (cuvettes.standard_error
            && (!std::isfinite(*cuvettes.standard_error) || *cuvettes.standard_error < 0.0))) {
        throw std::invalid_argument("invalid cuvette-pair calibration");
    }
    if (sample.points.empty() || sample.points.size() != blank.points.size()) {
        throw std::invalid_argument("sample and blank z grids differ");
    }
    if (sample.num_readings_per_point < 2 || blank.num_readings_per_point < 2) {
        throw std::invalid_argument("too few readings per point");
    }

    for (const auto& [sample_point, blank_point] : std::views::zip(sample.points, blank.points)) {
        if (!std::isfinite(sample_point.z_m) || !std::isfinite(blank_point.z_m)
            || sample_point.z_m != blank_point.z_m
            || sample_point.readings.size() != sample.num_readings_per_point
            || blank_point.readings.size() != blank.num_readings_per_point) {
            throw std::invalid_argument("invalid or mismatched scan point");
        }
    }
}

} // zfit::workup
