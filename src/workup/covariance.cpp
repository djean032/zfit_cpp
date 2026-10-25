#include "zfit/workup/covariance.hpp"

#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <vector>

namespace zfit::workup {

double MeasurementCovariance::operator[](std::size_t i, std::size_t j) const
{
    if (i >= point_count || j >= point_count
        || point_count > std::numeric_limits<std::size_t>::max() / point_count
        || values_row_major.size() != point_count * point_count) {
        throw std::out_of_range("invalid covariance index or storage");
    }

    return values_row_major[i * point_count + j];
}

double& MeasurementCovariance::operator[](std::size_t i, std::size_t j)
{
    if (i >= point_count || j >= point_count
        || point_count > std::numeric_limits<std::size_t>::max() / point_count
        || values_row_major.size() != point_count * point_count) {
        throw std::out_of_range("invalid covariance index or storage");
    }

    return values_row_major[i * point_count + j];
}

MeasurementCovariance build_measurement_covariance(
    const PreparedCurve& curve, const DarkSummary& dark, const CuvettePairCalibration& cuvettes)
{
    const std::size_t n = curve.transmission.size();

    if (n == 0 || n > std::numeric_limits<std::size_t>::max() / n
        || curve.scan_positions_m.size() != n || curve.transmission_standard_errors.size() != n
        || curve.dark_sensitivities.size() != n || curve.dark_record_id.empty()
        || curve.dark_record_id != dark.record_id || curve.session_id != dark.session_id
        || curve.cuvette_calibration_id.empty()
        || curve.cuvette_calibration_id != cuvettes.calibration_id
        || !std::isfinite(cuvettes.correction_factor) || cuvettes.correction_factor <= 0.0
        || !cuvettes.standard_error || !std::isfinite(*cuvettes.standard_error)
        || *cuvettes.standard_error < 0.0 || !std::isfinite(dark.ai0_mean_variance)
        || !std::isfinite(dark.ai1_mean_variance) || !std::isfinite(dark.mean_covariance)
        || dark.ai0_mean_variance < 0.0 || dark.ai1_mean_variance < 0.0) {
        throw std::invalid_argument("invalid covariance inputs");
    }

    const double dark_correlation_limit
        = std::sqrt(dark.ai0_mean_variance) * std::sqrt(dark.ai1_mean_variance);

    if (std::abs(dark.mean_covariance) > dark_correlation_limit) {
        throw std::invalid_argument("invalid dark covariance");
    }

    const double relative_c_error = *cuvettes.standard_error / cuvettes.correction_factor;
    if (!std::isfinite(relative_c_error)) {
        throw std::invalid_argument("invalid cuvette relative uncertainty");
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(curve.scan_positions_m[i]) || !std::isfinite(curve.transmission[i])
            || curve.transmission[i] <= 0.0 || !std::isfinite(curve.transmission_standard_errors[i])
            || curve.transmission_standard_errors[i] < 0.0
            || !std::isfinite(curve.dark_sensitivities[i][0])
            || !std::isfinite(curve.dark_sensitivities[i][1])) {
            throw std::invalid_argument("invalid prepared curve point");
        }
    }

    MeasurementCovariance result {
        .point_count = n,
        .values_row_major = std::vector<double>(n * n, 0.0),
    };

    for (std::size_t i = 0; i < n; ++i) {
        const auto& a = curve.dark_sensitivities[i];

        for (std::size_t j = i; j < n; ++j) {
            const auto& b = curve.dark_sensitivities[j];

            const double dark_term
                = a[0] * (dark.ai0_mean_variance * b[0] + dark.mean_covariance * b[1])
                + a[1] * (dark.mean_covariance * b[0] + dark.ai1_mean_variance * b[1]);

            const double cuvette_term = curve.transmission[i] * curve.transmission[j]
                * relative_c_error * relative_c_error;

            double value = dark_term + cuvette_term;
            if (i == j) {
                const double se = curve.transmission_standard_errors[i];
                value += se * se;
            }

            if (!std::isfinite(value) || (i == j && value < 0.0)) {
                throw std::invalid_argument("invalid measurement covariance");
            }

            result[i, j] = value;
            result[j, i] = value;
        }
    }
    return result;
}

CurveSetMeasurementCovariance assemble_curve_set_covariance(std::span<const PreparedCurve> curves,
    std::span<const DarkSummary> darks, std::span<const CuvettePairCalibration> cuvettes)
{
    const std::size_t curve_count = curves.size();
    if (curve_count == 0 || darks.size() != curve_count || cuvettes.size() != curve_count) {
        throw std::invalid_argument("invalid curve set covariance inputs");
    }

    // Per-curve blocks come through build_measurement_covariance, which
    // validates each curve against its own dark and cuvette records.
    std::vector<MeasurementCovariance> blocks;
    blocks.reserve(curve_count);

    CurveSetMeasurementCovariance result;
    result.curve_point_counts.reserve(curve_count);

    std::vector<std::size_t> offsets;
    offsets.reserve(curve_count + 1);
    offsets.push_back(0);

    std::size_t total = 0;
    for (std::size_t k = 0; k < curve_count; ++k) {
        blocks.push_back(build_measurement_covariance(curves[k], darks[k], cuvettes[k]));
        const std::size_t points = blocks.back().point_count;
        if (total > std::numeric_limits<std::size_t>::max() - points
            || total + points > std::numeric_limits<std::size_t>::max() / (total + points)) {
            throw std::invalid_argument("invalid curve set covariance inputs");
        }
        result.curve_point_counts.push_back(points);
        total += points;
        offsets.push_back(total);
    }

    result.joint.point_count = total;
    result.joint.values_row_major.assign(total * total, 0.0);

    for (std::size_t k = 0; k < curve_count; ++k) {
        const std::size_t offset_k = offsets[k];
        const auto& block = blocks[k];
        for (std::size_t i = 0; i < block.point_count; ++i) {
            for (std::size_t j = 0; j < block.point_count; ++j) {
                result.joint[offset_k + i, offset_k + j] = block[i, j];
            }
        }
    }

    // Cross-curve blocks from a shared cuvette correction. Identity rule:
    // equal cuvette calibration ids mean one characterized constant, and the
    // records must then agree exactly.
    for (std::size_t k = 0; k < curve_count; ++k) {
        for (std::size_t l = k + 1; l < curve_count; ++l) {
            const auto& record_k = cuvettes[k];
            const auto& record_l = cuvettes[l];
            if (record_k.calibration_id != record_l.calibration_id) {
                continue;
            }

            if (record_k.correction_factor != record_l.correction_factor
                || record_k.standard_error != record_l.standard_error) {
                throw std::invalid_argument("contradictory shared cuvette calibration records");
            }

            const double relative_error_k = *record_k.standard_error / record_k.correction_factor;
            const double relative_error_l = *record_l.standard_error / record_l.correction_factor;

            const std::size_t offset_k = offsets[k];
            const std::size_t offset_l = offsets[l];
            for (std::size_t i = 0; i < result.curve_point_counts[k]; ++i) {
                for (std::size_t j = 0; j < result.curve_point_counts[l]; ++j) {
                    const double cross = relative_error_k * relative_error_l
                        * curves[k].transmission[i] * curves[l].transmission[j];
                    if (!std::isfinite(cross)) {
                        throw std::invalid_argument("invalid shared cuvette covariance");
                    }

                    result.joint[offset_k + i, offset_l + j] = cross;
                    result.joint[offset_l + j, offset_k + i] = cross;
                }
            }
        }
    }

    return result;
}

}
