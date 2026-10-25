#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "dark.hpp"
#include "scan.hpp"

namespace zfit::workup {

struct MeasurementCovariance {
    std::size_t point_count { };
    std::vector<double> values_row_major;

    double operator[](std::size_t i, std::size_t j) const;
    double& operator[](std::size_t i, std::size_t j);
};

MeasurementCovariance build_measurement_covariance(
    const PreparedCurve& curve, const DarkSummary& dark, const CuvettePairCalibration& cuvettes);

// Joint covariance over several curves of one sample. Per-curve blocks come
// from build_measurement_covariance (reading noise, within-curve dark
// projection, within-curve cuvette term). Wherever two curves reference the
// **same characterized cuvette correction** (equal calibration ids), the
// shared-C uncertainty adds cross-curve blocks
// (dT_ki/dC)(dT_lj/dC) sigma_C^2 = relC_k relC_l T_ki T_lj, because
// T = S/(B C) makes dT/dC = -T/C. Curves with different calibration ids share
// nothing.
struct CurveSetMeasurementCovariance {
    std::vector<std::size_t> curve_point_counts; // aligned with the input curves
    MeasurementCovariance joint; // total point count
};

CurveSetMeasurementCovariance assemble_curve_set_covariance(std::span<const PreparedCurve> curves,
    std::span<const DarkSummary> darks, std::span<const CuvettePairCalibration> cuvettes);
} // namespace zfit::workup
