#pragma once

#include <cstddef>
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
}
