#pragma once

#include "scan.hpp"

namespace zfit::workup {

void validate_pair(
    const RawScan& sample, const RawScan& blank, const CuvettePairCalibration& cuvettes);
} // namespace zfit::workup
