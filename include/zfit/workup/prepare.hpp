#pragma once

#include "dark.hpp"
#include "zfit/workup/scan.hpp"

namespace zfit::workup {

PreparedCurve prepare_curve(const RawScan& sample, const RawScan& blank,
    const CuvettePairCalibration& cuvettes, const DarkSelection& dark);

}
