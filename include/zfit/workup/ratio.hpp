#pragma once

#include "dark.hpp"
#include "zfit/workup/scan.hpp"

namespace zfit::workup {

struct RatioSummary {
    double mean { };
    double mean_variance { };
    double dmean_d_dark_ai0 { };
    double dmean_d_dark_ai1 { };
};

RatioSummary summarize_ratios(const RawScanPoint& point, const DarkSummary& dark);
}
