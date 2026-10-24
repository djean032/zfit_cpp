#include "zfit/workup/ratio.hpp"

#include <cmath>
#include <stdexcept>

namespace zfit::workup {

RatioSummary summarize_ratios(const RawScanPoint& point, const DarkSummary& dark)
{
    if (point.readings.size() < 2 || !std::isfinite(dark.ai0_mean)
        || !std::isfinite(dark.ai1_mean)) {
        throw std::invalid_argument("invalid ratio inputs");
    }

    double mean = 0.0;
    double m2 = 0.0;
    double count = 0.0;
    double dmean_d_dark_ai0 = 0.0;
    double dmean_d_dark_ai1 = 0.0;

    for (const auto& reading : point.readings) {
        if (!std::isfinite(reading.ai0) || !std::isfinite(reading.ai1)) {
            throw std::invalid_argument("nonfinite detector reading");
        }

        const double reference = reading.ai0 - dark.ai0_mean;
        const double signal = reading.ai1 - dark.ai1_mean;
        if (!std::isfinite(reference) || reference <= 0.0 || !std::isfinite(signal)) {
            throw std::invalid_argument("invalid dark-corrected reading");
        }

        const double ratio = signal / reference;
        if (!std::isfinite(ratio)) {
            throw std::invalid_argument("nonfinite detector ratio");
        }

        count += 1.0;
        const double delta = ratio - mean;
        const double dr_d_dark_ai0 = ratio / reference;
        const double dr_d_dark_ai1 = -1.0 / reference;
        mean += delta / count;
        m2 += delta * (ratio - mean);
        dmean_d_dark_ai0 += (dr_d_dark_ai0 - dmean_d_dark_ai0) / count;
        dmean_d_dark_ai1 += (dr_d_dark_ai1 - dmean_d_dark_ai1) / count;
    }

    const double mean_variance = m2 / (count * (count - 1.0));
    if (!std::isfinite(mean) || !std::isfinite(mean_variance) || mean_variance < 0.0
        || !std::isfinite(dmean_d_dark_ai0) || !std::isfinite(dmean_d_dark_ai1)) {
        throw std::invalid_argument("invalid ratio statistics");
    }

    return {
        .mean = mean,
        .mean_variance = mean_variance,
        .dmean_d_dark_ai0 = dmean_d_dark_ai0,
        .dmean_d_dark_ai1 = dmean_d_dark_ai1,
    };
}
}
