#pragma once

#include "scan.hpp"
#include <chrono>
#include <span>

namespace zfit::workup {

struct DarkSummary {
    std::string record_id;
    std::string session_id;
    std::string instrument_configuration_id;
    std::chrono::system_clock::time_point acquired_at;

    double ai0_mean { };
    double ai1_mean { };
    double ai0_mean_variance { };
    double ai1_mean_variance { };
    double mean_covariance { };
};

struct DarkDriftLimits {
    double max_ai0_change { };
    double max_ai1_change { };
};

struct DarkSelection {
    DarkSummary selected;
    bool drift_checked { };
};

DarkSelection select_dark(
    const RawScan& scan, std::span<const DarkRecord> darks, DarkDriftLimits limits);

DarkSelection select_dark_for_pair(const RawScan& sample, const RawScan& blank,
    std::span<const DarkRecord> darks, DarkDriftLimits limits);

DarkSummary summarize_dark(const DarkRecord& record);

}
