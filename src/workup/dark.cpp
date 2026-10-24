#include "zfit/workup/dark.hpp"
#include "zfit/workup/scan.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>

namespace zfit::workup {

DarkSummary summarize_dark(const DarkRecord& record)
{
    if (record.record_id.empty() || record.instrument_configuration_id.empty()
        || record.readings.size() < 2) {
        throw std::invalid_argument("invalid dark record");
    }

    if (record.session_id.empty()) {
        throw std::invalid_argument("invalid session id");
    }

    double mean0 = 0.0;
    double mean1 = 0.0;
    double m20 = 0.0;
    double m21 = 0.0;
    double m201 = 0.0;
    double count = 0.0;

    for (const auto& reading : record.readings) {
        if (!std::isfinite(reading.ai0) || !std::isfinite(reading.ai1)) {
            throw std::invalid_argument("nonfinite dark reading");
        }

        count += 1.0;
        const double delta0 = reading.ai0 - mean0;
        const double delta1 = reading.ai1 - mean1;
        mean0 += delta0 / count;
        mean1 += delta1 / count;
        m20 += delta0 * (reading.ai0 - mean0);
        m21 += delta1 * (reading.ai1 - mean1);
        m201 += delta0 * (reading.ai1 - mean1);
    }

    if (!std::isfinite(mean0) || !std::isfinite(mean1) || !std::isfinite(m20) || !std::isfinite(m21)
        || !std::isfinite(m201)) {
        throw std::invalid_argument("nonfinite dark mean, variance, or covariance");
    }

    const double denominator = count * (count - 1.0);
    return {
        .record_id = record.record_id,
        .session_id = record.session_id,
        .instrument_configuration_id = record.instrument_configuration_id,
        .acquired_at = record.acquired_at,
        .ai0_mean = mean0,
        .ai1_mean = mean1,
        .ai0_mean_variance = m20 / denominator,
        .ai1_mean_variance = m21 / denominator,
        .mean_covariance = m201 / denominator,
    };
}

namespace {

    DarkSelection select_dark_in_interval(const std::string& session_id,
        const std::string& configuration_id, std::chrono::system_clock::time_point first,
        std::chrono::system_clock::time_point last, std::span<const DarkRecord> darks,
        DarkDriftLimits limits)
    {
        if (session_id.empty() || configuration_id.empty() || last < first
            || !std::isfinite(limits.max_ai0_change) || !std::isfinite(limits.max_ai1_change)
            || limits.max_ai0_change < 0.0 || limits.max_ai1_change < 0.0) {
            throw std::invalid_argument("invalid dark selection inputs");
        }

        const auto midpoint = first + (last - first) / 2;
        const DarkRecord* closest = nullptr;
        const DarkRecord* before = nullptr;
        const DarkRecord* after = nullptr;

        for (const auto& record : darks) {
            if (record.session_id != session_id
                || record.instrument_configuration_id != configuration_id) {
                continue;
            }

            const auto gap = record.acquired_at < midpoint ? midpoint - record.acquired_at
                                                           : record.acquired_at - midpoint;
            if (closest == nullptr) {
                closest = &record;
            } else {
                const auto closest_gap = closest->acquired_at < midpoint
                    ? midpoint - closest->acquired_at
                    : closest->acquired_at - midpoint;
                if (gap < closest_gap
                    || (gap == closest_gap && record.acquired_at < closest->acquired_at)) {
                    closest = &record;
                }
            }

            if (record.acquired_at <= first
                && (before == nullptr || record.acquired_at > before->acquired_at)) {
                before = &record;
            }
            // A dark at a single scan's timestamp is "before" only. For a pair,
            // a dark at the later scan time can bracket the entire interval.
            if ((record.acquired_at > last || (first < last && record.acquired_at == last))
                && (after == nullptr || record.acquired_at < after->acquired_at)) {
                after = &record;
            }
        }

        if (closest == nullptr) {
            throw std::invalid_argument("no matching dark record");
        }

        bool drift_checked = false;
        if (before != nullptr && after != nullptr && before->record_id != after->record_id) {
            const auto before_summary = summarize_dark(*before);
            const auto after_summary = summarize_dark(*after);
            if (std::abs(before_summary.ai0_mean - after_summary.ai0_mean) > limits.max_ai0_change
                || std::abs(before_summary.ai1_mean - after_summary.ai1_mean)
                    > limits.max_ai1_change) {
                throw std::invalid_argument("dark drift exceeds limit");
            }
            drift_checked = true;
        }

        return { .selected = summarize_dark(*closest), .drift_checked = drift_checked };
    }

} // namespace

DarkSelection select_dark(
    const RawScan& scan, std::span<const DarkRecord> darks, DarkDriftLimits limits)
{
    return select_dark_in_interval(scan.session_id, scan.instrument_configuration_id,
        scan.acquired_at, scan.acquired_at, darks, limits);
}

DarkSelection select_dark_for_pair(const RawScan& sample, const RawScan& blank,
    std::span<const DarkRecord> darks, DarkDriftLimits limits)
{
    if (sample.kind != ScanKind::Sample || blank.kind != ScanKind::Blank
        || sample.session_id.empty() || sample.session_id != blank.session_id
        || sample.instrument_configuration_id.empty()
        || sample.instrument_configuration_id != blank.instrument_configuration_id) {
        throw std::invalid_argument("sample and blank do not share a dark configuration");
    }

    return select_dark_in_interval(sample.session_id, sample.instrument_configuration_id,
        std::min(sample.acquired_at, blank.acquired_at),
        std::max(sample.acquired_at, blank.acquired_at), darks, limits);
}
}
