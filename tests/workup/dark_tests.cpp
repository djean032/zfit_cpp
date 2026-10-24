#include <chrono>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/workup/dark.hpp"

namespace {

zfit::workup::DarkRecord make_dark_record()
{
    zfit::workup::DarkRecord record;
    record.record_id = "dark-1";
    record.session_id = "session-1";
    record.instrument_configuration_id = "energy-1-filters-1";
    record.readings = { { .ai0 = 0.0, .ai1 = 0.0 }, { .ai0 = 2.0, .ai1 = 4.0 } };
    return record;
}

zfit::workup::DarkRecord dark_at(
    const std::string& id, int seconds, double ai0 = 1.0, double ai1 = 2.0)
{
    auto record = make_dark_record();
    record.record_id = id;
    record.acquired_at = std::chrono::system_clock::time_point { std::chrono::seconds { seconds } };
    record.readings = { { .ai0 = ai0, .ai1 = ai1 }, { .ai0 = ai0, .ai1 = ai1 } };
    return record;
}

zfit::workup::RawScan scan_at(int seconds)
{
    zfit::workup::RawScan scan;
    scan.session_id = "session-1";
    scan.instrument_configuration_id = "energy-1-filters-1";
    scan.acquired_at = std::chrono::system_clock::time_point { std::chrono::seconds { seconds } };
    return scan;
}

} // namespace

TEST_CASE("workup summarizes paired dark readings and covariance of their means")
{
    const auto record = make_dark_record();
    const auto summary = zfit::workup::summarize_dark(record);

    REQUIRE(summary.record_id == record.record_id);
    REQUIRE(summary.session_id == record.session_id);
    REQUIRE(summary.instrument_configuration_id == record.instrument_configuration_id);
    REQUIRE(summary.acquired_at == record.acquired_at);
    REQUIRE(summary.ai0_mean == Catch::Approx(1.0));
    REQUIRE(summary.ai1_mean == Catch::Approx(2.0));
    REQUIRE(summary.ai0_mean_variance == Catch::Approx(1.0));
    REQUIRE(summary.ai1_mean_variance == Catch::Approx(4.0));
    REQUIRE(summary.mean_covariance == Catch::Approx(2.0));
}

TEST_CASE("workup rejects dark records without a session ID")
{
    auto record = make_dark_record();
    record.session_id.clear();
    REQUIRE_THROWS_AS(zfit::workup::summarize_dark(record), std::invalid_argument);
}

TEST_CASE("workup rejects dark records with insufficient readings")
{
    auto record = make_dark_record();
    record.readings.pop_back();
    REQUIRE_THROWS_AS(zfit::workup::summarize_dark(record), std::invalid_argument);
}

TEST_CASE("workup rejects nonfinite dark channels")
{
    auto record = make_dark_record();
    record.readings[0].ai0 = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS(zfit::workup::summarize_dark(record), std::invalid_argument);

    record = make_dark_record();
    record.readings[1].ai1 = std::numeric_limits<double>::infinity();
    REQUIRE_THROWS_AS(zfit::workup::summarize_dark(record), std::invalid_argument);
}

TEST_CASE("workup rejects overflowing dark statistics")
{
    auto record = make_dark_record();
    record.readings[0].ai0 = -std::numeric_limits<double>::max();
    record.readings[1].ai0 = std::numeric_limits<double>::max();
    REQUIRE_THROWS_AS(zfit::workup::summarize_dark(record), std::invalid_argument);
}

TEST_CASE("workup selects the closest session-matched dark and checks drift")
{
    const auto scan = scan_at(100);
    auto wrong_session = dark_at("other-session", 99);
    wrong_session.session_id = "session-2";
    auto wrong_configuration = dark_at("other-config", 99);
    wrong_configuration.instrument_configuration_id = "energy-2-filters-1";
    const std::vector<zfit::workup::DarkRecord> darks {
        dark_at("older", 50),
        wrong_session,
        wrong_configuration,
        dark_at("before", 90),
        dark_at("after", 105),
        dark_at("later", 150),
    };

    const auto selected
        = zfit::workup::select_dark(scan, darks, { .max_ai0_change = 0.0, .max_ai1_change = 0.0 });
    REQUIRE(selected.selected.record_id == "after");
    REQUIRE(selected.selected.session_id == scan.session_id);
    REQUIRE(selected.drift_checked);
}

TEST_CASE("workup chooses the earlier dark on a time tie")
{
    const std::vector<zfit::workup::DarkRecord> darks {
        dark_at("before", 90),
        dark_at("after", 110),
    };
    const auto selected = zfit::workup::select_dark(
        scan_at(100), darks, { .max_ai0_change = 0.0, .max_ai1_change = 0.0 });
    REQUIRE(selected.selected.record_id == "before");
    REQUIRE(selected.drift_checked);
}

TEST_CASE("workup reports when dark drift could not be checked")
{
    const std::vector<zfit::workup::DarkRecord> darks { dark_at("before", 90) };
    const auto selected = zfit::workup::select_dark(
        scan_at(100), darks, { .max_ai0_change = 0.0, .max_ai1_change = 0.0 });
    REQUIRE(selected.selected.record_id == "before");
    REQUIRE_FALSE(selected.drift_checked);

    const std::vector<zfit::workup::DarkRecord> later { dark_at("after", 110) };
    const auto selected_later = zfit::workup::select_dark(
        scan_at(100), later, { .max_ai0_change = 0.0, .max_ai1_change = 0.0 });
    REQUIRE(selected_later.selected.record_id == "after");
    REQUIRE_FALSE(selected_later.drift_checked);
}

TEST_CASE("workup rejects missing matching darks and excessive drift")
{
    const auto scan = scan_at(100);
    auto wrong_session = dark_at("other-session", 90);
    wrong_session.session_id = "session-2";
    const std::vector<zfit::workup::DarkRecord> unrelated { wrong_session };
    REQUIRE_THROWS_AS(
        zfit::workup::select_dark(scan, unrelated, { 1.0, 1.0 }), std::invalid_argument);

    const std::vector<zfit::workup::DarkRecord> ai0_drift {
        dark_at("before", 90, 1.0, 2.0),
        dark_at("after", 110, 1.1, 2.0),
    };
    REQUIRE_THROWS_AS(
        zfit::workup::select_dark(scan, ai0_drift, { 0.05, 0.05 }), std::invalid_argument);

    const std::vector<zfit::workup::DarkRecord> ai1_drift {
        dark_at("before", 90, 1.0, 2.0),
        dark_at("after", 110, 1.0, 2.1),
    };
    REQUIRE_THROWS_AS(
        zfit::workup::select_dark(scan, ai1_drift, { 0.05, 0.05 }), std::invalid_argument);
}

TEST_CASE("workup rejects invalid dark drift limits")
{
    const std::vector<zfit::workup::DarkRecord> darks { dark_at("before", 90) };
    REQUIRE_THROWS_AS(
        zfit::workup::select_dark(scan_at(100), darks, { -1.0, 0.0 }), std::invalid_argument);
    REQUIRE_THROWS_AS(zfit::workup::select_dark(
                          scan_at(100), darks, { 0.0, std::numeric_limits<double>::quiet_NaN() }),
        std::invalid_argument);
}

TEST_CASE("workup selects one dark for a sample-blank pair and checks the whole interval")
{
    auto sample = scan_at(120);
    sample.kind = zfit::workup::ScanKind::Sample;
    auto blank = scan_at(100);
    blank.kind = zfit::workup::ScanKind::Blank;

    const std::vector<zfit::workup::DarkRecord> darks {
        dark_at("before", 90),
        dark_at("between", 109),
        dark_at("after", 130),
    };
    const auto selected = zfit::workup::select_dark_for_pair(sample, blank, darks, { 0.0, 0.0 });
    REQUIRE(selected.selected.record_id == "between");
    REQUIRE(selected.drift_checked);
}

TEST_CASE("workup counts distinct darks at the pair boundaries as a drift check")
{
    auto sample = scan_at(100);
    sample.kind = zfit::workup::ScanKind::Sample;
    auto blank = scan_at(120);
    blank.kind = zfit::workup::ScanKind::Blank;

    const std::vector<zfit::workup::DarkRecord> darks {
        dark_at("first", 100),
        dark_at("last", 120),
    };
    const auto selected = zfit::workup::select_dark_for_pair(sample, blank, darks, { 0.0, 0.0 });
    REQUIRE(selected.selected.record_id == "first"); // Midpoint time tie.
    REQUIRE(selected.drift_checked);
}

TEST_CASE("workup cannot check pair drift with only a dark between the scans")
{
    auto sample = scan_at(100);
    sample.kind = zfit::workup::ScanKind::Sample;
    auto blank = scan_at(120);
    blank.kind = zfit::workup::ScanKind::Blank;

    const std::vector<zfit::workup::DarkRecord> darks { dark_at("between", 110) };
    const auto selected = zfit::workup::select_dark_for_pair(sample, blank, darks, { 0.0, 0.0 });
    REQUIRE(selected.selected.record_id == "between");
    REQUIRE_FALSE(selected.drift_checked);
}

TEST_CASE("workup rejects pair-wide dark drift even with a nearby middle dark")
{
    auto sample = scan_at(100);
    sample.kind = zfit::workup::ScanKind::Sample;
    auto blank = scan_at(120);
    blank.kind = zfit::workup::ScanKind::Blank;

    const std::vector<zfit::workup::DarkRecord> darks {
        dark_at("before", 90, 1.0, 2.0),
        dark_at("between", 110, 1.0, 2.0),
        dark_at("after", 130, 1.0, 2.2),
    };
    REQUIRE_THROWS_AS(zfit::workup::select_dark_for_pair(sample, blank, darks, { 0.1, 0.1 }),
        std::invalid_argument);
}

TEST_CASE("workup rejects a sample-blank pair from different sessions")
{
    auto sample = scan_at(100);
    sample.kind = zfit::workup::ScanKind::Sample;
    auto blank = scan_at(120);
    blank.kind = zfit::workup::ScanKind::Blank;
    blank.session_id = "another-session";
    const std::vector<zfit::workup::DarkRecord> darks { dark_at("before", 90) };
    REQUIRE_THROWS_AS(zfit::workup::select_dark_for_pair(sample, blank, darks, { 0.0, 0.0 }),
        std::invalid_argument);
}

TEST_CASE("workup treats a dark at one scan's timestamp as before that scan")
{
    const std::vector<zfit::workup::DarkRecord> darks {
        dark_at("at-scan", 100),
        dark_at("later", 110),
    };
    const auto selected = zfit::workup::select_dark(scan_at(100), darks, { 0.0, 0.0 });
    REQUIRE(selected.selected.record_id == "at-scan");
    REQUIRE(selected.drift_checked);
}
