#include <limits>
#include <stdexcept>

#include <catch2/catch_test_macros.hpp>

#include "zfit/workup/matching.hpp"

namespace {

struct PairFixture {
    zfit::workup::RawScan sample;
    zfit::workup::RawScan blank;
    zfit::workup::CuvettePairCalibration cuvettes;
};

PairFixture make_pair()
{
    PairFixture pair;
    pair.sample.kind = zfit::workup::ScanKind::Sample;
    pair.sample.session_id = "session-1";
    pair.sample.sample_id = "sample-1";
    pair.sample.solvent_id = "solvent-1";
    pair.sample.cuvette_id = "sample-cuvette";
    pair.sample.acquisition_id = "sample-scan";
    pair.sample.instrument_configuration_id = "energy-1-filters-1";
    pair.sample.measured_pulse_energy_j = 50.0e-9;
    pair.sample.num_readings_per_point = 2;
    pair.sample.points = {
        { .z_m = -0.001, .readings = { { .ai0 = 1.0, .ai1 = 2.0 }, { .ai0 = 1.1, .ai1 = 2.1 } } },
        { .z_m = 0.001, .readings = { { .ai0 = 1.0, .ai1 = 1.5 }, { .ai0 = 0.9, .ai1 = 1.4 } } },
    };

    pair.blank = pair.sample;
    pair.blank.kind = zfit::workup::ScanKind::Blank;
    pair.blank.sample_id.clear();
    pair.blank.cuvette_id = "blank-cuvette";
    pair.blank.acquisition_id = "blank-scan";

    pair.cuvettes.calibration_id = "cuvette-calibration-1";
    pair.cuvettes.sample_cuvette_id = pair.sample.cuvette_id;
    pair.cuvettes.blank_cuvette_id = pair.blank.cuvette_id;
    pair.cuvettes.instrument_configuration_id = pair.sample.instrument_configuration_id;
    pair.cuvettes.correction_factor = 1.02;
    pair.cuvettes.standard_error = 0.003;
    return pair;
}

} // namespace

TEST_CASE("workup accepts a matched sample, blank, and cuvette calibration")
{
    const auto pair = make_pair();
    REQUIRE_NOTHROW(zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes));
}

TEST_CASE("workup rejects a mismatched sample and blank z grid")
{
    auto pair = make_pair();
    pair.blank.points[1].z_m += 1.0e-4;
    REQUIRE_THROWS_AS(
        zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes), std::invalid_argument);
}

TEST_CASE("workup rejects a different detector or filter configuration")
{
    auto pair = make_pair();
    pair.blank.instrument_configuration_id = "energy-1-filters-2";
    REQUIRE_THROWS_AS(
        zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes), std::invalid_argument);
}

TEST_CASE("workup rejects an unmatched cuvette")
{
    auto pair = make_pair();
    pair.blank.cuvette_id = "another-cuvette";
    REQUIRE_THROWS_AS(
        zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes), std::invalid_argument);
}

TEST_CASE("workup rejects an invalid cuvette correction")
{
    auto pair = make_pair();
    pair.cuvettes.correction_factor = 0.0;
    REQUIRE_THROWS_AS(
        zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes), std::invalid_argument);

    pair.cuvettes.correction_factor = 1.0;
    pair.cuvettes.standard_error = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS(
        zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes), std::invalid_argument);
}

TEST_CASE("workup rejects an empty scan point")
{
    auto pair = make_pair();
    pair.sample.points[0].readings.clear();
    REQUIRE_THROWS_AS(
        zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes), std::invalid_argument);
}

TEST_CASE("workup rejects a reading count inconsistent with its scan header")
{
    auto pair = make_pair();
    pair.sample.points[0].readings.pop_back();
    REQUIRE_THROWS_AS(
        zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes), std::invalid_argument);

    pair = make_pair();
    pair.blank.points[1].readings.push_back({ .ai0 = 1.0, .ai1 = 1.0 });
    REQUIRE_THROWS_AS(
        zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes), std::invalid_argument);
}

TEST_CASE("workup requires at least two readings per point")
{
    auto pair = make_pair();
    pair.sample.num_readings_per_point = 1;
    pair.sample.points[0].readings.pop_back();
    pair.sample.points[1].readings.pop_back();
    REQUIRE_THROWS_AS(
        zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes), std::invalid_argument);
}

TEST_CASE("workup allows different valid sample and blank reading counts")
{
    auto pair = make_pair();
    pair.blank.num_readings_per_point = 3;
    for (auto& point : pair.blank.points) {
        point.readings.push_back({ .ai0 = 1.0, .ai1 = 1.0 });
    }
    REQUIRE_NOTHROW(zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes));
}

TEST_CASE("workup rejects a nonfinite blank z position")
{
    auto pair = make_pair();
    pair.blank.points[0].z_m = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS(
        zfit::workup::validate_pair(pair.sample, pair.blank, pair.cuvettes), std::invalid_argument);
}
