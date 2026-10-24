#include <array>
#include <limits>
#include <stdexcept>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/workup/covariance.hpp"

namespace {

struct Fixture {
    zfit::workup::PreparedCurve curve;
    zfit::workup::DarkSummary dark;
    zfit::workup::CuvettePairCalibration cuvettes;
};

Fixture make_fixture()
{
    Fixture fixture;
    fixture.curve.session_id = "session-1";
    fixture.curve.dark_record_id = "dark-1";
    fixture.curve.cuvette_calibration_id = "cuvette-1";
    fixture.curve.scan_positions_m = { -0.001, 0.001 };
    fixture.curve.transmission = { 0.8, 0.5 };
    fixture.curve.transmission_standard_errors = { 0.1, 0.2 };
    fixture.curve.dark_sensitivities = { { { 1.0, 2.0 } }, { { 3.0, -1.0 } } };

    fixture.dark.record_id = "dark-1";
    fixture.dark.session_id = "session-1";
    fixture.dark.ai0_mean_variance = 4.0;
    fixture.dark.ai1_mean_variance = 9.0;
    fixture.dark.mean_covariance = 1.0;

    fixture.cuvettes.calibration_id = "cuvette-1";
    fixture.cuvettes.correction_factor = 2.0;
    fixture.cuvettes.standard_error = 0.2;
    return fixture;
}

} // namespace

TEST_CASE("workup builds symmetric covariance from reading, dark, and cuvette errors")
{
    const auto fixture = make_fixture();
    const auto covariance
        = zfit::workup::build_measurement_covariance(fixture.curve, fixture.dark, fixture.cuvettes);

    REQUIRE(covariance.point_count == 2);
    REQUIRE(covariance.values_row_major.size() == 4);
    // Dark terms: 44, -1, 39. Cuvette terms: .0064, .004, .0025.
    // Reading noise contributes .01 and .04 to the diagonal only.
    REQUIRE(covariance[0, 0] == Catch::Approx(44.0164));
    REQUIRE(covariance[0, 1] == Catch::Approx(-0.996));
    REQUIRE(covariance[1, 0] == covariance[0, 1]);
    REQUIRE(covariance[1, 1] == Catch::Approx(39.0425));
}

TEST_CASE("workup keeps independent reading error on the covariance diagonal")
{
    auto fixture = make_fixture();
    fixture.dark.ai0_mean_variance = 0.0;
    fixture.dark.ai1_mean_variance = 0.0;
    fixture.dark.mean_covariance = 0.0;
    fixture.cuvettes.standard_error = 0.0;
    const auto covariance
        = zfit::workup::build_measurement_covariance(fixture.curve, fixture.dark, fixture.cuvettes);

    REQUIRE(covariance[0, 0] == Catch::Approx(0.01));
    REQUIRE(covariance[1, 1] == Catch::Approx(0.04));
    REQUIRE(covariance[0, 1] == Catch::Approx(0.0).margin(1e-14));
}

TEST_CASE("workup rejects missing or invalid shared covariance sources")
{
    auto fixture = make_fixture();
    fixture.cuvettes.standard_error.reset();
    REQUIRE_THROWS_AS(
        zfit::workup::build_measurement_covariance(fixture.curve, fixture.dark, fixture.cuvettes),
        std::invalid_argument);

    fixture = make_fixture();
    fixture.curve.dark_record_id = "other-dark";
    REQUIRE_THROWS_AS(
        zfit::workup::build_measurement_covariance(fixture.curve, fixture.dark, fixture.cuvettes),
        std::invalid_argument);

    fixture = make_fixture();
    fixture.dark.mean_covariance = 7.0; // Exceeds sqrt(4 * 9).
    REQUIRE_THROWS_AS(
        zfit::workup::build_measurement_covariance(fixture.curve, fixture.dark, fixture.cuvettes),
        std::invalid_argument);

    fixture = make_fixture();
    fixture.curve.dark_sensitivities.pop_back();
    REQUIRE_THROWS_AS(
        zfit::workup::build_measurement_covariance(fixture.curve, fixture.dark, fixture.cuvettes),
        std::invalid_argument);
}

TEST_CASE("workup covariance indexing checks reads and writes")
{
    zfit::workup::MeasurementCovariance covariance {
        .point_count = 2,
        .values_row_major = std::vector<double>(4, 0.0),
    };
    covariance[0, 1] = 0.5;
    REQUIRE(covariance[0, 1] == 0.5);
    const auto& read_only = covariance;
    REQUIRE(read_only[0, 1] == 0.5);
    REQUIRE_THROWS_AS((covariance[2, 0] = 1.0), std::out_of_range);
    REQUIRE_THROWS_AS((read_only[0, 2]), std::out_of_range);

    covariance.values_row_major.pop_back();
    REQUIRE_THROWS_AS((covariance[0, 0]), std::out_of_range);
}
