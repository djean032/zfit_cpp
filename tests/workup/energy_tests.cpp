#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/workup/energy.hpp"

namespace {

using zfit::workup::EnergyCalibration;

zfit::workup::EnergyCalibration summary(std::vector<double> energies, const std::string& id)
{
    return zfit::workup::summarize_energy_measurement(energies, id);
}

} // namespace

TEST_CASE("workup summarizes a multi-shot energy burst")
{
    const auto burst = summary({ 1.0, 2.0 }, "meter-1");

    REQUIRE(burst.mean_j == Catch::Approx(1.5));
    // Sample variance of the mean for two readings: (0.5^2) / 1.
    REQUIRE(burst.standard_error_j == Catch::Approx(0.5));
    REQUIRE(burst.reading_count == 2);
    REQUIRE(burst.provenance_id == "meter-1");
}

TEST_CASE("workup energy SEM follows the standard error of the mean")
{
    const auto burst = summary({ 1.0, 2.0, 3.0 }, "meter-1");

    REQUIRE(burst.mean_j == Catch::Approx(2.0));
    REQUIRE(burst.standard_error_j == Catch::Approx(std::sqrt(2.0 / 6.0)));
    REQUIRE(burst.reading_count == 3);
}

TEST_CASE("workup energy SEM is zero only for repeated identical readings")
{
    const auto burst = summary({ 5.0e-9, 5.0e-9, 5.0e-9 }, "meter-2");

    REQUIRE(burst.mean_j == Catch::Approx(5.0e-9));
    REQUIRE(burst.standard_error_j == 0.0);
    REQUIRE(burst.reading_count == 3);
}

TEST_CASE("workup refuses degenerate energy measurements")
{
    REQUIRE_THROWS_AS(summary({ 1.0 }, "meter-1"), std::invalid_argument);
    REQUIRE_THROWS_AS(summary({ }, "meter-1"), std::invalid_argument);
    REQUIRE_THROWS_AS(summary({ 1.0, 0.0 }, "meter-1"), std::invalid_argument);
    REQUIRE_THROWS_AS(summary({ 1.0, -2.0 }, "meter-1"), std::invalid_argument);
    REQUIRE_THROWS_AS(summary({ 1.0, 2.0 }, ""), std::invalid_argument);

    const double nonfinite = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS(summary({ 1.0, nonfinite }, "meter-1"), std::invalid_argument);
}
