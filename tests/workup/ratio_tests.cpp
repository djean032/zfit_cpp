#include <limits>
#include <stdexcept>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/workup/ratio.hpp"

namespace {

zfit::workup::RawScanPoint make_point()
{
    return {
        .z_m = 0.0,
        .readings = { { .ai0 = 1.1, .ai1 = 1.2 }, { .ai0 = 2.1, .ai1 = 4.2 } },
    };
}

zfit::workup::DarkSummary make_dark()
{
    return {
        .ai0_mean = 0.1,
        .ai1_mean = 0.2,
        .ai0_mean_variance = 0.01,
        .ai1_mean_variance = 0.04,
        .mean_covariance = 0.02,
    };
}

} // namespace

TEST_CASE("workup summarizes dark-corrected ai1 over ai0 pulse ratios")
{
    // The two corrected ratios are 1 and 2. Their sample variance is 0.5;
    // the variance of their mean is 0.5 / 2 = 0.25.
    const auto result = zfit::workup::summarize_ratios(make_point(), make_dark());
    REQUIRE(result.mean == Catch::Approx(1.5));
    REQUIRE(result.mean_variance == Catch::Approx(0.25));
    REQUIRE(result.dmean_d_dark_ai0 == Catch::Approx(1.0));
    REQUIRE(result.dmean_d_dark_ai1 == Catch::Approx(-0.75));
}

TEST_CASE("workup dark sensitivities match finite differences of pulse-ratio means")
{
    constexpr double step = 1.0e-5;
    const auto point = make_point();
    const auto dark = make_dark();
    const auto result = zfit::workup::summarize_ratios(point, dark);

    auto plus = dark;
    auto minus = dark;
    plus.ai0_mean += step;
    minus.ai0_mean -= step;
    const double dmean_d_ai0 = (zfit::workup::summarize_ratios(point, plus).mean
                                   - zfit::workup::summarize_ratios(point, minus).mean)
        / (2.0 * step);
    REQUIRE(result.dmean_d_dark_ai0 == Catch::Approx(dmean_d_ai0).margin(1.0e-8));

    plus = dark;
    minus = dark;
    plus.ai1_mean += step;
    minus.ai1_mean -= step;
    const double dmean_d_ai1 = (zfit::workup::summarize_ratios(point, plus).mean
                                   - zfit::workup::summarize_ratios(point, minus).mean)
        / (2.0 * step);
    REQUIRE(result.dmean_d_dark_ai1 == Catch::Approx(dmean_d_ai1).margin(1.0e-8));
}

TEST_CASE("workup keeps shared dark uncertainty out of pulse-ratio variance")
{
    auto dark = make_dark();
    const auto baseline = zfit::workup::summarize_ratios(make_point(), dark);
    dark.ai0_mean_variance = 100.0;
    dark.ai1_mean_variance = 100.0;
    dark.mean_covariance = -50.0;

    const auto result = zfit::workup::summarize_ratios(make_point(), dark);
    REQUIRE(result.mean == baseline.mean);
    REQUIRE(result.mean_variance == baseline.mean_variance);
    REQUIRE(result.dmean_d_dark_ai0 == baseline.dmean_d_dark_ai0);
    REQUIRE(result.dmean_d_dark_ai1 == baseline.dmean_d_dark_ai1);
}

TEST_CASE("workup rejects insufficient and invalid detector readings")
{
    auto point = make_point();
    point.readings.pop_back();
    REQUIRE_THROWS_AS(zfit::workup::summarize_ratios(point, make_dark()), std::invalid_argument);

    point = make_point();
    point.readings[0].ai1 = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS(zfit::workup::summarize_ratios(point, make_dark()), std::invalid_argument);

    point = make_point();
    point.readings[0].ai0 = 0.1; // Dark-corrected reference is zero.
    REQUIRE_THROWS_AS(zfit::workup::summarize_ratios(point, make_dark()), std::invalid_argument);

    auto dark = make_dark();
    dark.ai0_mean = std::numeric_limits<double>::infinity();
    REQUIRE_THROWS_AS(zfit::workup::summarize_ratios(make_point(), dark), std::invalid_argument);
}

TEST_CASE("workup does not silently discard a negative corrected signal")
{
    auto point = make_point();
    point.readings = { { .ai0 = 1.1, .ai1 = -0.3 }, { .ai0 = 1.1, .ai1 = 0.7 } };
    const auto result = zfit::workup::summarize_ratios(point, make_dark());
    REQUIRE(result.mean == Catch::Approx(0.0).margin(1e-14));
    REQUIRE(result.mean_variance == Catch::Approx(0.25));
}
