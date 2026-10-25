#include <array>
#include <limits>
#include <stdexcept>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/workup/covariance.hpp"

namespace {

using zfit::workup::CurveSetMeasurementCovariance;
using zfit::workup::CuvettePairCalibration;
using zfit::workup::DarkSummary;
using zfit::workup::PreparedCurve;

struct CurveBlock {
    PreparedCurve curve;
    DarkSummary dark;
    CuvettePairCalibration cuvettes;
};

CurveBlock make_curve_block(const std::string& suffix, std::vector<double> transmission,
    const std::string& cuvette_id, double cuvette_standard_error)
{
    CurveBlock block;
    block.curve.session_id = "session-1";
    block.curve.dark_record_id = "dark-" + suffix;
    block.curve.cuvette_calibration_id = cuvette_id;
    block.curve.scan_positions_m = { -0.001, 0.001 };
    block.curve.transmission = std::move(transmission);
    block.curve.transmission_standard_errors = { 0.1, 0.2 };
    block.curve.dark_sensitivities = { { { 1.0, 2.0 } }, { { 3.0, -1.0 } } };

    block.dark.record_id = "dark-" + suffix;
    block.dark.session_id = "session-1";
    block.dark.ai0_mean_variance = 4.0;
    block.dark.ai1_mean_variance = 9.0;
    block.dark.mean_covariance = 1.0;

    block.cuvettes.calibration_id = cuvette_id;
    block.cuvettes.correction_factor = 2.0;
    block.cuvettes.standard_error = cuvette_standard_error;
    return block;
}

// Workup curve-set inputs as the bridge would hold them: parallel vectors,
// one entry per curve.
struct SetInputs {
    std::vector<PreparedCurve> curves;
    std::vector<DarkSummary> darks;
    std::vector<CuvettePairCalibration> cuvettes;
};

SetInputs make_set(int cuvette_count)
{
    SetInputs inputs;
    const auto first = make_curve_block("1", { 0.8, 0.5 }, "cuvette-a", 0.2);
    inputs.curves.push_back(first.curve);
    inputs.darks.push_back(first.dark);
    inputs.cuvettes.push_back(first.cuvettes);

    const bool share = cuvette_count == 1;
    const auto second = make_curve_block("2", { 0.9, 0.6 }, share ? "cuvette-a" : "cuvette-b", 0.2);
    inputs.curves.push_back(second.curve);
    inputs.darks.push_back(second.dark);
    inputs.cuvettes.push_back(second.cuvettes);
    return inputs;
}

} // namespace

TEST_CASE("workup assembles a block-diagonal curve set covariance for independent cuvettes")
{
    const auto inputs = make_set(2);
    const CurveSetMeasurementCovariance result
        = zfit::workup::assemble_curve_set_covariance(inputs.curves, inputs.darks, inputs.cuvettes);

    REQUIRE(result.curve_point_counts == std::vector<std::size_t> { 2, 2 });
    REQUIRE(result.joint.point_count == 4);
    REQUIRE(result.joint.values_row_major.size() == 16);

    // Off-curve entries vanish without a shared cuvette correction.
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            REQUIRE(result.joint[i, 2 + j] == 0.0);
            REQUIRE(result.joint[2 + j, i] == 0.0);
        }
    }

    // Diagonal blocks equal the per-curve covariance builders exactly. Dark
    // terms: 44, -1, 39; cuvette terms are relC^2 * T_i * T_j with each
    // curve's own transmissions; reading noise adds .01 and .04 on each
    // diagonal.
    REQUIRE(result.joint[0, 0] == Catch::Approx(44.0164));
    REQUIRE(result.joint[0, 1] == Catch::Approx(-0.996));
    REQUIRE(result.joint[1, 1] == Catch::Approx(39.0425));
    REQUIRE(result.joint[2, 2] == Catch::Approx(44.0 + .0081 + .01));
    REQUIRE(result.joint[2, 3] == Catch::Approx(-1.0 + .0054));
    REQUIRE(result.joint[3, 3] == Catch::Approx(39.0 + .0036 + .04));
}

TEST_CASE("workup adds shared-cuvette cross blocks to the curve set covariance")
{
    const auto inputs = make_set(1);
    const CurveSetMeasurementCovariance result
        = zfit::workup::assemble_curve_set_covariance(inputs.curves, inputs.darks, inputs.cuvettes);

    // dT/dC = -T/C: the shared 10% relative error gives
    // relC * relC * T_ki * T_lj between curves.
    const double relative_product = (0.2 / 2.0) * (0.2 / 2.0);
    REQUIRE(result.joint[0, 2] == Catch::Approx(relative_product * 0.8 * 0.9));
    REQUIRE(result.joint[0, 3] == Catch::Approx(relative_product * 0.8 * 0.6));
    REQUIRE(result.joint[1, 2] == Catch::Approx(relative_product * 0.5 * 0.9));
    REQUIRE(result.joint[1, 3] == Catch::Approx(relative_product * 0.5 * 0.6));

    // Symmetric across the whole matrix.
    for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            REQUIRE(result.joint[i, j] == result.joint[j, i]);
        }
    }
}

TEST_CASE("workup rejects contradictory or malformed curve set covariances")
{
    auto inputs = make_set(1);

    // Same id means the same characterized constant; disagreement is a
    // contradiction, not a per-curve difference.
    inputs.cuvettes[1].correction_factor = 2.1;
    REQUIRE_THROWS_AS(
        zfit::workup::assemble_curve_set_covariance(inputs.curves, inputs.darks, inputs.cuvettes),
        std::invalid_argument);

    // Curve and dark lists must align.
    inputs = make_set(1);
    inputs.darks.pop_back();
    REQUIRE_THROWS_AS(
        zfit::workup::assemble_curve_set_covariance(inputs.curves, inputs.darks, inputs.cuvettes),
        std::invalid_argument);

    // No cuvette calibration references at all.
    inputs = make_set(1);
    inputs.cuvettes.clear();
    REQUIRE_THROWS_AS(
        zfit::workup::assemble_curve_set_covariance(inputs.curves, inputs.darks, inputs.cuvettes),
        std::invalid_argument);

    // An empty curve set carries no covariance.
    REQUIRE_THROWS_AS(
        zfit::workup::assemble_curve_set_covariance(inputs.curves, inputs.darks, inputs.cuvettes),
        std::invalid_argument);
}
