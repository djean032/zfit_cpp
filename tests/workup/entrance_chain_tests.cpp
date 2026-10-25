#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

#include "zfit/workup/entrance_chain.hpp"

namespace {

using zfit::workup::compute_entrance_chain;
using zfit::workup::EntranceChainInputs;

// Fused silica window and toluene solvent at 532 nm, raw values with
// catalogue tolerances.
EntranceChainInputs make_inputs()
{
    EntranceChainInputs inputs;
    inputs.window_refractive_index = 1.4610;
    inputs.window_index_std = 0.005;
    inputs.solvent_refractive_index = 1.4969;
    inputs.solvent_index_std = 0.001;
    inputs.index_source_id = "index-catalogue";
    return inputs;
}

double fresnel_transmission(double from_index, double to_index)
{
    const double amplitude = (to_index - from_index) / (from_index + to_index);
    return 1.0 - amplitude * amplitude;
}

double chain_product(const EntranceChainInputs& inputs)
{
    return fresnel_transmission(1.0, inputs.window_refractive_index)
        * fresnel_transmission(inputs.window_refractive_index, inputs.solvent_refractive_index);
}

} // namespace

TEST_CASE("entrance chain computes Fresnel transmission of the entrance wall")
{
    const auto value = compute_entrance_chain(make_inputs());

    // One fused-silica wall into toluene: about 3.5% on the air surface and
    // 1.5e-4 on the window->solvent surface.
    REQUIRE(value.transmission == Catch::Approx(0.96474).epsilon(1e-4));
    REQUIRE(value.transmission > 0.0);
    REQUIRE(value.transmission <= 1.0);
    REQUIRE(value.index_source_id == "index-catalogue");
}

TEST_CASE("entrance chain propagates index tolerances through Fresnel")
{
    const auto value = compute_entrance_chain(make_inputs());

    // The analytic derivatives must agree with central finite differences of
    // the Fresnel product.
    auto inputs = make_inputs();
    constexpr double step = 1.0e-6;

    inputs.window_refractive_index += step;
    const double window_plus = chain_product(inputs);
    inputs.window_refractive_index -= 2.0 * step;
    const double window_minus = chain_product(inputs);
    inputs.window_refractive_index += step;
    const double kappa_dw = (window_plus - window_minus) / (2.0 * step);

    inputs.solvent_refractive_index += step;
    const double solvent_plus = chain_product(inputs);
    inputs.solvent_refractive_index -= 2.0 * step;
    const double solvent_minus = chain_product(inputs);
    const double kappa_ds = (solvent_plus - solvent_minus) / (2.0 * step);

    const double expected
        = std::sqrt(kappa_dw * kappa_dw * 0.005 * 0.005 + kappa_ds * kappa_ds * 0.001 * 0.001);
    REQUIRE(value.standard_error == Catch::Approx(expected).epsilon(1e-9));

    // For fused silica/toluene the window tolerance dominates by far.
    REQUIRE(value.standard_error == Catch::Approx(5.78e-4).epsilon(1e-2));
}

TEST_CASE("entrance chain transmission rises as the solvent index approaches the window")
{
    auto inputs = make_inputs();
    const double toluene = compute_entrance_chain(inputs).transmission;

    // Water sits farther from the window index: slightly more total loss.
    inputs.solvent_refractive_index = 1.3330;
    const double water = compute_entrance_chain(inputs).transmission;
    REQUIRE(water < toluene);
    REQUIRE(water == Catch::Approx(0.96286).epsilon(1e-4));
}

TEST_CASE("entrance chain rejects uncalibrated and malformed inputs")
{
    auto inputs = make_inputs();

    inputs.window_refractive_index = 0.0;
    REQUIRE_THROWS_AS(compute_entrance_chain(inputs), std::invalid_argument);
    inputs = make_inputs();
    inputs.solvent_refractive_index = -1.0;
    REQUIRE_THROWS_AS(compute_entrance_chain(inputs), std::invalid_argument);
    inputs = make_inputs();
    inputs.window_index_std = 0.0;
    REQUIRE_THROWS_AS(compute_entrance_chain(inputs), std::invalid_argument);
    inputs = make_inputs();
    inputs.solvent_index_std = 0.0;
    REQUIRE_THROWS_AS(compute_entrance_chain(inputs), std::invalid_argument);
    inputs = make_inputs();
    inputs.window_index_std = -0.005;
    REQUIRE_THROWS_AS(compute_entrance_chain(inputs), std::invalid_argument);
    inputs = make_inputs();
    inputs.window_refractive_index = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS(compute_entrance_chain(inputs), std::invalid_argument);
    inputs = make_inputs();
    inputs.index_source_id.clear();
    REQUIRE_THROWS_AS(compute_entrance_chain(inputs), std::invalid_argument);
}
