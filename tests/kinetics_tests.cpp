#include <array>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "zfit/kinetics.hpp"

namespace {

constexpr zfit::ModelParameters parameters {
    .sigma_g = 0.1,
    .sigma_s = 0.2,
    .sigma_t = 0.3,
    .tau10 = 2.0,
    .tau21 = 4.0,
    .tau13 = 5.0,
    .tau30 = 10.0,
    .tau43 = 20.0,
};

constexpr std::array<double, zfit::n_levels> populations { 10.0, 20.0, 30.0, 40.0, 50.0 };

double sum(const std::array<double, zfit::n_levels>& values)
{
    double total { };
    for (double value : values) {
        total += value;
    }
    return total;
}

} // namespace

TEST_CASE("optical kinetics transfers population between optical levels")
{
    std::array<double, zfit::n_levels> rate { };

    zfit::optical_rhs(7.0, populations, rate, parameters);

    REQUIRE(rate[0] == Catch::Approx(-7.0));
    REQUIRE(rate[1] == Catch::Approx(-21.0));
    REQUIRE(rate[2] == Catch::Approx(28.0));
    REQUIRE(rate[3] == Catch::Approx(-84.0));
    REQUIRE(rate[4] == Catch::Approx(84.0));
    REQUIRE(sum(rate) == Catch::Approx(0.0));
}

TEST_CASE("relaxation kinetics conserves population")
{
    std::array<double, zfit::n_levels> rate { };

    zfit::relaxation_rhs(populations, rate, parameters);

    REQUIRE(rate[0] == Catch::Approx(14.0));
    REQUIRE(rate[1] == Catch::Approx(-6.5));
    REQUIRE(rate[2] == Catch::Approx(-7.5));
    REQUIRE(rate[3] == Catch::Approx(2.5));
    REQUIRE(rate[4] == Catch::Approx(-2.5));
    REQUIRE(sum(rate) == Catch::Approx(0.0));
}

TEST_CASE("relaxation Jacobian matches the relaxation RHS")
{
    std::array<double, zfit::n_levels * zfit::n_levels> jacobian { };
    zfit::relaxation_jacobian(jacobian, parameters);

    constexpr std::array<double, zfit::n_levels * zfit::n_levels> expected {
        0.0,
        0.5,
        0.0,
        0.1,
        0.0,
        0.0,
        -0.7,
        0.25,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.25,
        0.0,
        0.0,
        0.0,
        0.2,
        0.0,
        -0.1,
        0.05,
        0.0,
        0.0,
        0.0,
        0.0,
        -0.05,
    };

    for (std::size_t index = 0; index < jacobian.size(); ++index) {
        REQUIRE(jacobian[index] == Catch::Approx(expected[index]));
    }

    constexpr double perturbation = 1.0e-6;
    for (std::size_t column = 0; column < zfit::n_levels; ++column) {
        auto plus = populations;
        auto minus = populations;
        plus[column] += perturbation;
        minus[column] -= perturbation;

        std::array<double, zfit::n_levels> plus_rate { };
        std::array<double, zfit::n_levels> minus_rate { };
        zfit::relaxation_rhs(plus, plus_rate, parameters);
        zfit::relaxation_rhs(minus, minus_rate, parameters);

        for (std::size_t row = 0; row < zfit::n_levels; ++row) {
            const double finite_difference
                = (plus_rate[row] - minus_rate[row]) / (2.0 * perturbation);

            REQUIRE(jacobian[row * zfit::n_levels + column]
                == Catch::Approx(finite_difference).margin(1.0e-8));
        }
    }
}

TEST_CASE("full kinetics equals the optical and relaxation split")
{
    std::array<double, zfit::n_levels> optical { };
    std::array<double, zfit::n_levels> relaxation { };
    std::array<double, zfit::n_levels> full { };

    zfit::optical_rhs(7.0, populations, optical, parameters);
    zfit::relaxation_rhs(populations, relaxation, parameters);
    zfit::local_kinetics_rhs(7.0, populations, full, parameters);

    for (std::size_t level = 0; level < zfit::n_levels; ++level) {
        REQUIRE(full[level] == Catch::Approx(optical[level] + relaxation[level]));
    }
    REQUIRE(sum(full) == Catch::Approx(0.0));
}

TEST_CASE("full Jacobian matches the full RHS")
{
    constexpr double photon_flux = 7.0;
    std::array<double, zfit::n_levels * zfit::n_levels> jacobian { };
    zfit::local_kinetics_jacobian(photon_flux, jacobian, parameters);

    constexpr std::array<double, zfit::n_levels * zfit::n_levels> expected {
        -0.7,
        0.5,
        0.0,
        0.1,
        0.0,
        0.7,
        -2.1,
        0.25,
        0.0,
        0.0,
        0.0,
        1.4,
        -0.25,
        0.0,
        0.0,
        0.0,
        0.2,
        0.0,
        -2.2,
        0.05,
        0.0,
        0.0,
        0.0,
        2.1,
        -0.05,
    };

    for (std::size_t index = 0; index < jacobian.size(); ++index) {
        REQUIRE(jacobian[index] == Catch::Approx(expected[index]));
    }

    constexpr double perturbation = 1.0e-6;
    for (std::size_t column = 0; column < zfit::n_levels; ++column) {
        auto plus = populations;
        auto minus = populations;
        plus[column] += perturbation;
        minus[column] -= perturbation;

        std::array<double, zfit::n_levels> plus_rate { };
        std::array<double, zfit::n_levels> minus_rate { };
        zfit::local_kinetics_rhs(photon_flux, plus, plus_rate, parameters);
        zfit::local_kinetics_rhs(photon_flux, minus, minus_rate, parameters);

        for (std::size_t row = 0; row < zfit::n_levels; ++row) {
            const double finite_difference
                = (plus_rate[row] - minus_rate[row]) / (2.0 * perturbation);

            REQUIRE(jacobian[row * zfit::n_levels + column]
                == Catch::Approx(finite_difference).margin(1.0e-8));
        }
    }
}

TEST_CASE("absorption uses ground singlet and triplet populations")
{
    REQUIRE(zfit::absorption_coefficient(populations, parameters) == Catch::Approx(17.0));
}
