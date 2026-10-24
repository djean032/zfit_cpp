#include "../include/zfit/kinetics.hpp"

#include <algorithm>
#include <array>

namespace zfit {

void relaxation_rhs(std::span<const double, n_levels> n, std::span<double, n_levels> dn,
    const ModelParameters& p) noexcept
{
    const auto r4 = 1.0 / p.tau10;
    const auto r5 = 1.0 / p.tau21;
    const auto r6 = 1.0 / p.tau13;
    const auto r7 = 1.0 / p.tau30;
    const auto r8 = 1.0 / p.tau43;

    dn[0] = n[1] * r4 + n[3] * r7;

    dn[1] = -n[1] * r4 + n[2] * r5 - n[1] * r6;

    dn[2] = -n[2] * r5;

    dn[3] = n[1] * r6 - n[3] * r7 + n[4] * r8;

    dn[4] = -n[4] * r8;
}

void relaxation_jacobian(
    std::span<double, n_levels * n_levels> jacobian, const ModelParameters& p) noexcept
{
    const auto r4 = 1.0 / p.tau10;
    const auto r5 = 1.0 / p.tau21;
    const auto r6 = 1.0 / p.tau13;
    const auto r7 = 1.0 / p.tau30;
    const auto r8 = 1.0 / p.tau43;

    const std::array<double, n_levels * n_levels> values {
        0.0,
        r4,
        0.0,
        r7,
        0.0,
        0.0,
        -r4 - r6,
        r5,
        0.0,
        0.0,
        0.0,
        0.0,
        -r5,
        0.0,
        0.0,
        0.0,
        r6,
        0.0,
        -r7,
        r8,
        0.0,
        0.0,
        0.0,
        0.0,
        -r8,
    };

    std::copy(values.begin(), values.end(), jacobian.begin());
}

void local_kinetics_jacobian(double photon_flux, std::span<double, n_levels * n_levels> jacobian,
    const ModelParameters& p) noexcept
{
    relaxation_jacobian(jacobian, p);

    const double r1 = p.sigma_g * photon_flux;
    const double r2 = p.sigma_s * photon_flux;
    const double r3 = p.sigma_t * photon_flux;

    jacobian[0 * n_levels + 0] -= r1;

    jacobian[1 * n_levels + 0] += r1;
    jacobian[1 * n_levels + 1] -= r2;

    jacobian[2 * n_levels + 1] += r2;

    jacobian[3 * n_levels + 3] -= r3;

    jacobian[4 * n_levels + 3] += r3;
}

void optical_rhs(double phi, std::span<const double, n_levels> n, std::span<double, n_levels> dn,
    const ModelParameters& p) noexcept
{

    const auto r1 = p.sigma_g * phi;
    const auto r2 = p.sigma_s * phi;
    const auto r3 = p.sigma_t * phi;

    dn[0] = -r1 * n[0];

    dn[1] = r1 * n[0] - r2 * n[1];

    dn[2] = r2 * n[1];

    dn[3] = -r3 * n[3];

    dn[4] = r3 * n[3];
}

void local_kinetics_rhs(double phi, std::span<const double, n_levels> n,
    std::span<double, n_levels> dn, const ModelParameters& p) noexcept
{
    std::array<double, n_levels> optical { };
    std::array<double, n_levels> relaxation { };

    optical_rhs(phi, n, optical, p);
    relaxation_rhs(n, relaxation, p);

    for (std::size_t i = 0; i < n_levels; ++i) {
        dn[i] = optical[i] + relaxation[i];
    }
}

double absorption_coefficient(
    std::span<const double, n_levels> n, const ModelParameters& p) noexcept
{
    return p.sigma_g * n[0] + p.sigma_s * n[1] + p.sigma_t * n[3];
}

} // namespace zfit
