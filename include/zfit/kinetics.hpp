#pragma once

#include <cstddef>
#include <span>

#include "parameters.hpp"

namespace zfit {

inline constexpr std::size_t n_levels = 5;

void optical_rhs(double photon_flux, std::span<const double, n_levels> n,
    std::span<double, n_levels> dn, const ModelParameters& params) noexcept;

void relaxation_rhs(std::span<const double, n_levels> n, std::span<double, n_levels> dn,
    const ModelParameters& params) noexcept;

void local_kinetics_rhs(double photon_flux, std::span<const double, n_levels> n,
    std::span<double, n_levels> dn, const ModelParameters& params) noexcept;

void relaxation_jacobian(
    std::span<double, n_levels * n_levels> jacobian, const ModelParameters& p) noexcept;

void local_kinetics_jacobian(double photon_flux, std::span<double, n_levels * n_levels> jacobian,
    const ModelParameters& params) noexcept;

[[nodiscard]]
double absorption_coefficient(
    std::span<const double, n_levels> n, const ModelParameters& params) noexcept;
} // namespace zfit
