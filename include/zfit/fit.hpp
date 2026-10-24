#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include "laser.hpp"
#include "model.hpp"
#include "parameters.hpp"

namespace zfit {

struct FitRequest {
    ModelParameters initial_parameters;
    SampleParameters sample;
    LaserParameters laser;
    SolverOptions solver_options { .waveform_knots = 129 };

    std::vector<double> scan_positions_m;
    std::vector<double> observed_transmission;
    std::vector<FitParameter> fitted_parameters;
    std::vector<double> measurement_covariance_row_major;

    int max_evaluations { 300 };
    double ftol { 1e-6 };
    double xtol { 1e-6 };
    double gtol { 1e-6 };
    double epsfcn { 1e-6 };
};

struct FitUncertainty {
    std::vector<FitParameter> parameter_order;
    std::size_t jacobian_rank { };

    std::vector<double> covariance_row_major;
    std::vector<double> standard_errors;
};

struct FitResult {
    ModelParameters parameters;
    std::vector<double> predicted_transmission;
    double residual_norm { };
    std::optional<double> chi_squared;
    std::optional<FitUncertainty> uncertainty;
    int function_evaluations { };
    int status { };
};

FitResult fit(const FitRequest& request);

} // namespace zfit
