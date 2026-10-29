#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

#include "laser.hpp"
#include "model.hpp"
#include "parameters.hpp"

namespace zfit {

// Optional calibration-uncertainty block for total (not conditional)
// parameter standard errors. Absent when every field is empty/zero; do not
// encode zero to mean "uncalibrated" — absence means conditional SEs,
// never an exact input. Fields listed here are treated as Gaussian inputs
// with the given standard deviations.
//
// Per-curve entries are independent across curves; shared entries are fully
// correlated across all curves (instrument-level characterization).
//
// Fixed-parameter entries carry a prior on a kinetics model parameter that is held
// out of fitted_parameters (e.g. a fixed lifetime whose value comes from an
// external characterization). A zero standard error claims the value exact
// (no propagation column — for placements verified insensitive, such as
// fs-class channels); a positive one is propagated like any other group,
// one extra joint-curve evaluation per nonzero entry per fit at the same
// epsfcn-limited relative step as the shared columns.
struct FixedParameterUncertainty {
    FitParameter parameter;
    double standard_error;
};

struct CalibrationUncertainties {
    // Per-curve standard errors of the measured multi-shot mean energy.
    std::vector<double> pulse_energy_stds_j;
    // Fully correlated across all curves.
    double energy_meter_std_j { };
    // Relative standard deviation of the shared entrance-chain energy
    // factor kappa (the fraction of the energymeter energy arriving at the
    // sample plane). Scales every curve's energy identically.
    double entrance_chain_std { };
    double waist_std_m { };
    double m2_std { };
    double concentration_std { };
    double path_length_std_m { };
    // Priors on held-out model parameters. Entries must not name a
    // fitted or optical parameter, must be unique, and must be finite nonnegative:
    // zero means exact.
    std::vector<FixedParameterUncertainty> fixed_parameters;
};

// One observed transmission curve at one laser configuration (in practice,
// a pulse energy). Selected waist/M2 values override the per-curve initial
// settings with one shared value during the fit.
struct EnergyCurve {
    LaserParameters laser; // Shared waist/M2 initial values when either is fitted.
    std::vector<double> scan_positions_m;
    std::vector<double> observed_transmission;
    // Optional m*m row-major measurement covariance: either every curve
    // supplies one, or (for shared-error structure) the joint covariance
    // below is supplied instead. Exactly one of the two modes may be used.
    std::vector<double> measurement_covariance_row_major;
};

struct FitRequest {
    ModelParameters initial_parameters;
    SampleParameters sample;
    SolverOptions solver_options { .waveform_knots = 129 };

    // One or more curves, jointly fit with shared selected parameters.
    // Waist and M2, when selected, are shared across curves; their initial
    // values must agree across the curves. M2 must start strictly above 1
    // when fitted (fixed M2=1 is allowed).
    std::vector<EnergyCurve> curves;
    std::vector<FitParameter> fitted_parameters;

    // Optional full joint measurement covariance: total residual count M x M
    // row-major with real cross-curve blocks. Mutually exclusive with
    // per-curve covariances.
    std::vector<double> joint_measurement_covariance_row_major;

    // Optional calibration propagation into total SEs (no effect when
    // absent). Requires a converged fit with a measurement covariance.
    CalibrationUncertainties calibration;

    int max_evaluations { 300 };
    double ftol { 1e-6 };
    double xtol { 1e-6 };
    double gtol { 1e-6 };
    double epsfcn { 1e-6 };
};

struct FitUncertainty {
    // Selected kinetics and optics together, in FitRequest::fitted_parameters order.
    std::vector<FitParameter> parameter_order;
    std::size_t jacobian_rank { };

    std::vector<double> covariance_row_major;
    std::vector<double> standard_errors;

    // Present only when calibration uncertainties were supplied and the
    // conditional covariance is full rank: linearized calibration
    // propagation added to the conditional covariance
    // (C_total = C_cond + N C_alpha N^T), same units and order.
    std::vector<double> total_covariance_row_major;
    std::vector<double> total_standard_errors;

    // Variance contributions per calibration group to each parameter's
    // total variance, in parameter_order sequence. Column groups:
    // 0 per-curve energy measurement (SEM), 1 energymeter calibration,
    // 2 entrance chain (kappa), 3 fixed waist, 4 fixed M2, 5 concentration,
    // 6 path length, 7 fixed-parameter priors.
    static constexpr std::size_t calibration_group_count = 8;
    std::vector<std::array<double, calibration_group_count>> calibration_variance_contributions;
};

// One predicted curve per request curve, in request order.
struct EnergyCurveResult {
    LaserParameters laser; // Effective laser settings used for this curve.
    std::vector<double> predicted_transmission;
};

struct FitResult {
    ModelParameters parameters;
    // Present only when the corresponding shared optical parameter was fitted.
    // The same fitted values also appear in every result curve's laser settings.
    std::optional<double> fitted_waist_m;
    std::optional<double> fitted_m2;
    std::vector<EnergyCurveResult> curves;
    double residual_norm { };

    // Raw (unwhitened) diagnostics; all reported norms stay raw.
    double initial_residual_norm { };
    double max_abs_residual { };
    bool converged { };

    std::optional<double> chi_squared;
    std::optional<FitUncertainty> uncertainty;
    // cminpack's optimization evaluation count only; the post-fit Jacobian
    // evaluations used for uncertainties are excluded.
    int function_evaluations { };
    int status { };
};

FitResult fit(const FitRequest& request);

} // namespace zfit
