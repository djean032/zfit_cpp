#pragma once

#include <span>

#include "zfit/fit.hpp"
#include "zfit/workup/covariance.hpp"
#include "zfit/workup/dark.hpp"
#include "zfit/workup/energy.hpp"
#include "zfit/workup/scan.hpp"

namespace zfit {

// Workup-to-fit boundary: assembles a fully wired FitRequest from workup
// outputs — one EnergyCurve per prepared curve (laser settings from the
// calibration, energies from the multi-shot means), the joint measurement
// covariance with shared-cuvette cross-curve blocks, and the calibration
// uncertainty block for total SEs (per-curve energy SEM plus instrument-level
// shared errors).
//
// Uncalibrated inputs are refused: every instrument-level and per-sample
// uncertainty must be present and positive. Missing estimates are never
// treated as exact. The bridge consumes workup values, uncertainties, and
// provenance only; it never reads calibration-record files.
//
// The returned request has no fitted_parameters or initial_parameters yet:
// the caller chooses which parameters to fit and where to start.
FitRequest make_fit_request(std::span<const workup::PreparedCurve> curves,
    std::span<const workup::EnergyCalibration> energies, std::span<const workup::DarkSummary> darks,
    std::span<const workup::CuvettePairCalibration> cuvettes,
    const workup::SampleCalibration& sample_calibration);

} // namespace zfit
