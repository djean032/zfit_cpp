# zfit Contributor Guide

## Current Architecture

Keep the physics, numerics, fitting, measurement workup, acquisition, and UI
separate. The current simulation and fitting path is:

```text
Laser + Kinetics -> MolModel -> CvodesIntegrator -> cminpack fit adapter

acquisition (future) -> zfit_workup (types only for now)
                     -> prepared absolute-transmission curves -> fit caller
GUI (future) -> acquisition, workup, fitting; never the other way around
```

- `zfit_core` (`Laser`, `Kinetics`, `MolModel`) has no SUNDIALS or optimizer
  dependency. `Laser` precomputes optical constants and on-axis irradiance;
  `Kinetics` implements five-level rates, analytic population Jacobians, and
  absorption. `MolModel` stores parameters, discretization, incident scan
  irradiances, and waveform support. Its whole-grid explicit/implicit routines
  remain available, but are **not** the production CVODES integration path.
- `zfit_cvodes` alone owns SUNDIALS resources. `CvodesIntegrator` solves one
  five-level spatial slice at a time with CVODES BDF and an analytic local
  Jacobian. It passes the outgoing time waveform to the next slice using
  log-intensity cubic-Hermite interpolation. It integrates transmitted energy
  by quadrature at the final slice and divides by analytic incident fluence.
- `zfit_fit` links cminpack v1.3.14. Its single-start `lmdif` fit accepts a
  runtime-selected subset of the eight kinetics parameters and **one or more
  energy curves jointly** with shared model parameters: `FitRequest::curves`
  holds one `EnergyCurve` per pulse energy (its own `LaserParameters`, scan
  grid, and observations). Residuals are flattened curve-order. Measurement
  covariance has two mutually exclusive modes: per-curve m x m blocks
  (all-or-none across the list, whitened block-diagonally with zero
  off-curve blocks), and `FitRequest::joint_measurement_covariance_row_major`
  (total count M x M with real cross-curve blocks, whitened as one matrix).
  It uses finite-difference optimizer Jacobians and log-parameter coordinates
  to keep selected cross sections and lifetimes positive. Unselected
  parameters remain fixed. Its residuals compare **absolute**, pre-calibrated
  transmission directly; there is no wing normalization, blank correction, or
  file parsing in the fit.
- Parameter uncertainties (when a measurement covariance is supplied and the
  status is 1-4): recompute the whitened residual Jacobian at the **final**
  fitted point with cminpack's `fdjac2` step rule (`sqrt(max(epsfcn, eps)) *
  |x|`, costing n extra joint curve evaluations per fit — never reuse
  `lmdif`'s possibly stale `fjac`), then Eigen `JacobiSVD` for rank and the
  log-parameter covariance, converted with
  `Cov(p_i, p_j) = p_i p_j Cov(log p_i, log p_j)`. Rank-deficient fits
  return the rank but **no** covariance/SEs; χ² is not rescaled into the
  SEs. Offline note: `epsfcn` below ~1e-9 makes the optimizer's internal
  step fall below CVODES noise and can converge to wrong parameters; keep
  the default 1e-6. **Total parameter SEs** add the linearized calibration
  propagation when `FitRequest::calibration` is supplied (any field set;
  absent = conditional SEs only): direction columns are finite differences
  of the **raw** residuals at the fitted point with respect to the model
  inputs (per-curve energy SEM from one all-curves-perturbed joint run;
  shared energymeter cal and entrance chain (relative), waist, M2,
  concentration, path length), whitened
  column-wise with the same factor, then
  `C_total_log = C_cond_log + N diag(sigma^2) N^T` with
  `N = C_cond_log (Jw^T A_w)`, scaled to physical units via `p_i p_j`. The
  fit returns `total_covariance_row_major`, `total_standard_errors`, and
  per-parameter variance contributions for **eight** calibration groups
  (0 energy SEM, 1 energymeter cal, 2 entrance chain, 3 waist, 4 M2,
  5 concentration, 6 path length, 7 fixed-parameter priors), so total
  variance = conditional variance + all contributions exactly. The
  energymeter column uses an absolute step (`h * E_k`); the entrance-chain
  column is relative (step `h`), and both come from the same
  all-curves-perturbed evaluation. Fixed parameters can carry priors in
  `CalibrationUncertainties::fixed_parameters`: a zero standard error
  claims the held-out value exact (no propagation column — used for
  placements verified insensitive, e.g. fs-class channels); a positive one
  adds one relative-step FD direction column per nonzero entry (same
  `epsfcn` caveat as the other columns, one extra joint-curve evaluation
  per fit), contributing to group 7. Priors must have unique valid
  parameter ids, finite nonnegative sigmas, and may never name a
  `fitted_parameters` entry (a fitted parameter carries its own SE; a prior
  would double-count). A calibration block may carry only priors; the
  per-curve energy list stays optional, and shared-only blocks (e.g.
  entrance chain alone) reach the propagation. `epsfcn` <
  ~1e-9 invalidates these FD direction columns too (same step rule as the
  Jacobian).
- `zfit_workup` is a headless library with no CVODES, cminpack, GUI, or
  hardware dependency, spanning raw detector-reading types through prepared
  curves, dark selection, energy statistics, and calibration covariance
  assembly: `EnergyCalibration`/`summarize_energy_measurement()`
  (multi-shot mean + SEM per curve), `compute_entrance_chain()` (the
  computed entrance-window Fresnel transmission kappa and its standard
  error, from raw material indices and catalogue tolerances),
  `SampleCalibration` (instrument-level
  laser/meter constants with provenance ids plus per-sample values, with
  missing/zero uncertainties refused by consumers), and
  `assemble_curve_set_covariance()` (per-curve blocks plus shared-cuvette
  cross blocks keyed on equal `calibration_id`; equal ids with differing
  records are rejected as contradictory). It still does not read files or
  control instruments. `build_measurement_covariance()` produces per-curve
  covariance the fit consumes; nothing crosses into physics libraries
  except the explicit bridge below.
- The workup-to-fit bridge lives in the fit layer
  (`include/zfit/fit_bridge.hpp`): `make_fit_request()` consumes workup
  values, uncertainties, and provenance only — never calibration-record
  files — and produces laser/sample inputs, the joint covariance with
  shared-C cross blocks, and the `FitRequest::calibration` block. Curve
  energies are physically local: the bridge scales the energymeter burst
  mean and its SEM by the computed entrance-chain kappa and forwards
  sigma_kappa/kappa as the relative `entrance_chain_std`. It refuses
  uncalibrated inputs (missing or zero instrument-level/sample uncertainties
  throw); the caller still chooses `initial_parameters` and
  `fitted_parameters`.
- Hardware control and GUI targets do not exist yet. Keep them optional when
  added. Do not expose SUNDIALS types to workup, fitting callers, or UI code.
  Each independent fit/simulation must own distinct CVODES state.

## Kinetics And Observable

The population order is `S0, S1, S2, T1, T2`. `tau21` is S2 -> S1 and
`tau13` is the S1 -> T1 intersystem-crossing time, as documented in
`parameters.hpp`; keep both the relaxation RHS and its Jacobian consistent.
Local kinetics uses fixed-extent spans, without a global population array.
Preserve population conservation and check significant negative populations.

Each scan position is independent. For each spatial slice, initialize
`S0 = number_density_per_cm3` and the other four levels to zero. Integrate the
entire pulse over

```cpp
t_start = options.pulse_start_fwhm * pulse_fwhm;
t_end = options.pulse_end_fwhm * pulse_fwhm;
```

Sample outgoing intensity at `waveform_knots` times and use those samples as
the next slice's incoming waveform. CVODES takes adaptive internal time steps;
waveform knots are not independent pulse solves. The observable is on-axis,
pulse-energy transmission: integrated outgoing irradiance divided by
integrated incident irradiance. Do not change that definition as part of data
calibration or optimizer work. There is no radial integration at present.

`SolverOptions` defaults to `nz = 64`, population relative tolerance `1e-6`,
quadrature relative tolerance `1e-7`, `atol = 1e-10`, 513 waveform knots, and
a pulse window of -4 to +4 FWHM. The fitting request currently defaults to
129 knots for a faster first pass; verify a final fit with 513 knots. Preserve
the transparent-sample transmission target of 1 within `1e-5`. The benchmark
compares 129, 257, and 513 knots against the 513-knot curve; its settings are
not a guarantee for every parameter range.

## Measurement And Workup Boundary

Production fits take calibrated absolute molecular transmission. Raw
transmitted/reference detector ratios can exceed one without physical sample
transmission exceeding one. A matched solvent/cuvette blank, with the same
energy, OD filters, detector configuration, and scan positions, is needed to
remove channel response while retaining the molecule's linear-absorption
baseline. For future workup, subtract dark signals and calculate sample/blank
ratios **outside the fit and physics libraries**. Do not force the wings to
one, clip measured values, or insert a gain factor into the physical solver.

Retain individual readings, settings, timestamps, calibration provenance,
sample/blank pairing, and units. Workup provides calibrated curves,
per-point measurement standard errors, multi-shot energy statistics, and
the curve-set covariance with shared-cuvette cross blocks. Record
uncertainties shared by an energy block, sample, or session separately from
per-point noise; do not silently treat them as independent readings. Beam
waist, M2, energymeter calibration, concentration, and path-length
uncertainties propagate into **total** parameter SEs through the fit's
linearized calibration layer, validated against seeded coherent Monte Carlo
refits. Remaining workup-to-fit work: coherent refit validation draws
against per-curve provenance the acquisition layer has not stored yet, and
per-curve C only for cases where the cuvettes actually differ between
energy levels.

Settled calibration conventions (record when implementing):

- Error tiers: per-point reading noise; per-curve errors (multi-shot
  energy measurement SEM via `EnergyCalibration`, within-curve dark
  projection); per-sample values (concentration, path length); and
  **instrument-level constants characterized once and entered into the
  system, not re-measured per session**: beam waist, M2, energymeter
  calibration uncertainty, and the cuvette pair correction C. The same
  cuvette pair is used at all energy levels of a sample scan, so C (and
  its uncertainty) is **shared across all of a sample's curves** and must
  appear as cross-curve covariance blocks
  (`(dT_k/dC)(dT_l/dC) sigma_C^2`); per-curve C applies only if the
  cuvettes actually differ.
- **Entrance chain (kappa): computed, never hand-entered.** The
  energymeter samples the beam in free air directly in front of the
  cuvette, so kappa — the fraction of the meter-measured energy reaching
  the sample plane — is exactly the single entrance-window Fresnel
  product `T(air -> window) * T(window -> solvent)` at normal incidence
  (polarization only enters at second order in the incidence angle;
  negligible for the focusing cone). `workup::compute_entrance_chain()`
  takes the raw window and solvent refractive indices at the laser
  wavelength plus their catalogue tolerances (both required, positive:
  claiming an exact index is refused) and returns kappa and its standard
  error from the Fresnel derivatives. There is **no characterization
  floor** and no per-session chain measurement; placement/drift effects
  ride inside the tolerated index tolerances by user choice. Material
  presets (name -> n, sigma_n) are a future GUI concern; the headless
  layer takes raw values with an index-source provenance id. Sample
  transmission beyond the sample plane — the exit wall and the sample/blank
  pairing — is handled separately by the cuvette pair correction C, not by
  kappa. Solvent linear absorption is settled negligible and is not part
  of kappa.
- A **dark measurement is taken immediately before each blank**. Dark
  selection keeps a **soft, validated** time window with drift-check
  flagging rather than hard rejection.
- Instrument-level constants each carry provenance pointing at their
  characterization record (laser, energymeter, cuvette campaign). The
  workup and fit layers consume only values, uncertainties, and
  provenance; they never read calibration-record files. Uncalibrated
  inputs must be refused (missing or zero uncertainty), not silently
  treated as exact.
- Calibration uncertainties propagated into **total** parameter SEs:
  per-curve **pulse energy** (SEM, independent across curves) plus shared
  **energymeter calibration**, **entrance chain** (kappa, relative std
  from material index tolerances; see above), **beam waist**, **M2**,
  **concentration**, and **cuvette path length** (included from the start;
  the uncertainty source for L is still to be supplied).

The legacy comparison test and standalone benchmark contain test-only wing
shifts; these are **not** preprocessing performed by `zfit_fit`.

## Validation And Next Steps

As of the fixed-parameter-prior baseline, **110 CTest cases** pass with
fitting enabled. They cover kinetics/Jacobians, laser and waveform behavior,
transparent transmission, scan reset, spatial/time/waveform convergence,
a historical legacy curve fixture, single- and two-parameter synthetic fit
recovery (including a 513-knot verification), fit diagnostics (raw norms,
converged flag), covariance/SE reporting (order symmetry, FD-step and
513-knot SE stability, rank-deficient handling), joint multi-energy fits
(shared-parameter recovery, SE reduction, per-curve weighting, and
all-or-none covariance validation), workup energy statistics, curve-set
covariance assembly (block-diagonal and shared-C cross blocks),
`make_fit_request()` bridge validation and an end-to-end synthetic
readings-to-workup-to-bridge-to-joint-fit recovery, and total-SE tests
(no-calibration absence, conditional-plus-contributions identity, a
seeded Monte-Carlo coherent-refit match of the linearized propagation,
entrance-chain-only propagation (regression: the gate used to drop it),
fixed-parameter prior propagation (8-group identity, monotone widening),
zero-prior-equals-exact, and prior validation refusals). The
historical legacy comparison
is not a gate for calibrating new instrument data: it uses its own numerical
grid and test-only wing convention. Keep the suite green when changing layers.

Priorities, in order:

1. ~~Build the headless workup/acquisition boundary and get calibrated sample
   and blank data into the absolute-transmission fit.~~ Done at the boundary
   value level (types + statistics + covariance assembly + bridge), with
   tests including raw-ratio-above-one handling in the workup.
2. ~~Propagate shared calibration errors coherently through the joint fit.~~
   Done: the bridge assembles the joint covariance with shared-C cross-curve
   blocks, and the fit computes total parameter SEs by linearized
   calibration propagation over alpha = { per-curve energy SEM, energymeter
   calibration, entrance chain (computed kappa), waist, M2, concentration,
   path length }, validated against
   seeded Monte-Carlo coherent refits. Remaining: noise-realization
   validation and provenance fields carried by a future acquisition layer
   into per-curve records, plus per-curve C (only when the cuvettes
   actually differ between energy levels).
3. Add instrument control and a GUI behind optional build targets, keeping
   acquisition, workup, and numerical physics independent.
4. **Then** investigate analytic sensitivities as a measured improvement to
   fitting. They would replace how the whitened Jacobian is obtained (same
   coordinate convention, whiten each column, unchanged SVD/SE path) and
   potentially cminpack's optimizer Jacobian; validate sensitivities against
   the finite-difference reference for the selected parameters and show the
   improvement in benchmarks. Do not revive the old ARKODE whole-grid plan
   without a new accuracy/performance justification.

Profile before changing waveform interpolation, quadrature behavior, buffer
reuse, or local kinetics. Changes must preserve transparent-sample, spatial,
waveform-grid, tolerance, and fitting acceptance tests. Multistart and
coarse-to-fine polish are future fitting improvements, not existing features.
