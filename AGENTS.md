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
  runtime-selected subset of the eight kinetics parameters. It uses
  finite-difference Jacobians and log-parameter coordinates to keep selected
  cross sections and lifetimes positive. Unselected parameters remain fixed.
  Its residuals compare **absolute**, pre-calibrated transmission directly;
  there is no wing normalization, blank correction, or file parsing in the fit.
- `zfit_workup` is presently a headless, header-only boundary containing raw
  detector-reading and prepared-curve types. It does not yet perform workup,
  read files, or control instruments. Add processing as a separate library
  implementation and keep it independent of CVODES, cminpack, and GUI APIs.
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
sample/blank pairing, and units. Workup will eventually provide calibrated
curves and per-point measurement standard errors. Record uncertainties shared
by an energy block, sample, or session separately from per-point noise; do not
silently treat them as independent readings. Beam waist from a known-sample
calibration, measured pulse energies, concentration, and cuvette path length
will eventually need coherent propagation through joint fits of a sample's
five energy curves. Neither joint fitting nor uncertainty propagation is
implemented yet.

The legacy comparison test and standalone benchmark contain test-only wing
shifts; these are **not** preprocessing performed by `zfit_fit`.

## Validation And Next Steps

As of the CVODES fitting baseline, the 28 CTest cases pass with fitting
enabled. They cover kinetics/Jacobians, laser and waveform behavior,
transparent transmission, scan reset, spatial/time/waveform convergence,
a historical legacy curve fixture, and one- and two-parameter synthetic fit
recovery (including a 513-knot verification). The historical legacy comparison
is not a gate for calibrating new instrument data: it uses its own numerical
grid and test-only wing convention. Keep the suite green when changing layers.

Priorities, in order:

1. Build the headless workup/acquisition boundary and get calibrated sample
   and blank data into the absolute-transmission fit. Implement workup tests
   with synthetic detector ratios, including raw ratios above one.
2. Add fit diagnostics, optional per-point standard errors, and defensible
   covariance/uncertainty reporting. Use weighted residuals only when their
   uncertainties are defined; carry shared calibration errors coherently.
   Extend fitting to jointly analyze multiple energy curves with shared
   kinetics. Validate with synthetic noise, rank checks, and 513-knot results.
3. Add instrument control and a GUI behind optional build targets, keeping
   acquisition, workup, and numerical physics independent.
4. **Then** investigate analytic sensitivities as a measured improvement to
   fitting. Finite differences remain the working gradient path until
   sensitivities are validated against them for the selected parameters and
   shown to be worthwhile in benchmarks. Do not revive the old ARKODE
   whole-grid plan without a new accuracy/performance justification.

Profile before changing waveform interpolation, quadrature behavior, buffer
reuse, or local kinetics. Changes must preserve transparent-sample, spatial,
waveform-grid, tolerance, and fitting acceptance tests. Multistart and
coarse-to-fine polish are future fitting improvements, not existing features.
