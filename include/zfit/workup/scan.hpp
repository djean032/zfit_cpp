#pragma once

#include <array>
#include <chrono>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace zfit::workup {

enum class ScanKind { Blank, Sample };

// Preserve individual detector readings; calibration and correction happen
// later, outside the physics model and the optimizer.
struct DetectorReading {
    double ai0 { };
    double ai1 { };
};

struct RawScanPoint {
    double z_m { };
    std::vector<DetectorReading> readings;
};

struct CuvettePairCalibration {
    std::string calibration_id;
    std::string sample_cuvette_id;
    std::string blank_cuvette_id;
    std::string instrument_configuration_id;

    double correction_factor { 1.0 };
    std::optional<double> standard_error;
};

struct DarkRecord {
    std::string record_id;
    std::string session_id;
    std::string instrument_configuration_id;
    std::chrono::system_clock::time_point acquired_at;
    std::vector<DetectorReading> readings;
};

struct RawScan {
    ScanKind kind { };
    std::string session_id;
    std::string sample_id;
    std::string solvent_id;
    std::string cuvette_id;
    std::chrono::system_clock::time_point acquired_at;
    double measured_pulse_energy_j { };
    std::string acquisition_id;
    std::string instrument_configuration_id;
    std::size_t num_readings_per_point { };
    std::vector<RawScanPoint> points;
};

// Calibrated absolute molecular transmission and per-point reading-noise
// standard errors. Shared dark and cuvette uncertainties remain separate,
// referenced by their provenance IDs.
struct PreparedCurve {
    std::string session_id;
    std::string sample_id;
    std::string sample_acquisition_id;
    std::string blank_acquisition_id;
    std::string dark_record_id;
    std::string cuvette_calibration_id;
    std::vector<double> scan_positions_m;
    std::vector<double> transmission;
    std::vector<double> transmission_standard_errors;
    std::vector<std::array<double, 2>> dark_sensitivities;

    bool dark_drift_checked { };
    double measured_pulse_energy_j { };
};

// Calibration values the fit consumes beyond the detector readings.
// Instrument-level constants are characterized once and entered into the
// system (never re-measured per session); per-sample values come from the
// sample manifest. The workup and fit consume values, uncertainties, and
// provenance only — never calibration-record files. Consumers must refuse
// missing or zero uncertainties instead of treating uncalibrated values as
// exact.
struct SampleCalibration {
    // Laser characterization (instrument-level).
    std::string laser_characterization_id;
    double wavelength_m { };
    double waist_m { };
    double waist_standard_error_m { };
    double pulse_fwhm_s { };
    double m2 { };
    double m2_standard_error { };

    // Energymeter calibration (instrument-level), shared across all curves
    // and sessions; seen as a fully correlated multiplier uncertainty.
    std::string energy_meter_calibration_id;
    double energy_meter_standard_error_j { };

    // Entrance optical chain: the fraction of the energymeter-measured
    // energy that actually arrives at the sample plane. Never hand-entered:
    // produced by workup::compute_entrance_chain() from the window and
    // solvent refractive indices and their catalogue tolerances (the meter
    // samples in free air just before the cuvette, so kappa is the single
    // entrance-window Fresnel product). The fit consumes locally scaled
    // energies kappa * measured_energy; the standard error is the (absolute)
    // standard deviation of kappa itself. The calibration id names the
    // index source (catalogue or datasheet) backing the values.
    std::string entrance_chain_calibration_id;
    double entrance_chain_transmission { };
    double entrance_chain_standard_error { };

    // Cuvette pair correction is typed separately (CuvettePairCalibration)
    // and characterized once with the cuvette campaign.

    // Per-sample values.
    double concentration_mol_per_l { };
    double concentration_standard_error_mol_per_l { };
    double path_length_m { };
    double path_length_standard_error_m { };
};

} // namespace zfit::workup
