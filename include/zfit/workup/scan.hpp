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

} // namespace zfit::workup
