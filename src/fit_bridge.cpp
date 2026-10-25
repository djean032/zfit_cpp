#include "zfit/fit_bridge.hpp"

#include <cmath>
#include <stdexcept>

namespace zfit {
namespace {

    void validate_sample_calibration(const workup::SampleCalibration& calibration)
    {
        if (calibration.laser_characterization_id.empty()
            || calibration.energy_meter_calibration_id.empty()
            || calibration.entrance_chain_calibration_id.empty()
            || !std::isfinite(calibration.wavelength_m) || calibration.wavelength_m <= 0.0
            || !std::isfinite(calibration.pulse_fwhm_s) || calibration.pulse_fwhm_s <= 0.0
            || !std::isfinite(calibration.m2) || calibration.m2 <= 0.0) {
            throw std::invalid_argument("invalid sample calibration");
        }

        // The entrance chain is a passive transmission fraction: kappa in
        // (0, 1], with the standard error of kappa itself.
        if (!std::isfinite(calibration.entrance_chain_transmission)
            || calibration.entrance_chain_transmission <= 0.0
            || calibration.entrance_chain_transmission > 1.0) {
            throw std::invalid_argument("invalid sample calibration");
        }

        // Instrument-level and per-sample uncertainties: characterized values
        // arrive with a standard error, or the input is refused.
        for (double value : { calibration.waist_m, calibration.waist_standard_error_m,
                 calibration.m2_standard_error, calibration.energy_meter_standard_error_j,
                 calibration.entrance_chain_standard_error, calibration.concentration_mol_per_l,
                 calibration.concentration_standard_error_mol_per_l, calibration.path_length_m,
                 calibration.path_length_standard_error_m }) {
            if (!std::isfinite(value) || value <= 0.0) {
                throw std::invalid_argument("invalid sample calibration uncertainty");
            }
        }
    }

} // namespace

FitRequest make_fit_request(std::span<const workup::PreparedCurve> curves,
    std::span<const workup::EnergyCalibration> energies, std::span<const workup::DarkSummary> darks,
    std::span<const workup::CuvettePairCalibration> cuvettes,
    const workup::SampleCalibration& sample_calibration)
{
    const std::size_t curve_count = curves.size();
    if (curve_count == 0 || energies.size() != curve_count || darks.size() != curve_count
        || cuvettes.size() != curve_count) {
        throw std::invalid_argument("invalid fit request bridge inputs");
    }

    validate_sample_calibration(sample_calibration);

    FitRequest request;
    request.sample = SampleParameters {
        .thickness_m = sample_calibration.path_length_m,
        .concentration_mol_per_l = sample_calibration.concentration_mol_per_l,
    };
    request.curves.reserve(curve_count);
    request.calibration.pulse_energy_stds_j.reserve(curve_count);

    for (std::size_t k = 0; k < curve_count; ++k) {
        const auto& curve = curves[k];
        const auto& energy = energies[k];

        if (curve.transmission.size() != curve.scan_positions_m.size()
            || !std::isfinite(energy.mean_j) || energy.mean_j <= 0.0
            || !std::isfinite(energy.standard_error_j) || energy.standard_error_j < 0.0
            || energy.reading_count < 2 || energy.provenance_id.empty()) {
            throw std::invalid_argument("invalid fit request bridge inputs");
        }

        request.curves.push_back(EnergyCurve {
            .laser = LaserParameters {
                .wavelength_m = sample_calibration.wavelength_m,
                .waist_m = sample_calibration.waist_m,
                // Physically local energy at the sample plane: the
                // energymeter burst mean scaled by the entrance chain.
                .pulse_energy_j = sample_calibration.entrance_chain_transmission * energy.mean_j,
                .pulse_fwhm_s = sample_calibration.pulse_fwhm_s,
                .m2 = sample_calibration.m2,
            },
            .scan_positions_m = curve.scan_positions_m,
            .observed_transmission = curve.transmission,
            .measurement_covariance_row_major = { },
        });
        request.calibration.pulse_energy_stds_j.push_back(
            sample_calibration.entrance_chain_transmission * energy.standard_error_j);
    }

    const auto joint = assemble_curve_set_covariance(curves, darks, cuvettes);
    request.joint_measurement_covariance_row_major = std::move(joint.joint.values_row_major);

    request.calibration.energy_meter_std_j = sample_calibration.energy_meter_standard_error_j;
    // Entrance-chain uncertainty enters the propagation as a relative
    // standard deviation of the shared factor kappa.
    request.calibration.entrance_chain_std = sample_calibration.entrance_chain_standard_error
        / sample_calibration.entrance_chain_transmission;
    request.calibration.waist_std_m = sample_calibration.waist_standard_error_m;
    request.calibration.m2_std = sample_calibration.m2_standard_error;
    request.calibration.concentration_std
        = sample_calibration.concentration_standard_error_mol_per_l;
    request.calibration.path_length_std_m = sample_calibration.path_length_standard_error_m;

    return request;
}

} // namespace zfit
