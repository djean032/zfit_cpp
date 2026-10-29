#include "zfit/fit.hpp"

#include <cminpackcpp.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Eigen/Core"
#include "Eigen/SVD"
#include "zfit/cvodes_integrator.hpp"
#include "zfit/model.hpp"
#include "zfit/parameters.hpp"
#include "zfit/whitening.hpp"

namespace zfit {
namespace {

    double& parameter_ref(ModelParameters& p, FitParameter id)
    {
        switch (id) {
        case FitParameter::Sigma_G:
            return p.sigma_g;
        case FitParameter::Sigma_S:
            return p.sigma_s;
        case FitParameter::Sigma_T:
            return p.sigma_t;
        case FitParameter::Tau10:
            return p.tau10;
        case FitParameter::Tau21:
            return p.tau21;
        case FitParameter::Tau13:
            return p.tau13;
        case FitParameter::Tau30:
            return p.tau30;
        case FitParameter::Tau43:
            return p.tau43;
        case FitParameter::Waist:
        case FitParameter::M2:
        case FitParameter::Count:
            break;
        }
        throw std::invalid_argument("invalid fit parameter");
    }

    const double& parameter_ref(const ModelParameters& p, FitParameter id)
    {
        switch (id) {
        case FitParameter::Sigma_G:
            return p.sigma_g;
        case FitParameter::Sigma_S:
            return p.sigma_s;
        case FitParameter::Sigma_T:
            return p.sigma_t;
        case FitParameter::Tau10:
            return p.tau10;
        case FitParameter::Tau21:
            return p.tau21;
        case FitParameter::Tau13:
            return p.tau13;
        case FitParameter::Tau30:
            return p.tau30;
        case FitParameter::Tau43:
            return p.tau43;
        case FitParameter::Waist:
        case FitParameter::M2:
        case FitParameter::Count:
            break;
        }
        throw std::invalid_argument("invalid fit parameter");
    }

    bool selected(const FitRequest& request, FitParameter id)
    {
        return std::find(request.fitted_parameters.begin(), request.fitted_parameters.end(), id)
            != request.fitted_parameters.end();
    }

    double initial_value(const FitRequest& request, FitParameter id)
    {
        switch (id) {
        case FitParameter::Waist:
            return request.curves.front().laser.waist_m;
        case FitParameter::M2:
            return request.curves.front().laser.m2;
        default:
            return parameter_ref(request.initial_parameters, id);
        }
    }

    struct DecodedParameters {
        ModelParameters model;
        std::optional<double> waist_m;
        std::optional<double> m2;
    };

    LaserParameters effective_laser(const LaserParameters& input, const DecodedParameters& decoded)
    {
        auto laser = input;
        if (decoded.waist_m) {
            laser.waist_m = *decoded.waist_m;
        }
        if (decoded.m2) {
            laser.m2 = *decoded.m2;
        }
        return laser;
    }

    // dp/dx for the optimizer coordinate x. M2 = 1 + (M2_initial - 1) exp(x).
    double physical_scale(const DecodedParameters& decoded, FitParameter id)
    {
        switch (id) {
        case FitParameter::Waist:
            return *decoded.waist_m;
        case FitParameter::M2:
            return *decoded.m2 - 1.0;
        default:
            return parameter_ref(decoded.model, id);
        }
    }

    // Flattened residual ordering: curve 0's points, then curve 1's, and so on.
    // Residual buffers, whitening, and Jacobian rows all use this ordering.
    std::size_t total_residual_count(const FitRequest& request)
    {
        std::size_t total = 0;
        for (const EnergyCurve& curve : request.curves) {
            total += curve.scan_positions_m.size();
        }
        return total;
    }

    void validate_request(const FitRequest& request)
    {
        const auto n = request.fitted_parameters.size();

        if (request.curves.empty() || n == 0) {
            throw std::invalid_argument("invalid fit data or parameter count");
        }

        bool any_covariance = false;
        bool every_covariance = true;
        std::size_t m = 0;
        for (const EnergyCurve& curve : request.curves) {
            const auto k = curve.scan_positions_m.size();
            if (k < n || k != curve.observed_transmission.size()) {
                throw std::invalid_argument("invalid fit data or parameter count");
            }

            if (m > std::numeric_limits<std::size_t>::max() - k) {
                throw std::invalid_argument("invalid fit data or parameter count");
            }
            m += k;

            const auto& covariance = curve.measurement_covariance_row_major;
            if (covariance.empty()) {
                every_covariance = false;
            } else if (k > std::numeric_limits<std::size_t>::max() / k
                || covariance.size() != k * k) {
                throw std::invalid_argument("invalid measurement covariance size");
            } else {
                any_covariance = true;
            }
        }

        // A missing block would be noiseless data, which is a different fit
        // than a weighted one; do not silently mix them. The joint
        // covariance (with cross-curve blocks) is a separate, mutually
        // exclusive mode.
        if (any_covariance != every_covariance) {
            throw std::invalid_argument(
                "every curve or no curve must supply a measurement covariance");
        }
        const bool joint_supplied = !request.joint_measurement_covariance_row_major.empty();
        if (joint_supplied) {
            if (any_covariance) {
                throw std::invalid_argument(
                    "per-curve and joint measurement covariances are mutually exclusive");
            }
            if (m > std::numeric_limits<std::size_t>::max() / m
                || request.joint_measurement_covariance_row_major.size() != m * m) {
                throw std::invalid_argument("invalid joint measurement covariance size");
            }
        }
        if (m > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument("fit residual count exceeds cminpack limits");
        }

        if (request.max_evaluations <= 0 || !std::isfinite(request.ftol) || request.ftol <= 0.0
            || !std::isfinite(request.xtol) || request.xtol <= 0.0 || !std::isfinite(request.gtol)
            || request.gtol < 0.0 || !std::isfinite(request.epsfcn) || request.epsfcn < 0.0) {
            throw std::invalid_argument("invalid fit controls");
        }

        for (const EnergyCurve& curve : request.curves) {
            for (std::size_t i = 0; i < curve.scan_positions_m.size(); ++i) {
                if (!std::isfinite(curve.scan_positions_m[i])
                    || !std::isfinite(curve.observed_transmission[i])) {
                    throw std::invalid_argument("nonfinite fit data");
                }
            }
        }

        std::array<bool, selectable_fit_parameter_count> seen { };
        for (FitParameter id : request.fitted_parameters) {
            const auto index = static_cast<std::size_t>(id);
            if (index >= selectable_fit_parameter_count || seen[index]) {
                throw std::invalid_argument("invalid or duplicate fit parameter");
            }
            seen[index] = true;

            const double value = initial_value(request, id);
            if (!std::isfinite(value) || (id == FitParameter::M2 ? value <= 1.0 : value <= 0.0)) {
                throw std::invalid_argument("fitted parameters must start positive (M2 > 1)");
            }
            if (id == FitParameter::Waist || id == FitParameter::M2) {
                for (const auto& curve : request.curves) {
                    const double current = id == FitParameter::Waist ? curve.laser.waist_m
                                                                    : curve.laser.m2;
                    if (current != value) {
                        throw std::invalid_argument("fitted optics must share initial values");
                    }
                }
            }
        }

        // Calibration-uncertainty block: absent when all empty/zero, or a
        // fully-specified finite nonnegative set. Zero never means
        // "uncalibrated"; omit the field instead. The exception is the
        // fixed-parameter prior list, where zero is meaningful: it claims
        // the held-out value exact.
        const auto& calibration = request.calibration;
        const bool calibration_absent = calibration.pulse_energy_stds_j.empty()
            && calibration.energy_meter_std_j == 0.0 && calibration.entrance_chain_std == 0.0
            && calibration.waist_std_m == 0.0 && calibration.m2_std == 0.0
            && calibration.concentration_std == 0.0 && calibration.path_length_std_m == 0.0
            && calibration.fixed_parameters.empty();

        if (!calibration_absent) {
            // Per-curve energy uncertainties are either absent entirely or
            // must match the curve exactly; the shared scalars use zero =
            // not supplied. The block may carry only fixed-parameter
            // priors.
            if (!calibration.pulse_energy_stds_j.empty()
                && calibration.pulse_energy_stds_j.size() != request.curves.size()) {
                throw std::invalid_argument(
                    "calibration energy uncertainties must match the curve count");
            }
            for (double sigma : calibration.pulse_energy_stds_j) {
                if (!std::isfinite(sigma) || sigma < 0.0) {
                    throw std::invalid_argument("invalid calibration uncertainty");
                }
            }
            for (double sigma : { calibration.energy_meter_std_j, calibration.entrance_chain_std,
                     calibration.waist_std_m, calibration.m2_std, calibration.concentration_std,
                     calibration.path_length_std_m }) {
                if (!std::isfinite(sigma) || sigma < 0.0) {
                    throw std::invalid_argument("invalid calibration uncertainty");
                }
            }

            // Fixed-parameter priors on held-out slots: valid parameter ids,
            // finite nonnegative sigmas (zero = exact), no duplicates, and
            // never a fitted parameter (a fitted parameter already carries
            // its own SE; a prior on top would double-count).
            std::array<bool, selectable_fit_parameter_count> seen_fixed { };
            for (const auto& fixed : calibration.fixed_parameters) {
                const auto index = static_cast<std::size_t>(fixed.parameter);
                if (!std::isfinite(fixed.standard_error) || fixed.standard_error < 0.0) {
                    throw std::invalid_argument("invalid fixed-parameter uncertainty");
                }
                if (index >= static_cast<std::size_t>(FitParameter::Waist)
                    || seen_fixed[index]) {
                    throw std::invalid_argument("invalid or duplicate fixed parameter");
                }
                seen_fixed[index] = true;
                if (seen[index]) {
                    throw std::invalid_argument(
                        "a fixed-parameter uncertainty may not name a fitted parameter");
                }
            }
        }
    }

    DecodedParameters decode_parameters(const FitRequest& request, const double* x)
    {
        DecodedParameters decoded { .model = request.initial_parameters };

        for (std::size_t j = 0; j < request.fitted_parameters.size(); ++j) {
            const auto id = request.fitted_parameters[j];
            const double reference = initial_value(request, id)
                - (id == FitParameter::M2 ? 1.0 : 0.0);
            const double value = reference * std::exp(x[j]);

            if (!std::isfinite(value) || value <= 0.0) {
                throw std::runtime_error("fit produced an invalid parameter");
            }
            switch (id) {
            case FitParameter::Waist:
                decoded.waist_m = value;
                break;
            case FitParameter::M2:
                decoded.m2 = 1.0 + value;
                if (!std::isfinite(*decoded.m2) || *decoded.m2 <= 1.0) {
                    throw std::runtime_error("fit produced an invalid M2");
                }
                break;
            default:
                parameter_ref(decoded.model, id) = value;
                break;
            }
        }
        return decoded;
    }

    std::vector<double> simulate_curve(const SampleParameters& sample, const LaserParameters& laser,
        const SolverOptions& options, const std::vector<double>& positions,
        const ModelParameters& parameters)
    {
        MolModel model { parameters, sample, Laser { laser }, options, positions };
        CvodesIntegrator integrator { model };

        std::vector<double> predicted;
        predicted.reserve(positions.size());
        for (std::size_t i = 0; i < positions.size(); ++i) {
            const double transmission = integrator.simulate(i);
            if (!std::isfinite(transmission)) {
                throw std::runtime_error("fit produced nonfinite transmission");
            }
            predicted.push_back(transmission);
        }
        return predicted;
    }

    void evaluate_raw_residuals(
        const FitRequest& request, std::span<const double> log_parameters, std::span<double> output)
    {
        if (log_parameters.size() != request.fitted_parameters.size()
            || output.size() != total_residual_count(request)) {
            throw std::invalid_argument("invalid residual dimensions");
        }

        const auto decoded = decode_parameters(request, log_parameters.data());

        std::size_t offset = 0;
        for (const EnergyCurve& curve : request.curves) {
            const auto predicted = simulate_curve(request.sample, effective_laser(curve.laser, decoded),
                request.solver_options, curve.scan_positions_m, decoded.model);
            for (std::size_t i = 0; i < predicted.size(); ++i) {
                output[offset + i] = predicted[i] - curve.observed_transmission[i];
            }
            offset += predicted.size();
        }
    }

    void evaluate_residuals(const FitRequest& request, std::span<const double> log_parameters,
        const detail::ResidualWhitener* whitener, std::span<double> output)
    {
        evaluate_raw_residuals(request, log_parameters, output);
        if (whitener) {
            whitener->whiten(output, output);
        }
    }

    Eigen::MatrixXd finite_difference_whitened_jacobian(const FitRequest& request,
        std::span<const double> log_parameters, const detail::ResidualWhitener& whitener,
        std::span<const double> baseline_whitened_residual)
    {
        const std::size_t m = total_residual_count(request);
        const std::size_t n = request.fitted_parameters.size();

        if (log_parameters.size() != n || baseline_whitened_residual.size() != m) {
            throw std::invalid_argument("invalid final Jacobian dimensions");
        }

        const double eps
            = std::sqrt(std::max(request.epsfcn, std::numeric_limits<double>::epsilon()));

        Eigen::MatrixXd jacobian(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(n));

        std::vector<double> perturbed_x(log_parameters.begin(), log_parameters.end());
        std::vector<double> perturbed_residual(m);

        for (std::size_t j = 0; j < n; ++j) {
            const double original = perturbed_x[j];
            double h = eps * std::abs(original);
            if (h == 0.0)
                h = eps;

            perturbed_x[j] = original + h;
            const double actual_h = perturbed_x[j] - original;
            if (!std::isfinite(actual_h) || actual_h <= 0.0) {
                throw std::runtime_error("invalid final Jacobian step");
            }

            evaluate_residuals(request, perturbed_x, &whitener, perturbed_residual);
            perturbed_x[j] = original;

            for (std::size_t i = 0; i < m; ++i) {
                const double derivative
                    = (perturbed_residual[i] - baseline_whitened_residual[i]) / actual_h;
                if (!std::isfinite(derivative)) {
                    throw std::runtime_error("nonfinite final Jacobian");
                }

                jacobian(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = derivative;
            }
        }

        return jacobian;
    }

    FitUncertainty estimate_uncertainty(const FitRequest& request,
        const DecodedParameters& decoded, const Eigen::MatrixXd& jacobian)
    {
        const std::size_t n = request.fitted_parameters.size();
        if (jacobian.rows() != static_cast<Eigen::Index>(total_residual_count(request))
            || jacobian.cols() != static_cast<Eigen::Index>(n)) {
            throw std::invalid_argument("invalid uncertainty Jacobian");
        }

        FitUncertainty result;
        result.parameter_order = request.fitted_parameters;

        // Rank uses Eigen's default machine-precision tolerance (JacobiSVD::rank);
        // the whitened Jacobian's singular values set that scale.
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(jacobian, Eigen::ComputeThinV);
        result.jacobian_rank = static_cast<std::size_t>(svd.rank());
        if (result.jacobian_rank != n) {
            return result;
        }

        const auto inverse_squared = svd.singularValues().array().inverse().square().matrix();
        if (!inverse_squared.allFinite()) {
            throw std::runtime_error("nonfinite uncertainty singular values");
        }

        const Eigen::MatrixXd covariance_log
            = svd.matrixV() * inverse_squared.asDiagonal() * svd.matrixV().transpose();
        if (!covariance_log.allFinite()) {
            throw std::runtime_error("nonfinite log-parameter covariance");
        }

        std::vector<double> scales;
        scales.reserve(n);
        for (FitParameter id : request.fitted_parameters) {
            scales.push_back(physical_scale(decoded, id));
        }

        std::vector<double> covariance(n * n);
        std::vector<double> standard_errors(n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i; j < n; ++j) {
                const double value = scales[i] * scales[j]
                    * covariance_log(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
                if (!std::isfinite(value)) {
                    throw std::runtime_error("nonfinite physical covariance");
                }
                covariance[i * n + j] = value;
                covariance[j * n + i] = value;
            }

            if (covariance[i * n + i] <= 0.0) {
                throw std::runtime_error("invalid parameter variance");
            }
            standard_errors[i] = std::sqrt(covariance[i * n + i]);
        }

        result.covariance_row_major = std::move(covariance);
        result.standard_errors = std::move(standard_errors);
        return result;
    }

    // Adds the linearized calibration propagation to a FitUncertainty that
    // already carries a full-rank conditional covariance. Sensitivity
    // directions are finite differences of the raw residuals at the fitted
    // point with respect to the model inputs themselves, then whitened
    // column-wise with the same factor used for the residuals and Jacobian.
    void apply_calibration_propagation(const FitRequest& request,
        std::span<const double> fitted_log_parameters, const DecodedParameters& decoded,
        const Eigen::MatrixXd& whitened_jacobian, const detail::ResidualWhitener& whitener,
        std::span<const double> baseline_raw_residuals, FitUncertainty& uncertainty)
    {
        const std::size_t m = total_residual_count(request);
        const std::size_t n = request.fitted_parameters.size();
        if (fitted_log_parameters.size() != n || baseline_raw_residuals.size() != m
            || whitened_jacobian.rows() != static_cast<Eigen::Index>(m)
            || whitened_jacobian.cols() != static_cast<Eigen::Index>(n)
            || uncertainty.covariance_row_major.size() != n * n
            || uncertainty.standard_errors.size() != n) {
            throw std::invalid_argument("invalid calibration propagation dimensions");
        }

        const double h
            = std::sqrt(std::max(request.epsfcn, std::numeric_limits<double>::epsilon()));
        const std::size_t curve_count = request.curves.size();
        const auto& calibration = request.calibration;

        // Direction slots. The group value indexes
        // FitUncertainty::calibration_variance_contributions columns.
        enum Group {
            EnergyMeasurements = 0,
            EnergyMeter = 1,
            EntranceChain = 2,
            Waist = 3,
            M2 = 4,
            Concentration = 5,
            PathLength = 6,
            FixedParameters = 7,
        };
        struct Slot {
            std::size_t group;
            double sigma;
        };
        std::vector<Slot> slots;
        std::vector<std::optional<std::size_t>> slot_of_curve(curve_count);
        std::optional<std::size_t> meter_slot;
        std::optional<std::size_t> entrance_chain_slot;
        std::optional<std::size_t> waist_slot;
        std::optional<std::size_t> m2_slot;
        std::optional<std::size_t> concentration_slot;
        std::optional<std::size_t> length_slot;

        // Per-curve energy columns exist only when the list is supplied
        // (full length, validated); shared-only blocks otherwise start
        // their slots at zero without them.
        const bool per_curve_energy = !calibration.pulse_energy_stds_j.empty();
        for (std::size_t k = 0; k < curve_count; ++k) {
            const double sigma = per_curve_energy ? calibration.pulse_energy_stds_j[k] : 0.0;
            if (sigma > 0.0) {
                slot_of_curve[k] = slots.size();
                slots.push_back({ EnergyMeasurements, sigma });
            }
        }
        if (calibration.energy_meter_std_j > 0.0) {
            meter_slot = slots.size();
            slots.push_back({ EnergyMeter, calibration.energy_meter_std_j });
        }
        if (calibration.entrance_chain_std > 0.0) {
            entrance_chain_slot = slots.size();
            slots.push_back({ EntranceChain, calibration.entrance_chain_std });
        }
        if (calibration.waist_std_m > 0.0 && !selected(request, FitParameter::Waist)) {
            waist_slot = slots.size();
            slots.push_back({ Waist, calibration.waist_std_m });
        }
        if (calibration.m2_std > 0.0 && !selected(request, FitParameter::M2)) {
            m2_slot = slots.size();
            slots.push_back({ M2, calibration.m2_std });
        }
        if (calibration.concentration_std > 0.0) {
            concentration_slot = slots.size();
            slots.push_back({ Concentration, calibration.concentration_std });
        }
        if (calibration.path_length_std_m > 0.0) {
            length_slot = slots.size();
            slots.push_back({ PathLength, calibration.path_length_std_m });
        }

        // Fixed-parameter priors: one relative-step direction column each.
        // Validated already: unique ids, not fitted, sigma >= 0.
        struct FixedSlot {
            FitParameter parameter;
            std::size_t slot;
        };
        std::vector<FixedSlot> fixed_slots;
        for (const auto& fixed : calibration.fixed_parameters) {
            if (fixed.standard_error <= 0.0) {
                continue;
            }
            const double value = parameter_ref(request.initial_parameters, fixed.parameter);
            if (!std::isfinite(value) || value <= 0.0) {
                throw std::runtime_error("invalid fixed parameter for calibration propagation");
            }
            fixed_slots.push_back({ fixed.parameter, slots.size() });
            slots.push_back({ FixedParameters, fixed.standard_error });
        }

        const std::size_t alpha_count = slots.size();
        if (alpha_count == 0) {
            return;
        }

        std::vector<std::size_t> offsets;
        offsets.reserve(curve_count + 1);
        offsets.push_back(0);
        for (const EnergyCurve& curve : request.curves) {
            offsets.push_back(offsets.back() + curve.scan_positions_m.size());
        }

        // Sensitivity columns, one per slot: (r(perturbed) - r(baseline)) / step.
        std::vector<std::vector<double>> columns(alpha_count);

        auto fill_block = [&](const std::optional<std::size_t>& slot, std::size_t curve,
                              const std::vector<double>& perturbed_raw, double step) {
            auto& column = columns[*slot];
            if (column.empty()) {
                column.assign(m, 0.0);
            }
            for (std::size_t i = offsets[curve]; i < offsets[curve + 1]; ++i) {
                const double derivative = (perturbed_raw[i] - baseline_raw_residuals[i]) / step;
                if (!std::isfinite(derivative)) {
                    throw std::runtime_error("nonfinite calibration sensitivity");
                }
                column[i] = derivative;
            }
        };

        auto compute_raw = [&](auto&& mutate) {
            FitRequest modified = request;
            mutate(modified);
            std::vector<double> raw(m);
            evaluate_raw_residuals(modified, fitted_log_parameters, raw);
            return raw;
        };

        // One joint evaluation with every curve's energy perturbed supplies
        // the per-curve energy columns, the shared energymeter calibration
        // column, and the entrance-chain column at once (curves are
        // independent; the last two scale the same energies differently).
        const bool need_energy_column = meter_slot.has_value() || entrance_chain_slot.has_value()
            || std::any_of(slot_of_curve.begin(), slot_of_curve.end(),
                [](const auto& slot) { return slot.has_value(); });
        if (need_energy_column) {
            const auto energy_raw = compute_raw([&h](FitRequest& modified) {
                for (auto& curve : modified.curves) {
                    curve.laser.pulse_energy_j *= 1.0 + h;
                    if (!std::isfinite(curve.laser.pulse_energy_j)
                        || curve.laser.pulse_energy_j <= 0.0) {
                        throw std::runtime_error("invalid calibrated energy perturbation");
                    }
                }
            });

            for (std::size_t k = 0; k < curve_count; ++k) {
                const double step = h * request.curves[k].laser.pulse_energy_j;
                if (!std::isfinite(step) || step <= 0.0) {
                    throw std::runtime_error("invalid calibrated energy step");
                }

                if (slot_of_curve[k]) {
                    fill_block(slot_of_curve[k], k, energy_raw, step);
                }
                if (meter_slot) {
                    // dT/d(beta) with beta a shared absolute energy bias
                    // equals the curve's own dT/dE in every block.
                    fill_block(meter_slot, k, energy_raw, step);
                }
                if (entrance_chain_slot) {
                    // The entrance chain scales every energy identically:
                    // the column is the sensitivity to the relative factor
                    // itself, so the step is h, not h * E_k.
                    fill_block(entrance_chain_slot, k, energy_raw, h);
                }
            }
        }

        if (waist_slot) {
            const auto raw = compute_raw([&h](FitRequest& modified) {
                for (auto& curve : modified.curves) {
                    curve.laser.waist_m *= 1.0 + h;
                    if (!std::isfinite(curve.laser.waist_m) || curve.laser.waist_m <= 0.0) {
                        throw std::runtime_error("invalid calibrated waist perturbation");
                    }
                }
            });
            for (std::size_t k = 0; k < curve_count; ++k) {
                fill_block(waist_slot, k, raw, h * request.curves[k].laser.waist_m);
            }
        }

        if (m2_slot) {
            const auto raw = compute_raw([&h](FitRequest& modified) {
                for (auto& curve : modified.curves) {
                    curve.laser.m2 *= 1.0 + h;
                    if (!std::isfinite(curve.laser.m2) || curve.laser.m2 <= 0.0) {
                        throw std::runtime_error("invalid calibrated M2 perturbation");
                    }
                }
            });
            for (std::size_t k = 0; k < curve_count; ++k) {
                fill_block(m2_slot, k, raw, h * request.curves[k].laser.m2);
            }
        }

        if (concentration_slot) {
            const double value = request.sample.concentration_mol_per_l;
            if (!std::isfinite(value) || value <= 0.0) {
                throw std::runtime_error("invalid calibrated concentration");
            }
            // One shared sample value: a single perturbation moves every
            // curve, so the entire column is filled in one pass.
            const auto raw = compute_raw([&h, value](FitRequest& modified) {
                modified.sample.concentration_mol_per_l = value * (1.0 + h);
            });
            std::vector<double>& column = columns[*concentration_slot];
            column.assign(m, 0.0);
            for (std::size_t i = 0; i < m; ++i) {
                const double derivative = (raw[i] - baseline_raw_residuals[i]) / (h * value);
                if (!std::isfinite(derivative)) {
                    throw std::runtime_error("nonfinite calibration sensitivity");
                }
                column[i] = derivative;
            }
        }

        if (length_slot) {
            const double value = request.sample.thickness_m;
            if (!std::isfinite(value) || value <= 0.0) {
                throw std::runtime_error("invalid calibrated path length");
            }
            const auto raw = compute_raw([&h, value](FitRequest& modified) {
                modified.sample.thickness_m = value * (1.0 + h);
            });
            std::vector<double>& column = columns[*length_slot];
            column.assign(m, 0.0);
            for (std::size_t i = 0; i < m; ++i) {
                const double derivative = (raw[i] - baseline_raw_residuals[i]) / (h * value);
                if (!std::isfinite(derivative)) {
                    throw std::runtime_error("nonfinite calibration sensitivity");
                }
                column[i] = derivative;
            }
        }

        // One joint evaluation per nonzero prior: the held-out parameter
        // moves every curve at once, so each column fills in one pass,
        // with the same relative step rule as the shared columns.
        for (const auto& fixed : fixed_slots) {
            const double value = parameter_ref(request.initial_parameters, fixed.parameter);
            const auto raw = compute_raw([&h, fixed, value](FitRequest& modified) {
                parameter_ref(modified.initial_parameters, fixed.parameter)
                    = value * (1.0 + h);
            });
            std::vector<double>& column = columns[fixed.slot];
            column.assign(m, 0.0);
            for (std::size_t i = 0; i < m; ++i) {
                const double derivative = (raw[i] - baseline_raw_residuals[i]) / (h * value);
                if (!std::isfinite(derivative)) {
                    throw std::runtime_error("nonfinite calibration sensitivity");
                }
                column[i] = derivative;
            }
        }

        // Whiten each sensitivity column with the same factor as the
        // residuals and Jacobian.
        Eigen::MatrixXd a_w(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(alpha_count));
        for (std::size_t a = 0; a < alpha_count; ++a) {
            auto& column = columns[a];
            if (column.size() != m) {
                column.resize(m, 0.0);
            }
            whitener.whiten(
                column, std::span<double>(a_w.col(static_cast<Eigen::Index>(a)).data(), m));
        }
        if (!a_w.allFinite()) {
            throw std::runtime_error("nonfinite whitened calibration sensitivities");
        }

        // Parameter covariance: conditional plus N C_alpha N^T, all in log
        // coordinates, then scaled back to physical parameters.
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(whitened_jacobian, Eigen::ComputeThinV);
        const auto inverse_squared = svd.singularValues().array().inverse().square().matrix();
        const Eigen::MatrixXd conditional_log
            = svd.matrixV() * inverse_squared.asDiagonal() * svd.matrixV().transpose();
        const Eigen::MatrixXd n_matrix = conditional_log * (whitened_jacobian.transpose() * a_w);

        Eigen::VectorXd sigma2(static_cast<Eigen::Index>(alpha_count));
        for (std::size_t a = 0; a < alpha_count; ++a) {
            sigma2[static_cast<Eigen::Index>(a)] = slots[a].sigma * slots[a].sigma;
        }
        const Eigen::MatrixXd total_log
            = conditional_log + n_matrix * sigma2.asDiagonal() * n_matrix.transpose();
        if (!n_matrix.allFinite() || !total_log.allFinite()) {
            throw std::runtime_error("nonfinite calibration propagation");
        }

        std::vector<double> scales(n);
        for (std::size_t i = 0; i < n; ++i) {
            scales[i] = physical_scale(decoded, request.fitted_parameters[i]);
        }

        // Variance contributed to each parameter per calibration group:
        // p_i^2 * sum over the group's directions of (N(i, a) sigma_a)^2,
        // with total variance = conditional variance + all contributions.
        uncertainty.calibration_variance_contributions.assign(n, { });
        for (std::size_t a = 0; a < alpha_count; ++a) {
            const std::size_t group = slots[a].group;
            for (std::size_t i = 0; i < n; ++i) {
                const double scaled
                    = n_matrix(i, static_cast<Eigen::Index>(a)) * slots[a].sigma * scales[i];
                uncertainty.calibration_variance_contributions[i][group] += scaled * scaled;
            }
        }

        std::vector<double> total_covariance(n * n);
        std::vector<double> total_standard_errors(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                total_covariance[i * n + j] = scales[i] * scales[j] * total_log(i, j);
                if (!std::isfinite(total_covariance[i * n + j])) {
                    throw std::runtime_error("nonfinite total parameter covariance");
                }
            }
            const double variance = total_covariance[i * n + i];
            if (variance < 0.0) {
                throw std::runtime_error("negative total parameter variance");
            }
            total_standard_errors[i] = std::sqrt(variance);
        }

        uncertainty.total_covariance_row_major = std::move(total_covariance);
        uncertainty.total_standard_errors = std::move(total_standard_errors);
    }

} // namespace

FitResult fit(const FitRequest& request)
{
    validate_request(request);

    const int m = static_cast<int>(total_residual_count(request));
    const int n = static_cast<int>(request.fitted_parameters.size());

    std::optional<detail::ResidualWhitener> whitener;
    if (!request.joint_measurement_covariance_row_major.empty()) {
        // Real cross-curve blocks (shared errors): whiten with the full
        // joint covariance as supplied.
        whitener.emplace(
            static_cast<std::size_t>(m), request.joint_measurement_covariance_row_major);
    } else if (!request.curves.front().measurement_covariance_row_major.empty()) {
        // All-or-none weighted curves (validated above): assemble one
        // block-diagonal covariance with zero off-curve blocks. SPD blocks
        // make the block-diagonal matrix SPD, so the whitener needs no
        // new code path.
        std::vector<double> blocks(static_cast<std::size_t>(m) * static_cast<std::size_t>(m), 0.0);
        std::size_t offset = 0;
        for (const EnergyCurve& curve : request.curves) {
            const auto k = curve.scan_positions_m.size();
            const auto& covariance = curve.measurement_covariance_row_major;
            for (std::size_t i = 0; i < k; ++i) {
                for (std::size_t j = 0; j < k; ++j) {
                    blocks[(offset + i) * static_cast<std::size_t>(m) + (offset + j)]
                        = covariance[i * k + j];
                }
            }
            offset += k;
        }
        whitener.emplace(static_cast<std::size_t>(m), blocks);
    }

    std::vector<double> x(n, 0.0);
    std::vector<double> residuals(m);
    std::vector<double> raw_residuals(m);
    std::exception_ptr callback_error;

    bool captured_initial_norm = false;
    double initial_norm_squared = 0.0;

    auto residual_callback = [&](int callback_m, int callback_n, const double* trial_x,
                                 double* fvec, int /*iflag*/) -> int {
        try {
            if (callback_m != m || callback_n != n) {
                throw std::logic_error("unexpected fit dimensions");
            }

            const std::span<double> fvec_span { fvec, static_cast<std::size_t>(m) };
            if (whitener) {
                evaluate_raw_residuals(request,
                    std::span<const double> { trial_x, static_cast<std::size_t>(n) },
                    raw_residuals);
                whitener->whiten(raw_residuals, fvec_span);
            } else {
                evaluate_residuals(request,
                    std::span<const double> { trial_x, static_cast<std::size_t>(n) }, nullptr,
                    fvec_span);
            }

            if (!captured_initial_norm) {
                // lmdif's first evaluation is exactly at the initial parameters.
                captured_initial_norm = true;
                const std::span<const double> initial = whitener
                    ? std::span<const double>(raw_residuals)
                    : std::span<const double>(fvec, static_cast<std::size_t>(m));
                for (double value : initial) {
                    initial_norm_squared += value * value;
                }
                if (!std::isfinite(initial_norm_squared)) {
                    initial_norm_squared = 0.0;
                }
            }

            return 0;
        } catch (...) {
            callback_error = std::current_exception();
            return -1;
        }
    };

    std::vector<double> diag(n, 1.0);
    std::vector<double> fjac(static_cast<std::size_t>(m) * n);
    std::vector<int> ipvt(n);
    std::vector<double> qtf(n);
    std::vector<double> wa1(n), wa2(n), wa3(n), wa4(m);
    int function_evaluations = 0;

    const int status = cminpack::lmdif(residual_callback, m, n, x.data(), residuals.data(),
        request.ftol, request.xtol, request.gtol, request.max_evaluations, request.epsfcn,
        diag.data(), 1, 1.0, 0, &function_evaluations, fjac.data(), m, ipvt.data(), qtf.data(),
        wa1.data(), wa2.data(), wa3.data(), wa4.data());

    if (callback_error) {
        std::rethrow_exception(callback_error);
    }

    if (status <= 0) {
        throw std::runtime_error("cminpack fit failed");
    }

    const auto decoded = decode_parameters(request, x.data());

    double sum_squares = 0.0;
    double max_abs = 0.0;
    std::vector<double> final_residuals(static_cast<std::size_t>(m), 0.0);
    std::vector<EnergyCurveResult> curve_results;
    curve_results.reserve(request.curves.size());

    std::size_t offset = 0;
    for (const EnergyCurve& curve : request.curves) {
        const auto laser = effective_laser(curve.laser, decoded);
        auto predicted = simulate_curve(request.sample, laser, request.solver_options,
            curve.scan_positions_m, decoded.model);

        for (std::size_t i = 0; i < predicted.size(); ++i) {
            const double difference = predicted[i] - curve.observed_transmission[i];
            sum_squares += difference * difference;
            if (max_abs < std::abs(difference)) {
                max_abs = std::abs(difference);
            }
            final_residuals[offset + i] = difference;
        }

        offset += predicted.size();
        curve_results.push_back({ laser, std::move(predicted) });
    }

    std::optional<double> chi_squared;
    std::vector<double> raw_final_residuals;
    if (whitener) {
        // The whitener is applied in place; keep the raw residuals for the
        // calibration sensitivity baselines afterwards.
        raw_final_residuals = final_residuals;
        whitener->whiten(final_residuals, final_residuals);

        double value = 0.0;
        for (double residual : final_residuals) {
            value += residual * residual;
        }

        if (!std::isfinite(value)) {
            throw std::runtime_error("nonfinite fit chi-squared");
        }

        chi_squared = value;
    }

    std::optional<FitUncertainty> uncertainty;
    if (whitener && status >= 1 && status <= 4) {
        const auto jacobian
            = finite_difference_whitened_jacobian(request, x, *whitener, final_residuals);
        uncertainty = estimate_uncertainty(request, decoded, jacobian);

        // Total SEs: linearized calibration propagation on top of a
        // full-rank conditional covariance. Absent otherwise — the
        // conditional SEs already stand on their own.
        const auto& calibration = request.calibration;
        const bool calibration_supplied = !calibration.pulse_energy_stds_j.empty()
            || calibration.energy_meter_std_j > 0.0 || calibration.entrance_chain_std > 0.0
            || calibration.waist_std_m > 0.0 || calibration.m2_std > 0.0
            || calibration.concentration_std > 0.0 || calibration.path_length_std_m > 0.0
            || !calibration.fixed_parameters.empty();
        if (calibration_supplied && !uncertainty->covariance_row_major.empty()) {
            apply_calibration_propagation(
                request, x, decoded, jacobian, *whitener, raw_final_residuals, *uncertainty);
        }
    }

    return FitResult {
        .parameters = decoded.model,
        .fitted_waist_m = decoded.waist_m,
        .fitted_m2 = decoded.m2,
        .curves = std::move(curve_results),
        .residual_norm = std::sqrt(sum_squares),
        .initial_residual_norm = std::sqrt(initial_norm_squared),
        .max_abs_residual = max_abs,
        .converged = (status >= 1 && status <= 4),
        .chi_squared = chi_squared,
        .uncertainty = std::move(uncertainty),
        .function_evaluations = function_evaluations,
        .status = status,
    };
}
} // namespace zfit
