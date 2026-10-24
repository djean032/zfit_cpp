#include "zfit/fit.hpp"

#include <cminpackcpp.hpp>

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
        case FitParameter::Count:
            break;
        }
        throw std::invalid_argument("invalid fit parameter");
    }

    void validate_request(const FitRequest& request)
    {
        const auto m = request.scan_positions_m.size();
        const auto n = request.fitted_parameters.size();

        if (n == 0 || m < n || m != request.observed_transmission.size()
            || m > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument("invalid fit data or parameter count");
        }
        const auto& covariance = request.measurement_covariance_row_major;
        if (!covariance.empty()
            && (m > std::numeric_limits<std::size_t>::max() / m || covariance.size() != m * m)) {
            throw std::invalid_argument("invalid measurement covariance size");
        }

        if (request.max_evaluations <= 0 || !std::isfinite(request.ftol) || request.ftol <= 0.0
            || !std::isfinite(request.xtol) || request.xtol <= 0.0 || !std::isfinite(request.gtol)
            || request.gtol < 0.0 || !std::isfinite(request.epsfcn) || request.epsfcn < 0.0) {
            throw std::invalid_argument("invalid fit controls");
        }

        for (std::size_t i = 0; i < m; ++i) {
            if (!std::isfinite(request.scan_positions_m[i])
                || !std::isfinite(request.observed_transmission[i])) {
                throw std::invalid_argument("nonfinite fit data");
            }
        }

        std::array<bool, fit_parameter_count> seen { };
        auto initial = request.initial_parameters;
        for (FitParameter id : request.fitted_parameters) {
            const auto index = static_cast<std::size_t>(id);
            if (index >= fit_parameter_count || seen[index]) {
                throw std::invalid_argument("invalid or duplicate fit parameter");
            }
            seen[index] = true;

            const double value = parameter_ref(initial, id);
            if (!std::isfinite(value) || value <= 0.0) {
                throw std::invalid_argument("fitted parameters must start positive");
            }
        }
    }

    ModelParameters decode_parameters(const FitRequest& request, const double* x)
    {
        auto parameters = request.initial_parameters;

        for (std::size_t j = 0; j < request.fitted_parameters.size(); ++j) {
            auto id = request.fitted_parameters[j];
            auto initial = request.initial_parameters;
            const double reference = parameter_ref(initial, id);
            const double value = reference * std::exp(x[j]);

            if (!std::isfinite(value) || value <= 0.0) {
                throw std::runtime_error("fit produced an invalid parameter");
            }
            parameter_ref(parameters, id) = value;
        }
        return parameters;
    }

    std::vector<double> simulate_curve(const FitRequest& request, const ModelParameters& parameters)
    {
        MolModel model { parameters, request.sample, Laser { request.laser },
            request.solver_options, request.scan_positions_m };
        CvodesIntegrator integrator { model };

        std::vector<double> predicted;
        predicted.reserve(request.scan_positions_m.size());
        for (std::size_t i = 0; i < request.scan_positions_m.size(); ++i) {
            const double transmission = integrator.simulate(i);
            if (!std::isfinite(transmission)) {
                throw std::runtime_error("fit produced nonfinite transmission");
            }
            predicted.push_back(transmission);
        }
        return predicted;
    }

    void evaluate_residuals(const FitRequest& request, std::span<const double> log_parameters,
        const detail::ResidualWhitener* whitener, std::span<double> output)
    {
        if (log_parameters.size() != request.fitted_parameters.size()
            || output.size() != request.scan_positions_m.size()) {
            throw std::invalid_argument("invalid residual dimensions");
        }

        const auto parameters = decode_parameters(request, log_parameters.data());
        const auto predicted = simulate_curve(request, parameters);

        for (std::size_t i = 0; i < output.size(); ++i) {
            output[i] = predicted[i] - request.observed_transmission[i];
        }
        if (whitener) {
            whitener->whiten(output, output);
        }
    }

    Eigen::MatrixXd finite_difference_whitened_jacobian(const FitRequest& request,
        std::span<const double> log_parameters, const detail::ResidualWhitener& whitener,
        std::span<const double> baseline_whitened_residual)
    {
        const std::size_t m = request.scan_positions_m.size();
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
        const ModelParameters& fitted_parameters, const Eigen::MatrixXd& jacobian)
    {
        const std::size_t n = request.fitted_parameters.size();
        if (jacobian.rows() != static_cast<Eigen::Index>(request.scan_positions_m.size())
            || jacobian.cols() != static_cast<Eigen::Index>(n)) {
            throw std::invalid_argument("invalid uncertainty Jacobian");
        }

        FitUncertainty result;
        result.parameter_order = request.fitted_parameters;

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

        auto parameters_copy = fitted_parameters;
        std::vector<double> scales;
        scales.reserve(n);
        for (FitParameter id : request.fitted_parameters) {
            scales.push_back(parameter_ref(parameters_copy, id));
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

} // namespace

FitResult fit(const FitRequest& request)
{
    validate_request(request);

    std::optional<detail::ResidualWhitener> whitener;
    if (!request.measurement_covariance_row_major.empty()) {
        whitener.emplace(request.scan_positions_m.size(), request.measurement_covariance_row_major);
    }

    const int m = static_cast<int>(request.scan_positions_m.size());
    const int n = static_cast<int>(request.fitted_parameters.size());

    std::vector<double> x(n, 0.0);
    std::vector<double> residuals(m);
    std::exception_ptr callback_error;

    auto residual_callback = [&](int callback_m, int callback_n, const double* trial_x,
                                 double* fvec, int /*iflag*/) -> int {
        try {
            if (callback_m != m || callback_n != n) {
                throw std::logic_error("unexpected fit dimensions");
            }

            evaluate_residuals(request,
                std::span<const double> { trial_x, static_cast<std::size_t>(n) },
                whitener ? &*whitener : nullptr,
                std::span<double> { fvec, static_cast<std::size_t>(m) });
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

    const auto parameters = decode_parameters(request, x.data());
    auto predicted = simulate_curve(request, parameters);

    double sum_squares = 0.0;
    std::vector<double> final_residuals;
    if (whitener) {
        final_residuals.resize(static_cast<std::size_t>(m));
    }
    for (int i = 0; i < m; ++i) {
        const double difference = predicted[i] - request.observed_transmission[i];
        sum_squares += difference * difference;
        if (whitener) {
            final_residuals[i] = difference;
        }
    }

    std::optional<double> chi_squared;
    if (whitener) {
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
        const auto jacobian = finite_difference_whitened_jacobian(
            request, x, *whitener, final_residuals);
        uncertainty = estimate_uncertainty(request, parameters, jacobian);
    }

    return FitResult {
        .parameters = parameters,
        .predicted_transmission = std::move(predicted),
        .residual_norm = std::sqrt(sum_squares),
        .chi_squared = chi_squared,
        .uncertainty = std::move(uncertainty),
        .function_evaluations = function_evaluations,
        .status = status,
    };
}
} // namespace zfit
