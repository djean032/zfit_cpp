#include "../include/zfit/cvodes_integrator.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <sundials/sundials_errors.h>
#include <sundials/sundials_linearsolver.h>
#include <sundials/sundials_matrix.h>
#include <sundials/sundials_nvector.h>
#include <sundials/sundials_types.h>
#include <type_traits>

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_ls.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_context.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

namespace zfit {

static_assert(std::is_same_v<sunrealtype, double>);

namespace {

    void check_cvodes(int status, const char* operation)
    {
        if (status != CV_SUCCESS) {
            throw std::runtime_error(operation);
        }
    };

    void check_cvodes_evolve(int status, const char* operation)
    {
        if (status != CV_SUCCESS && status != CV_TSTOP_RETURN) {
            throw std::runtime_error(operation);
        }
    };

    void check_sundials(SUNErrCode status, const char* operation)
    {
        if (status != SUN_SUCCESS) {
            throw std::runtime_error(operation);
        }
    };

    int cvodes_rhs(sunrealtype time, N_Vector population, N_Vector population_rate, void* user_data)
    {
        const auto* context = static_cast<const CvodesContext*>(user_data);
        if (context == nullptr || context->model == nullptr || context->incoming_waveform == nullptr
            || !std::isfinite(static_cast<double>(time))) {
            return 1;
        }

        if (N_VGetLength(population) != static_cast<sunindextype>(n_levels)
            || N_VGetLength(population_rate) != static_cast<sunindextype>(n_levels)) {
            return 1;
        }

        const auto* population_data = N_VGetArrayPointer(population);
        auto* population_rate_data = N_VGetArrayPointer(population_rate);
        if (population_data == nullptr || population_rate_data == nullptr) {
            return 1;
        }

        const double intensity = context->incoming_waveform->intensity(static_cast<double>(time));
        if (!(intensity >= 0.0) || !std::isfinite(intensity)) {
            return 1;
        }
        const std::span<const double, n_levels> n { population_data, n_levels };
        const std::span<double, n_levels> dn { population_rate_data, n_levels };

        local_kinetics_rhs(intensity * context->model->inverse_photon_energy_per_j, n, dn,
            context->model->model_params);
        return 0;
    }

    int cvodes_quadrature_rhs(
        sunrealtype time, N_Vector population, N_Vector transmitted_energy_rate, void* user_data)
    {
        const auto* context = static_cast<const CvodesContext*>(user_data);
        if (context == nullptr || transmitted_energy_rate == nullptr
            || N_VGetLength(transmitted_energy_rate) != 1) {
            return 1;
        }

        auto* energy_rate_data = N_VGetArrayPointer(transmitted_energy_rate);
        if (energy_rate_data == nullptr) {
            return 1;
        }

        if (!context->accumulate_transmitted_energy) {
            energy_rate_data[0] = 0.0;
            return 0;
        }

        if (context->model == nullptr || context->incoming_waveform == nullptr
            || population == nullptr
            || N_VGetLength(population) != static_cast<sunindextype>(n_levels)
            || !std::isfinite(static_cast<double>(time))) {
            return 1;
        }

        const auto* population_data = N_VGetArrayPointer(population);
        if (population_data == nullptr) {
            return 1;
        }

        const double incoming_intensity
            = context->incoming_waveform->intensity(static_cast<double>(time));
        if (!(incoming_intensity >= 0.0) || !std::isfinite(incoming_intensity)) {
            return 1;
        }

        const std::span<const double, n_levels> n { population_data, n_levels };
        const double alpha = absorption_coefficient(n, context->model->model_params);
        const double outgoing_intensity
            = incoming_intensity * std::exp(-alpha * context->model->dz_cm);

        if (!(outgoing_intensity >= 0.0) || !std::isfinite(outgoing_intensity)) {
            return 1;
        }

        energy_rate_data[0] = outgoing_intensity;
        return 0;
    }

    int cvodes_jacobian(sunrealtype time, N_Vector, N_Vector, SUNMatrix jacobian, void* user_data,
        N_Vector, N_Vector, N_Vector)
    {
        const auto* context = static_cast<const CvodesContext*>(user_data);
        if (jacobian == nullptr || context == nullptr || context->model == nullptr
            || context->incoming_waveform == nullptr || !std::isfinite(static_cast<double>(time))
            || SUNMatGetID(jacobian) != SUNMATRIX_DENSE) {
            return 1;
        }

        const double intensity = context->incoming_waveform->intensity(static_cast<double>(time));
        if (!(intensity >= 0.0) || !std::isfinite(intensity)) {
            return 1;
        }

        std::array<double, n_levels * n_levels> local { };
        local_kinetics_jacobian(intensity * context->model->inverse_photon_energy_per_j, local,
            context->model->model_params);

        if (SUNMatZero(jacobian) != SUN_SUCCESS) {
            return 1;
        }

        for (std::size_t row = 0; row < n_levels; ++row) {
            for (std::size_t col = 0; col < n_levels; ++col) {
                SM_ELEMENT_D(jacobian, static_cast<sunindextype>(row),
                    static_cast<sunindextype>(col)) = local[row * n_levels + col];
            }
        }

        return 0;
    }

} // namespace

struct CvodesIntegrator::Impl {
    MolModel* model { };
    CvodesContext context { };

    SUNContext sun_context { };
    N_Vector population { };
    N_Vector transmitted_energy { };
    void* cvodes { };
    SUNMatrix matrix { };
    SUNLinearSolver linear_solver { };

    ~Impl();
};

CvodesIntegrator::Impl::~Impl()
{
    if (cvodes != nullptr) {
        CVodeFree(&cvodes);
    }

    if (linear_solver != nullptr) {
        SUNLinSolFree(linear_solver);
    }

    if (matrix != nullptr) {
        SUNMatDestroy(matrix);
    }

    if (population != nullptr) {
        N_VDestroy(population);
    }

    if (transmitted_energy != nullptr) {
        N_VDestroy(transmitted_energy);
    }

    if (sun_context != nullptr) {
        SUNContext_Free(&sun_context);
    }
}

CvodesIntegrator::CvodesIntegrator(MolModel& model)
    : impl(std::make_unique<Impl>())
{

    const double t_start = model.solver_options.pulse_start_fwhm * model.laser.params.pulse_fwhm_s;
    const double t_end = model.solver_options.pulse_end_fwhm * model.laser.params.pulse_fwhm_s;

    if (!(t_end > t_start) || !std::isfinite(t_start) || !std::isfinite(t_end)) {
        throw std::invalid_argument("pulse end time must follow pulse start time");
    }

    if (!(model.solver_options.population_rtol > 0.0)
        || !(model.solver_options.quadrature_rtol > 0.0) || !(model.solver_options.atol > 0.0)
        || !std::isfinite(model.solver_options.population_rtol)
        || !std::isfinite(model.solver_options.quadrature_rtol)
        || !std::isfinite(model.solver_options.atol)) {
        throw std::invalid_argument("CVODES tolerances must be positive and finite");
    }

    if (model.solver_options.max_steps
        > static_cast<std::size_t>(std::numeric_limits<long int>::max())) {
        throw std::length_error("CVODES maximum step count exceeds long int");
    }

    const double population_atol = model.solver_options.atol * model.number_density_per_cm3;
    if (!(population_atol > 0.0) || !std::isfinite(population_atol)) {
        throw std::invalid_argument(
            "CVODES scaled population tolerance must be positive and finite");
    }

    impl->model = &model;
    impl->context.model = &model;

    check_sundials(
        SUNContext_Create(SUN_COMM_NULL, &impl->sun_context), "SUNContext_Create failed");

    impl->population = N_VNew_Serial(static_cast<sunindextype>(n_levels), impl->sun_context);
    if (impl->population == nullptr) {
        throw std::runtime_error("N_VNew_Serial failed");
    }

    impl->transmitted_energy = N_VNew_Serial(1, impl->sun_context);
    if (impl->transmitted_energy == nullptr) {
        throw std::runtime_error("N_VNew_Serial failed");
    }

    N_VConst(0.0, impl->transmitted_energy);

    auto* population_data = N_VGetArrayPointer(impl->population);
    if (population_data == nullptr) {
        throw std::runtime_error("N_VGetArrayPointer failed");
    }

    std::span<double, n_levels> population { population_data, n_levels };
    model.initialize_cell_state(population);

    impl->matrix = SUNDenseMatrix(static_cast<sunindextype>(n_levels),
        static_cast<sunindextype>(n_levels), impl->sun_context);
    if (impl->matrix == nullptr) {
        throw std::runtime_error("SUNDenseMatrix failed");
    }

    impl->linear_solver = SUNLinSol_Dense(impl->population, impl->matrix, impl->sun_context);
    if (impl->linear_solver == nullptr) {
        throw std::runtime_error("SUNLinSol_Dense failed");
    }

    impl->cvodes = CVodeCreate(CV_BDF, impl->sun_context);
    if (impl->cvodes == nullptr) {
        throw std::runtime_error("CVodeCreate failed");
    }

    check_cvodes(
        CVodeInit(impl->cvodes, cvodes_rhs, static_cast<sunrealtype>(t_start), impl->population),
        "CVodeInit failed");

    check_cvodes(CVodeSetUserData(impl->cvodes, &impl->context), "CVodeSetUserData failed");

    check_cvodes(CVodeQuadInit(impl->cvodes, cvodes_quadrature_rhs, impl->transmitted_energy),
        "CVodeQuadInit failed");

    check_cvodes(CVodeSetQuadErrCon(impl->cvodes, SUNFALSE), "CVodeSetQuadErrCon failed");

    check_cvodes(
        CVodeSStolerances(impl->cvodes, model.solver_options.population_rtol, population_atol),
        "CVodeSStolerances failed");

    check_cvodes(
        CVodeSetMaxNumSteps(impl->cvodes, static_cast<long int>(model.solver_options.max_steps)),
        "CVodeSetMaxNumSteps failed");

    check_cvodes(CVodeSetLinearSolver(impl->cvodes, impl->linear_solver, impl->matrix),
        "CVodeSetLinearSolver failed");

    check_cvodes(CVodeSetJacFn(impl->cvodes, cvodes_jacobian), "CVodeSetJacFn failed");

    check_cvodes(
        CVodeSetStopTime(impl->cvodes, static_cast<sunrealtype>(t_end)), "CVodeSetStopTime failed");
}

CvodesIntegrator::~CvodesIntegrator() = default;

double CvodesIntegrator::simulate(std::size_t scan_index)
{
    MolModel& model = *impl->model;

    if (model.solver_options.nz == 0) {
        throw std::invalid_argument("CVODES requires at least one spatial slice");
    }

    PulseWaveform incoming = model.incident_waveform(scan_index);
    if (!incoming.valid()) {
        throw std::logic_error("invalid incident waveform");
    }

    const double t_start = incoming.times.front();
    const double t_end = incoming.times.back();
    const double temporal_decay = model.laser.temporal_decay_per_s2;
    if (!(temporal_decay > 0.0) || !std::isfinite(temporal_decay)) {
        throw std::invalid_argument("laser temporal decay must be positive");
    }

    const double root_decay = std::sqrt(temporal_decay);
    const double temporal_integral = std::sqrt(std::numbers::pi) / (2.0 * root_decay)
        * (std::erf(root_decay * t_end) - std::erf(root_decay * t_start));
    const double incident_fluence = model.incident_peak_w_per_cm2[scan_index] * temporal_integral;
    if (!(incident_fluence > 0.0) || !std::isfinite(incident_fluence)) {
        throw std::invalid_argument("incident fluence must be positive and finite");
    }

    const double quadrature_atol = model.solver_options.atol * incident_fluence;
    if (!(quadrature_atol > 0.0) || !std::isfinite(quadrature_atol)) {
        throw std::invalid_argument("CVODES quadrature tolerance must be positive and finite");
    }

    check_cvodes(
        CVodeQuadSStolerances(impl->cvodes, model.solver_options.quadrature_rtol, quadrature_atol),
        "CVodeQuadSStolerances failed");

    auto* population_data = N_VGetArrayPointer(impl->population);
    auto* energy_data = N_VGetArrayPointer(impl->transmitted_energy);
    if (population_data == nullptr || energy_data == nullptr) {
        throw std::runtime_error("N_VGetArrayPointer failed");
    }

    std::span<double, n_levels> population { population_data, n_levels };

    for (std::size_t slice = 0; slice < model.solver_options.nz; ++slice) {
        const bool final_slice = slice + 1 == model.solver_options.nz;

        if (!incoming.valid()) {
            throw std::logic_error("invalid incoming waveform");
        }

        impl->context.incoming_waveform = &incoming;
        impl->context.accumulate_transmitted_energy = final_slice;

        model.initialize_cell_state(population);
        check_cvodes(CVodeReInit(impl->cvodes, static_cast<sunrealtype>(t_start), impl->population),
            "CVodeReInit failed");

        N_VConst(0.0, impl->transmitted_energy);
        check_cvodes(
            CVodeQuadReInit(impl->cvodes, impl->transmitted_energy), "CVodeQuadReInit failed");

        check_cvodes(CVodeSetQuadErrCon(impl->cvodes, final_slice ? SUNTRUE : SUNFALSE),
            "CVodeSetQuadErrCon failed");

        check_cvodes(CVodeSetStopTime(impl->cvodes, static_cast<sunrealtype>(t_end)),
            "CVodeSetStopTime failed");

        PulseWaveform outgoing;
        outgoing.times = incoming.times;
        outgoing.log_intensity.resize(incoming.times.size());
        outgoing.log_sensitivities.assign(incoming.times.size() * fit_parameter_count, 0.0);

        for (std::size_t knot = 0; knot < incoming.times.size(); ++knot) {
            if (knot != 0) {
                sunrealtype t_reached { };
                check_cvodes_evolve(
                    CVode(impl->cvodes, static_cast<sunrealtype>(incoming.times[knot]),
                        impl->population, &t_reached, CV_NORMAL),
                    "CVode failed");
            }

            const std::span<const double, n_levels> n { population_data, n_levels };
            const double alpha = absorption_coefficient(n, model.model_params);
            const double log_intensity = incoming.log_intensity[knot] - alpha * model.dz_cm;
            if (!std::isfinite(log_intensity)) {
                throw std::runtime_error("transmitted log intensity is not finite");
            }

            outgoing.log_intensity[knot] = log_intensity;
        }

        if (!outgoing.valid()) {
            throw std::logic_error("invalid transmitted waveform");
        }

        incoming = std::move(outgoing);
    }

    sunrealtype t_reached { };
    check_cvodes(
        CVodeGetQuad(impl->cvodes, &t_reached, impl->transmitted_energy), "CVodeGetQuad failed");

    energy_data = N_VGetArrayPointer(impl->transmitted_energy);
    if (energy_data == nullptr || !(energy_data[0] >= 0.0) || !std::isfinite(energy_data[0])) {
        throw std::runtime_error("transmitted energy is invalid");
    }

    return energy_data[0] / incident_fluence;
}

} // namespace zfit
