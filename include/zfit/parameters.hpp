#pragma once

#include <cstddef>
namespace zfit {

struct ModelParameters {
    double sigma_g { }; // S0 -> S1
    double sigma_s { }; // S1 -> S2
    double sigma_t { }; // T1 -> T2

    double tau10 { }; // S1 -> S0
    double tau21 { }; // S2 -> S1
    double tau13 { }; // S1 -> T1
    double tau30 { }; // T1 -> S0
    double tau43 { }; // T2 -> T1
};

enum class FitParameter {
    Sigma_G,
    Sigma_S,
    Sigma_T,
    Tau10,
    Tau21,
    Tau13,
    Tau30,
    Tau43,
    Count,
};

inline constexpr std::size_t fit_parameter_count = static_cast<std::size_t>(FitParameter::Count);

} // namespace zfit
