#pragma once

#include <vector>

namespace zfit {

struct Experiment {
    std::vector<double> z;
    std::vector<double> transmission;
    // std::vector<double> sigma;
};
} // namespace zfit
