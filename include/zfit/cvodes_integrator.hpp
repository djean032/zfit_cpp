#pragma once

#include <cstddef>
#include <memory>

#include "model.hpp"

namespace zfit {

class CvodesIntegrator;

struct CvodesContext {
    const MolModel* model { };
    const PulseWaveform* incoming_waveform { };
    bool accumulate_transmitted_energy { };
};

class CvodesIntegrator {
public:
    struct Impl;

    explicit CvodesIntegrator(MolModel& model);
    ~CvodesIntegrator();

    CvodesIntegrator(const CvodesIntegrator&) = delete;

    CvodesIntegrator& operator=(const CvodesIntegrator&) = delete;

    double simulate(std::size_t scan_index);

    std::unique_ptr<Impl> impl;
};

} // namespace zfit
