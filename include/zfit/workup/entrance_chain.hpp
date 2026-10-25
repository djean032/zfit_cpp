#ifndef ZFIT_WORKUP_ENTRANCE_CHAIN_HPP
#define ZFIT_WORKUP_ENTRANCE_CHAIN_HPP

#include <string>

namespace zfit::workup {

// Computed entrance optical chain: the fraction of the energymeter-measured
// energy that reaches the molecular sample plane. The meter samples the beam
// in free air directly in front of the cuvette, so the chain is exactly the
// single entrance window at normal incidence: T(air -> window) times
// T(window -> solvent) from Fresnel's law (polarization only enters at
// second order in the incidence angle, negligible for a focusing cone).
// Solvent linear absorption is settled negligible and is not included;
// refraction through the window is settled exactly absorbed into the
// effective (waist, M2) calibration.
//
// Users never enter kappa: the headless layer takes raw refractive indices
// and their catalogue tolerances (material presets with n and tolerance are
// a future GUI concern). The uncertainty comes only from those tolerances
// propagated through the Fresnel derivatives; a characterization floor for
// placement/drift effects is explicitly not modeled (no floor term).
struct EntranceChainInputs {
    // Refractive indices at the laser wavelength, raw values.
    double window_refractive_index { };
    double solvent_refractive_index { };
    // Catalogue tolerances of those indices, as standard deviations
    // (absolute, not relative). Both required and positive: claiming an
    // exact material index is refused, never silently exact.
    double window_index_std { };
    double solvent_index_std { };
    // Source of the index values and tolerances (catalogue or datasheet
    // reference). Required.
    std::string index_source_id;
};

struct EntranceChainValue {
    double transmission { }; // kappa, in (0, 1]
    double standard_error { }; // standard deviation of kappa, > 0
    std::string index_source_id;
};

// Computes kappa and its uncertainty. Throws std::invalid_argument on
// non-finite/non-positive indices, non-finite/non-positive tolerances, or
// an empty index source id.
EntranceChainValue compute_entrance_chain(const EntranceChainInputs& inputs);

} // namespace zfit::workup

#endif
