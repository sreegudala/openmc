#ifndef OPENMC_TALLIES_FILTER_NOISE_PHASE_H
#define OPENMC_TALLIES_FILTER_NOISE_PHASE_H

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Separates the static and noise phases of a noise run with a sampled noise
//! source (frequency-domain noise mode).
//!
//! Bin 0 ("static") scores the static sampling generation, which gives the
//! static flux phi_0 in the same normalisation as the noise result; bin 1
//! ("noise") scores the noise transport. Only tallies with this filter score
//! during the static sampling generation; no tally scores during the plain
//! static generations.
//==============================================================================

class NoisePhaseFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  NoisePhaseFilter() { n_bins_ = 2; }
  ~NoisePhaseFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "noisephase"; }
  FilterType type() const override { return FilterType::NOISE_PHASE; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_NOISE_PHASE_H
