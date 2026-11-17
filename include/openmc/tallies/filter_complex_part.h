#ifndef OPENMC_TALLIES_FILTER_COMPLEX_PART_H
#define OPENMC_TALLIES_FILTER_COMPLEX_PART_H

#include "openmc/tallies/filter.h"

namespace openmc {

//==============================================================================
//! Splits scores into the real and imaginary parts of a complex particle
//! weight (frequency-domain noise mode).
//!
//! The filter always has two bins, 0 = real and 1 = imaginary. It matches
//! both bins with unit filter weight; the tally drivers select which weight
//! component a given bin combination scores (see Particle::wgt_tally()).
//==============================================================================

class ComplexPartFilter : public Filter {
public:
  //----------------------------------------------------------------------------
  // Constructors, destructors

  ComplexPartFilter() { n_bins_ = 2; }
  ~ComplexPartFilter() = default;

  //----------------------------------------------------------------------------
  // Methods

  std::string type_str() const override { return "complexpart"; }
  FilterType type() const override { return FilterType::COMPLEX_PART; }

  void from_xml(pugi::xml_node node) override;

  void get_all_bins(const Particle& p, TallyEstimator estimator,
    FilterMatch& match) const override;

  void to_statepoint(hid_t filter_group) const override;

  std::string text_label(int bin) const override;
};

} // namespace openmc
#endif // OPENMC_TALLIES_FILTER_COMPLEX_PART_H
