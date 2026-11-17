#include "openmc/tallies/filter_complex_part.h"

#include "openmc/hdf5_interface.h"
#include "openmc/vector.h"

namespace openmc {

void ComplexPartFilter::from_xml(pugi::xml_node node)
{
  n_bins_ = 2;
}

void ComplexPartFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  match.bins_.push_back(0);
  match.weights_.push_back(1.0);
  match.bins_.push_back(1);
  match.weights_.push_back(1.0);
}

void ComplexPartFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  vector<std::string> bins {"real", "imag"};
  write_dataset(filter_group, "bins", bins);
}

std::string ComplexPartFilter::text_label(int bin) const
{
  return bin == 0 ? "Complex part: real" : "Complex part: imag";
}

} // namespace openmc
