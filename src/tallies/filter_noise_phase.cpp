#include "openmc/tallies/filter_noise_phase.h"

#include "openmc/hdf5_interface.h"
#include "openmc/noise.h"
#include "openmc/vector.h"

namespace openmc {

void NoisePhaseFilter::from_xml(pugi::xml_node node)
{
  n_bins_ = 2;
}

void NoisePhaseFilter::get_all_bins(
  const Particle& p, TallyEstimator estimator, FilterMatch& match) const
{
  match.bins_.push_back(noise::static_phase() ? 0 : 1);
  match.weights_.push_back(1.0);
}

void NoisePhaseFilter::to_statepoint(hid_t filter_group) const
{
  Filter::to_statepoint(filter_group);
  vector<std::string> bins {"static", "noise"};
  write_dataset(filter_group, "bins", bins);
}

std::string NoisePhaseFilter::text_label(int bin) const
{
  return bin == 0 ? "Noise phase: static" : "Noise phase: noise";
}

} // namespace openmc
