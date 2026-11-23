//! \file noise.h
//! Frequency-domain neutron noise transport with complex particle weights
//! (multigroup, fixed source). Collision-based scheme of Rouchon et al.:
//! flights use Sigma_t + eta * omega / v and a "copy" reaction multiplies the
//! weight by (1 - i / eta).

#ifndef OPENMC_NOISE_H
#define OPENMC_NOISE_H

#include <complex>

#include "hdf5.h"
#include "pugixml.hpp"

#include "openmc/vector.h"

namespace openmc {

class Particle;

namespace noise {

//==============================================================================
// Sampled noise source (design doc Stage 6)
//==============================================================================

//! Planar interface vibrating as x0 + eps sin(omega_0 t) along one axis
//! (design doc Stage 8, amendment 2026-11-23). Harmonic n has
//! dSigma_x = (Sigma_x^L - Sigma_x^R) c_n(u), u = (x - x0) / eps, with L the
//! material on the negative side and c_n = (2/n) sin(n arccos u) e^{-i n pi/2}.
struct Vibration {
  int axis {-1};            //!< 0, 1, 2 for x, y, z; -1: not a vibration
  double position {0.0};    //!< interface position x0 [cm]
  double amplitude {0.0};   //!< eps [cm]
  int harmonic {0};         //!< n >= 1
  int32_t negative_id {-1}; //!< material ID on the negative side (L)
  int32_t positive_id {-1}; //!< material ID on the positive side (R)
  int negative {-1};        //!< material index of L
  int positive {-1};        //!< material index of R
};

//! Cross-section perturbation of a set of cells, sampled as a noise source at
//! static collisions (Rouchon et al.; table of design doc Section 2.4).
//! Relative complex amplitudes per group; an empty vector means zero. A
//! vibration instead derives its amplitudes from two materials' data.
struct Perturbation {
  vector<int32_t> cell_ids;             //!< cell IDs from the input
  vector<std::complex<double>> total;   //!< dSigma_t / Sigma_t per group
  vector<std::complex<double>> scatter; //!< dSigma_s / Sigma_s per group
  vector<std::complex<double>> fission; //!< dSigma_f / Sigma_f per group
  Vibration vibration;                  //!< used when vibration.axis >= 0

  bool is_vibration() const { return vibration.axis >= 0; }
};

//! Harmonic coefficient c_n(u) of a vibrating interface, zero for |u| >= 1
std::complex<double> vibration_coefficient(int n, double u);

extern vector<Perturbation> perturbations; //!< perturbations sampled as source
extern int nskip;           //!< static generations per noise batch (>= 1)
extern int static_inactive; //!< static generations before the first batch

//! Read a <perturbation> element of the <noise> settings
void read_perturbation(pugi::xml_node node);

//! Reset the sampled-source settings to their defaults
void reset_settings();

//! True while the static generations of a sampled-source run are transported
bool static_phase();

//! Validate the model for noise mode and precompute frequency-dependent data.
//! Called from openmc_simulation_init when settings::noise_mode is on.
void initialize();

//! Noise-mode replacement for sample_reaction in multigroup transport
//!
//! Sequence at a collision (following the published collision-based scheme):
//! fission branching, copy reaction with probability Sigma_w / Sigma_t~,
//! implicit capture w <- w Sigma_s / Sigma_t~, norm-based roulette, scatter.
void sample_reaction(Particle& p);

//! Norm-based Russian roulette on |w| using the survival-biasing cutoffs.
//! With weight windows on, it is applied only outside the window mesh
//! (from collision_mg); inside, the windows take its place.
void roulette(Particle& p);

//! Weight scale of the current noise source: 1 for an external source,
//! the batch's mean source modulus for a sampled source. Roulette
//! thresholds and weight windows are in units of it.
double source_weight_scale();

//! Bank noise fission sites for the next fission generation
//!
//! n = floor(nuSigma_f / (Sigma_t~ k) + xi) sites, each inheriting the
//! parent's complex weight; delayed sites are multiplied by
//! lambda / (lambda + i omega). No emission time is sampled.
void create_fission_sites(Particle& p);

//! Transport one batch as a loop over fission generations
//!
//! Generation 0 is the external source or, when perturbations are given, the
//! noise source sampled at the collisions of one static generation (after
//! nskip - 1 plain static generations, and static_inactive more before the
//! first batch). Fission sites of generation g form
//! generation g + 1; copies are processed within the generation. The loop
//! ends when the fission bank is empty. It aborts when the number of
//! generations exceeds settings::noise_max_generations or when the sum of
//! |w| over the bank exceeds settings::noise_growth_limit times its value in
//! the first fission generation.
void transport();

//! Write per-generation statistics to a statepoint group "noise"
void write_statepoint(hid_t file_id);

//! Reset accumulated generation statistics
void reset();

} // namespace noise
} // namespace openmc

#endif // OPENMC_NOISE_H
