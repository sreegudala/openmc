//! \file noise.h
//! Frequency-domain neutron noise transport with complex particle weights
//! (multigroup, fixed source). Collision-based scheme of Rouchon et al.:
//! flights use Sigma_t + eta * omega / v and a "copy" reaction multiplies the
//! weight by (1 - i / eta).

#ifndef OPENMC_NOISE_H
#define OPENMC_NOISE_H

#include "hdf5.h"

namespace openmc {

class Particle;

namespace noise {

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
void roulette(Particle& p);

//! Bank noise fission sites for the next fission generation
//!
//! n = floor(nuSigma_f / (Sigma_t~ k) + xi) sites, each inheriting the
//! parent's complex weight; delayed sites are multiplied by
//! lambda / (lambda + i omega). No emission time is sampled.
void create_fission_sites(Particle& p);

//! Transport one batch as a loop over fission generations
//!
//! Generation 0 is the external source. Fission sites of generation g form
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
