#include "openmc/physics_common.h"

#include "openmc/random_lcg.h"
#include "openmc/settings.h"

namespace openmc {

//==============================================================================
// RUSSIAN_ROULETTE
//==============================================================================

void russian_roulette(Particle& p, double weight_survive)
{
  // Positive real weights use the original arithmetic. Complex weights, and
  // negative real weights that arise in noise mode (e.g. (1 - i)^4 = -4),
  // use the norm-based rule, which preserves the phase/sign.
  if (p.wgt_im() == 0.0 && p.wgt() > 0.0) {
    if (weight_survive * prn(p.current_seed()) < p.wgt()) {
      p.wgt() = weight_survive;
    } else {
      p.kill();
    }
  } else {
    // Complex weight (noise mode): norm-based roulette that preserves the
    // phase of the weight (Fauvel et al. 2025). The particle survives with
    // probability |w| / weight_survive and takes weight w * weight_survive/|w|.
    double w_abs = p.wgt_abs();
    if (weight_survive * prn(p.current_seed()) < w_abs) {
      p.scale_wgt(weight_survive / w_abs);
    } else {
      p.kill();
    }
  }
}

void apply_russian_roulette(Particle& p)
{
  // Exit if survival biasing is turned off
  if (!settings::survival_biasing)
    return;

  // if survival normalization is on, use normalized weight cutoff and
  // normalized weight survive
  // For complex weights (noise mode) roulette acts on the modulus
  double w = (p.wgt_im() == 0.0 && p.wgt() > 0.0) ? p.wgt() : p.wgt_abs();
  if (settings::survival_normalization) {
    if (w < settings::weight_cutoff * p.wgt_born()) {
      russian_roulette(p, settings::weight_survive * p.wgt_born());
    }
  } else if (w < settings::weight_cutoff) {
    russian_roulette(p, settings::weight_survive);
  }
}
} // namespace openmc
