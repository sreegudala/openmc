#include "openmc/noise.h"

#include <algorithm>
#include <cmath>
#include <complex>

#include <fmt/core.h>

#include "openmc/bank.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/math_functions.h"
#include "openmc/mesh.h"
#include "openmc/message_passing.h"
#include "openmc/mgxs_interface.h"
#include "openmc/openmp_interface.h"
#include "openmc/particle.h"
#include "openmc/physics_common.h"
#include "openmc/physics_mg.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/source.h"
#include "openmc/tallies/filter_noise_phase.h"
#include "openmc/tallies/tally.h"
#include "openmc/vector.h"
#include "openmc/xml_interface.h"

namespace openmc {
namespace noise {

//==============================================================================
// Global variables
//==============================================================================

vector<Perturbation> perturbations;
int nskip {20};
int static_inactive {50};

//==============================================================================
// Module state
//==============================================================================

namespace {

// Per-thread banks of fission sites for the next fission generation
vector<vector<SourceSite>> thread_fission_banks;

// Offset added to particle ids of fission-generation tracks so that their
// random number seeds never coincide with those of source particles or of
// tracks in other generations or batches
int64_t track_offset {0};

// Per-generation statistics accumulated over batches: sum of |w| of the bank
// entering each generation, number of sites, and number of batches reaching
// that generation; generations per batch
vector<double> gen_sum_abs_wgt;
vector<double> gen_sum_abs_wgt_sq;
vector<double> gen_sum_wgt_re;
vector<double> gen_sum_wgt_im;
vector<double> gen_n_sites;
vector<int64_t> gen_n_batches;
vector<int> batch_generations;

// Cancellation diagnostics accumulated over all generations and batches:
// rendezvous count, sites entering cancellation, occupied bins, sum of |w|
// before and after averaging, and sites removed by the bank roulette
double cancel_n_rendezvous {0.0};
double cancel_n_sites {0.0};
double cancel_n_occupied {0.0};
double cancel_abs_before {0.0};
double cancel_abs_after {0.0};
double cancel_n_removed {0.0};

// Sampled noise source (Stage 6). The static phase runs the noise collision
// physics at omega = 0 (real weights), which is a static fixed-point game with
// fission divided by noise_keff; the static bank is renormalised to
// n_particles sites every generation.
vector<int> cell_perturbation;                  // cell index -> perturbation
vector<vector<SourceSite>> thread_source_banks; // sampled noise source sites
vector<SourceSite> static_bank;                 // static source sites
bool in_static {false};                         // static generation running
bool sampling {false};                          // static sampling generation
double omega_noise {0.0};  // noise frequency while the static phase runs
double weight_scale {1.0}; // roulette thresholds are scaled by this factor
vector<double> static_k_gen;      // static k estimate of every generation
vector<double> batch_n_source;    // sampled source sites per batch
vector<double> batch_scale;       // mean source modulus per batch
vector<double> batch_source_sum;  // [batch][perturbation][re, im]
vector<double> source_sum_thread; // [thread][perturbation][re, im]

} // namespace

//==============================================================================
// Settings for the sampled noise source
//==============================================================================

namespace {

vector<std::complex<double>> read_complex_array(
  pugi::xml_node node, const char* name)
{
  vector<std::complex<double>> out;
  if (!check_for_node(node, name))
    return out;
  auto v = get_node_array<double>(node, name);
  if (v.size() % 2 != 0) {
    fatal_error(fmt::format("Noise perturbation '{}' must list (real, "
                            "imaginary) pairs.",
      name));
  }
  for (size_t i = 0; i < v.size(); i += 2)
    out.emplace_back(v[i], v[i + 1]);
  return out;
}

} // namespace

void read_perturbation(pugi::xml_node node)
{
  Perturbation pt;
  if (!check_for_node(node, "cells"))
    fatal_error("A noise perturbation must list its cells.");
  pt.cell_ids = get_node_array<int32_t>(node, "cells");
  pt.total = read_complex_array(node, "total");
  pt.scatter = read_complex_array(node, "scatter");
  pt.fission = read_complex_array(node, "fission");
  perturbations.push_back(pt);
}

void reset_settings()
{
  perturbations.clear();
  nskip = 20;
  static_inactive = 50;
}

bool static_phase()
{
  return in_static;
}

void reset()
{
  track_offset = static_cast<int64_t>(settings::n_batches) *
                 settings::gen_per_batch * settings::n_particles;
  gen_sum_abs_wgt.clear();
  gen_sum_abs_wgt_sq.clear();
  gen_sum_wgt_re.clear();
  gen_sum_wgt_im.clear();
  gen_n_sites.clear();
  gen_n_batches.clear();
  batch_generations.clear();
  cancel_n_rendezvous = 0.0;
  cancel_n_sites = 0.0;
  cancel_n_occupied = 0.0;
  cancel_abs_before = 0.0;
  cancel_abs_after = 0.0;
  cancel_n_removed = 0.0;
  static_bank.clear();
  static_k_gen.clear();
  batch_n_source.clear();
  batch_scale.clear();
  batch_source_sum.clear();
  in_static = false;
  sampling = false;
  weight_scale = 1.0;
}

void initialize()
{
  reset();

  if (settings::event_based) {
    warning("Event-based transport is not available in noise mode; "
            "history-based transport will be used.");
    settings::event_based = false;
  }

  // Weight windows compare and split the real part of the weight only: the
  // cutoff kills every particle with Re w <= cutoff and split copies lose
  // the imaginary part (noise_dev/stage7/stage7_prep_log.md, N3c).
  if (settings::weight_windows_on) {
    fatal_error("Weight windows are not yet supported in noise mode (they "
                "act on the real part of complex weights only).");
  }

  // Inverse velocity must be supplied: OpenMC otherwise silently substitutes
  // an estimate from the group boundaries.
  if (settings::noise_omega > 0.0) {
    for (int i = 0; i < data::mg.macro_xs_.size(); ++i) {
      if (!data::mg.macro_xs_[i].has_inverse_velocity()) {
        fatal_error(fmt::format("Noise mode at a non-zero frequency requires "
                                "inverse-velocity data, which are missing for "
                                "material '{}'.",
          data::mg.macro_xs_[i].name));
      }
    }
  }

  // Void regions carry no group speed, so the i omega / v term cannot be
  // applied there.
  for (const auto& c : model::cells) {
    if (c->type_ == Fill::MATERIAL) {
      for (auto m : c->material_) {
        if (m == MATERIAL_VOID && settings::noise_omega > 0.0) {
          warning(fmt::format("Cell {} is void; the i omega / v term of the "
                              "noise operator is not applied in void.",
            c->id_));
        }
      }
    }
  }

  // Sampled noise source: map cells to perturbations and broadcast scalar
  // amplitudes to every group
  cell_perturbation.assign(model::cells.size(), -1);
  int G = data::mg.num_energy_groups_;
  for (int ip = 0; ip < perturbations.size(); ++ip) {
    auto& pt = perturbations[ip];
    for (auto id : pt.cell_ids) {
      auto it = model::cell_map.find(id);
      if (it == model::cell_map.end()) {
        fatal_error(fmt::format(
          "Cell {} of a noise perturbation does not exist.", id));
      }
      if (cell_perturbation[it->second] >= 0) {
        fatal_error(fmt::format(
          "Cell {} belongs to more than one noise perturbation.", id));
      }
      cell_perturbation[it->second] = ip;
    }
    for (auto* v : {&pt.total, &pt.scatter, &pt.fission}) {
      if (v->empty()) {
        v->assign(G, 0.0);
      } else if (v->size() == 1) {
        v->assign(G, (*v)[0]);
      } else if (v->size() != G) {
        fatal_error(fmt::format("Noise perturbation amplitudes must have 1 "
                                "or {} groups.",
          G));
      }
    }
  }
  if (!perturbations.empty() && mpi::n_procs > 1) {
    fatal_error("A sampled noise source (perturbations) is not yet "
                "available with more than one MPI process.");
  }
}

namespace {

// Roulette thresholds. A sampled noise source has weights of order
// |dSigma / Sigma|; the thresholds are scaled by the batch's mean source
// modulus, which is equivalent to normalising the source weights by it and
// multiplying the tallies back (roulette is homogeneous in the weight scale).
double cutoff()
{
  return settings::weight_cutoff * weight_scale;
}
double survive()
{
  return settings::weight_survive * weight_scale;
}

} // namespace

void roulette(Particle& p)
{
  if (p.wgt_abs() < cutoff()) {
    russian_roulette(p, survive());
  }
}

void create_fission_sites(Particle& p)
{
  const auto& xs = p.macro_xs();
  double nu_t = xs.nu_fission / (xs.total * settings::noise_keff);
  int nu = static_cast<int>(nu_t);
  if (prn(p.current_seed()) <= (nu_t - nu))
    nu++;
  if (nu == 0)
    return;

  auto& macro = data::mg.macro_xs_[p.material()];
  int t = p.mg_xs_cache().t;
  int a = p.mg_xs_cache().a;
  double omega = settings::noise_omega;
  auto& bank = thread_fission_banks[thread_num()];
  int64_t parent = simulation::work_index[mpi::rank] + p.current_work();

  for (int i = 0; i < nu; ++i) {
    SourceSite site;
    site.r = p.r();
    site.particle = ParticleType::neutron();
    site.time = p.time();
    site.wgt = p.wgt();
    site.wgt_im = p.wgt_im();

    // Isotropic emission
    double mu = 2. * prn(p.current_seed()) - 1.;
    double phi = 2. * PI * prn(p.current_seed());
    site.u.x = mu;
    site.u.y = std::sqrt(1. - mu * mu) * std::cos(phi);
    site.u.z = std::sqrt(1. - mu * mu) * std::sin(phi);

    // Prompt or delayed family and outgoing group (real probabilities)
    int dg;
    int gout;
    macro.sample_fission_energy(p.g(), dg, gout, p.current_seed(), t, a);
    site.E = gout;
    site.delayed_group = dg + 1;

    // Delayed neutrons carry lambda / (lambda + i omega)
    if (dg != -1 && omega != 0.0) {
      double lam =
        macro.get_xs(MgxsType::DECAY_RATE, p.g(), nullptr, nullptr, &dg, t, a);
      double d = lam * lam + omega * omega;
      double f_re = lam * lam / d;
      double f_im = -lam * omega / d;
      double wr = site.wgt;
      double wi = site.wgt_im;
      site.wgt = wr * f_re - wi * f_im;
      site.wgt_im = wr * f_im + wi * f_re;
    }

    // Norm-based roulette of the site. Sites inherit the parent's weight, so
    // without this a lineage with |nu Sigma_f / (Sigma_t~ k)| ~ 1 never dies
    // and its weight only shrinks with each delayed emission.
    double w_abs = std::hypot(site.wgt, site.wgt_im);
    if (w_abs < cutoff()) {
      if (survive() * prn(p.current_seed()) < w_abs) {
        double f = survive() / w_abs;
        site.wgt *= f;
        site.wgt_im *= f;
      } else {
        continue;
      }
    }

    site.parent_id = parent;
    site.progeny_id = p.n_progeny()++;
    bank.push_back(site);
  }
}

namespace {

//! Bank noise source sites at a static collision of weight w inside a
//! perturbed cell (design doc Section 2.4): a copy with -dSigma_t/Sigma_t w,
//! fission sites with w dSigma_f/Sigma_f (delayed ones times
//! lambda/(lambda + i omega)), and a sampled scatter with
//! w yield (Sigma_s/Sigma_t) dSigma_s/Sigma_s.
void sample_source_sites(Particle& p, int ip)
{
  using cplx = std::complex<double>;
  const auto& pt = perturbations[ip];
  int g = p.g();
  double w = p.wgt();
  const auto& xs = p.macro_xs();
  auto& macro = data::mg.macro_xs_[p.material()];
  int t = p.mg_xs_cache().t;
  int a = p.mg_xs_cache().a;
  auto& bank = thread_source_banks[thread_num()];
  int64_t parent = simulation::work_index[mpi::rank] + p.current_work();
  int64_t i_sum = thread_num() * perturbations.size() + ip;
  double* sum = &source_sum_thread[2 * i_sum];

  auto push = [&](Position r, Direction u, int gout, int dg, cplx wc) {
    SourceSite site;
    site.r = r;
    site.u = u;
    site.E = gout;
    site.particle = ParticleType::neutron();
    site.time = p.time();
    site.wgt = wc.real();
    site.wgt_im = wc.imag();
    site.delayed_group = dg + 1;
    site.parent_id = parent;
    site.progeny_id = p.n_progeny()++;
    bank.push_back(site);
    sum[0] += wc.real();
    sum[1] += wc.imag();
  };

  // Copy: -dSigma_t psi_0
  if (pt.total[g] != 0.0) {
    push(p.r(), p.u(), g, -1, -pt.total[g] * w);
  }

  // Fission: + chi dnuSigma_f phi_0 / k
  if (pt.fission[g] != 0.0 && macro.fissionable) {
    double nu_t = xs.nu_fission / (xs.total * settings::noise_keff);
    int nu = static_cast<int>(nu_t);
    if (prn(p.current_seed()) <= (nu_t - nu))
      nu++;
    for (int i = 0; i < nu; ++i) {
      int dg;
      int gout;
      macro.sample_fission_energy(g, dg, gout, p.current_seed(), t, a);
      cplx wc = pt.fission[g] * w;
      if (dg != -1 && omega_noise != 0.0) {
        double lam =
          macro.get_xs(MgxsType::DECAY_RATE, g, nullptr, nullptr, &dg, t, a);
        wc *= lam / cplx(lam, omega_noise);
      }
      double mu = 2. * prn(p.current_seed()) - 1.;
      double phi = 2. * PI * prn(p.current_seed());
      Direction u {mu, std::sqrt(1. - mu * mu) * std::cos(phi),
        std::sqrt(1. - mu * mu) * std::sin(phi)};
      push(p.r(), u, gout, dg, wc);
    }
  }

  // Scatter: + dSigma_s (kernel) psi_0
  double sigma_s = xs.total - xs.absorption;
  if (pt.scatter[g] != 0.0 && sigma_s > 0.0) {
    int gout;
    double mu;
    double mult = 1.0;
    macro.sample_scatter(g, gout, mu, mult, p.current_seed(), t, a);
    Direction u = rotate_angle(p.u(), mu, nullptr, p.current_seed());
    push(p.r(), u, gout, -1, pt.scatter[g] * (w * mult * sigma_s / xs.total));
  }
}

} // namespace

void sample_reaction(Particle& p)
{
  // Sampled noise source at static collisions in perturbed cells
  if (sampling) {
    int ip = cell_perturbation[p.lowest_coord().cell()];
    if (ip >= 0)
      sample_source_sites(p, ip);
  }

  // Fission branching: sites go to the next fission generation
  if (data::mg.macro_xs_[p.material()].fissionable &&
      settings::create_fission_neutrons) {
    noise::create_fission_sites(p);
  }

  const auto& xs = p.macro_xs();
  double sigma_t = xs.total; // includes sigma_w
  double sigma_w = xs.noise;

  // Copy reaction: same position, direction and group; weight w (1 - i/eta)
  if (sigma_w > 0.0 && prn(p.current_seed()) * sigma_t < sigma_w) {
    // Guard against a divergent within-generation game (extended collision
    // radius >= 1, i.e. omega above the collision-based divergence cut-off)
    if (p.local_secondary_bank().size() >= settings::max_secondaries) {
      fatal_error(fmt::format(
        "Noise mode: the secondary bank of particle {} exceeded "
        "max_secondaries ({}). The within-generation collision game is "
        "likely divergent at this frequency; increase eta.",
        p.id(), settings::max_secondaries));
    }
    double eta = settings::noise_eta;
    double wr = p.wgt();
    double wi = p.wgt_im();
    p.create_secondary(
      wr + wi / eta, p.u(), p.E(), p.type(), wi - wr / eta);
  }

  // Implicit capture: the continuing particle carries the scattered weight
  double sigma_s = sigma_t - sigma_w - xs.absorption;
  p.scale_wgt(sigma_s / sigma_t);

  // Norm-based roulette preserving the phase
  roulette(p);
  if (!p.alive())
    return;

  scatter(p);
}

//==============================================================================
// Fission-generation loop
//==============================================================================

namespace {

void transport_tracks(int64_t n, bool from_bank)
{
#pragma omp parallel
  {
    Particle p;
#pragma omp for schedule(runtime)
    for (int64_t i = 1; i <= n; ++i) {
      initialize_particle_track(p, i, from_bank);
      if (from_bank) {
        // Unique id and seed for this fission-generation track
        p.id() = track_offset + simulation::work_index[mpi::rank] + i;
        init_particle_seeds(p.id(), p.seeds());
        p.event_revive_from_secondary(
          simulation::shared_secondary_bank_read[i - 1]);
        p.wgt_born() = p.wgt();
      }
      transport_history_based_single_particle(p);
    }
  }
}

//! Gather per-thread banks into one vector in canonical order (parent global
//! index, then progeny index), so results do not depend on the thread count
vector<SourceSite> gather_canonical(vector<vector<SourceSite>>& banks)
{
  int64_t n = 0;
  for (const auto& b : banks)
    n += b.size();
  vector<SourceSite> all;
  all.reserve(n);
  for (auto& b : banks) {
    all.insert(all.end(), b.begin(), b.end());
    b.clear();
  }
  std::sort(all.begin(), all.end(), [](const SourceSite& a, const SourceSite& b) {
    return a.parent_id < b.parent_id ||
           (a.parent_id == b.parent_id && a.progeny_id < b.progeny_id);
  });
  return all;
}

void set_bank(SharedArray<SourceSite>& bank, const vector<SourceSite>& sites)
{
  bank.resize(0);
  bank.extend_uninitialized(sites.size());
  for (int64_t i = 0; i < sites.size(); ++i)
    bank[i] = sites[i];
}

//! Gather thread fission banks into shared_secondary_bank_write
void collect_fission_banks()
{
  set_bank(simulation::shared_secondary_bank_write,
    gather_canonical(thread_fission_banks));
}

//==============================================================================
// Static phase of a sampled-source run
//==============================================================================

//! Active tally lists, saved while the static phase restricts them
struct TallyLists {
  vector<int> analog, tracklength, timed, collision, meshsurf, surface, pulse;
};

TallyLists save_tally_lists()
{
  return {model::active_analog_tallies, model::active_tracklength_tallies,
    model::active_timed_tracklength_tallies, model::active_collision_tallies,
    model::active_meshsurf_tallies, model::active_surface_tallies,
    model::active_pulse_height_tallies};
}

void restore_tally_lists(const TallyLists& l)
{
  model::active_analog_tallies = l.analog;
  model::active_tracklength_tallies = l.tracklength;
  model::active_timed_tracklength_tallies = l.timed;
  model::active_collision_tallies = l.collision;
  model::active_meshsurf_tallies = l.meshsurf;
  model::active_surface_tallies = l.surface;
  model::active_pulse_height_tallies = l.pulse;
}

//! Keep only tallies with a NoisePhaseFilter (all = false: keep none)
void restrict_tally_lists(const TallyLists& l, bool phase_tallies)
{
  auto keep = [&](const vector<int>& in) {
    vector<int> out;
    if (phase_tallies) {
      for (int i : in) {
        if (model::tallies[i]->get_filter<NoisePhaseFilter>())
          out.push_back(i);
      }
    }
    return out;
  };
  model::active_analog_tallies = keep(l.analog);
  model::active_tracklength_tallies = keep(l.tracklength);
  model::active_timed_tracklength_tallies = keep(l.timed);
  model::active_collision_tallies = keep(l.collision);
  model::active_meshsurf_tallies = keep(l.meshsurf);
  model::active_surface_tallies = keep(l.surface);
  model::active_pulse_height_tallies = keep(l.pulse);
}

//! Transport one static generation from static_bank and renormalise its
//! fission sites to n_particles sites of unit weight (systematic resampling
//! by weight, canonical order).
void static_generation()
{
  int64_t n = static_bank.size();
  set_bank(simulation::shared_secondary_bank_read, static_bank);
  calculate_work(n);
  transport_tracks(n, true);
  track_offset += n;
  simulation::shared_secondary_bank_read.clear();

  auto sites = gather_canonical(thread_fission_banks);
  double W = 0.0;
  for (const auto& s : sites)
    W += s.wgt;
  static_k_gen.push_back(W / n * settings::noise_keff);
  if (sites.empty() || W <= 0.0) {
    fatal_error("Noise mode: the static fission source died out.");
  }

  int64_t n_new = settings::n_particles;
  uint64_t seed = init_seed(track_offset, STREAM_SOURCE);
  double step = W / n_new;
  double next = prn(&seed) * step;
  double cum = 0.0;
  static_bank.clear();
  static_bank.reserve(n_new);
  int64_t j = 0;
  for (int64_t i = 0; i < n_new; ++i) {
    while (cum + sites[j].wgt <= next && j + 1 < sites.size()) {
      cum += sites[j].wgt;
      ++j;
    }
    SourceSite s = sites[j];
    s.wgt = 1.0;
    s.wgt_im = 0.0;
    static_bank.push_back(s);
    next += step;
  }
}

//! Run the static phase for one batch and leave the sampled noise source in
//! shared_secondary_bank_read. Returns the number of source sites.
int64_t sample_noise_source()
{
  omega_noise = settings::noise_omega;
  settings::noise_omega = 0.0;
  weight_scale = 1.0;
  in_static = true;
  auto lists = save_tally_lists();
  restrict_tally_lists(lists, false);

  int n_plain = nskip - 1;
  if (static_bank.empty()) {
    // Initial static source from the external source distribution
    int64_t n = settings::n_particles;
    static_bank.resize(n);
    for (int64_t i = 0; i < n; ++i) {
      uint64_t seed = init_seed(track_offset + i + 1, STREAM_SOURCE);
      static_bank[i] = sample_external_source(&seed);
      static_bank[i].wgt = 1.0;
      static_bank[i].wgt_im = 0.0;
    }
    track_offset += n;
    n_plain += static_inactive;
  }
  for (int i = 0; i < n_plain; ++i)
    static_generation();

  // Sampling generation: static tallies score, source sites are banked
  int n_pert = perturbations.size();
  thread_source_banks.assign(num_threads(), {});
  source_sum_thread.assign(2 * num_threads() * n_pert, 0.0);
  restrict_tally_lists(lists, true);
  sampling = true;
  static_generation();
  sampling = false;

  restore_tally_lists(lists);
  in_static = false;
  settings::noise_omega = omega_noise;

  // Source bank and its statistics
  auto src = gather_canonical(thread_source_banks);
  for (int ip = 0; ip < n_pert; ++ip) {
    double re = 0.0;
    double im = 0.0;
    for (int t = 0; t < num_threads(); ++t) {
      re += source_sum_thread[2 * (t * n_pert + ip)];
      im += source_sum_thread[2 * (t * n_pert + ip) + 1];
    }
    batch_source_sum.push_back(re);
    batch_source_sum.push_back(im);
  }
  double sum_abs = 0.0;
  for (const auto& s : src)
    sum_abs += std::hypot(s.wgt, s.wgt_im);
  weight_scale = src.empty() ? 1.0 : sum_abs / src.size();
  batch_n_source.push_back(src.size());
  batch_scale.push_back(weight_scale);
  set_bank(simulation::shared_secondary_bank_read, src);

  // Tallies are normalised per static particle of the sampling generation
  simulation::total_weight = settings::n_particles;
  return src.size();
}

//! Mesh-averaged weight cancellation of the local fission bank followed by a
//! norm-based roulette of the bank (design doc Stage 5).
//!
//! Sites are binned by mesh cell (and outgoing group when
//! settings::noise_cancellation_group, needed because sites are banked after
//! the outgoing group is sampled). The per-bin complex sum and count are
//! reduced over MPI ranks and every site in a bin takes the bin's mean weight;
//! the total complex weight per bin is preserved. Sites outside the mesh are
//! left unchanged. Returns the number of sites removed by roulette.
int64_t cancel_and_roulette(SharedArray<SourceSite>& bank, int64_t id_base)
{
  const auto& mesh = *model::meshes[settings::noise_cancellation_mesh];
  int n_groups = settings::noise_cancellation_group
                   ? data::mg.num_energy_groups_
                   : 1;
  int64_t n_bins = static_cast<int64_t>(mesh.n_bins()) * n_groups;

  // Per-bin sums: re, im, count
  vector<double> sums(3 * n_bins, 0.0);
  vector<int64_t> bin(bank.size());
  double abs_before = 0.0;
  double abs_outside = 0.0;
  for (int64_t i = 0; i < bank.size(); ++i) {
    double a = std::hypot(bank[i].wgt, bank[i].wgt_im);
    abs_before += a;
    int b = mesh.get_bin(bank[i].r);
    if (b < 0) {
      abs_outside += a;
      bin[i] = -1;
      continue;
    }
    int64_t k = static_cast<int64_t>(b) * n_groups +
                (n_groups > 1 ? static_cast<int>(bank[i].E) : 0);
    bin[i] = k;
    sums[3 * k] += bank[i].wgt;
    sums[3 * k + 1] += bank[i].wgt_im;
    sums[3 * k + 2] += 1.0;
  }
#ifdef OPENMC_MPI
  MPI_Allreduce(MPI_IN_PLACE, sums.data(), sums.size(), MPI_DOUBLE, MPI_SUM,
    mpi::intracomm);
#endif

  // Global index of the first local site (for reproducible roulette seeds)
  int64_t offset = 0;
#ifdef OPENMC_MPI
  int64_t local = bank.size();
  MPI_Exscan(&local, &offset, 1, MPI_INT64_T, MPI_SUM, mpi::intracomm);
  if (mpi::rank == 0)
    offset = 0;
#endif

  // Assign bin means and roulette
  int64_t n_keep = 0;
  for (int64_t i = 0; i < bank.size(); ++i) {
    SourceSite site = bank[i];
    if (bin[i] >= 0) {
      int64_t k = bin[i];
      site.wgt = sums[3 * k] / sums[3 * k + 2];
      site.wgt_im = sums[3 * k + 1] / sums[3 * k + 2];
    }
    double w_abs = std::hypot(site.wgt, site.wgt_im);
    if (settings::noise_cancellation_roulette && w_abs < cutoff()) {
      uint64_t seed =
        init_seed(id_base + offset + i + 1, STREAM_URR_PTABLE);
      if (w_abs > 0.0 && survive() * prn(&seed) < w_abs) {
        double f = survive() / w_abs;
        site.wgt *= f;
        site.wgt_im *= f;
      } else {
        continue;
      }
    }
    bank[n_keep++] = site;
  }
  int64_t n_removed = bank.size() - n_keep;

  // Diagnostics. Bin totals are already global; local counts are reduced.
  double diag[4] = {abs_before, static_cast<double>(bank.size()),
    static_cast<double>(n_removed), abs_outside};
#ifdef OPENMC_MPI
  MPI_Allreduce(MPI_IN_PLACE, diag, 4, MPI_DOUBLE, MPI_SUM, mpi::intracomm);
#endif
  double abs_after = diag[3];
  double occupied = 0.0;
  for (int64_t k = 0; k < n_bins; ++k) {
    if (sums[3 * k + 2] > 0.0) {
      abs_after += std::hypot(sums[3 * k], sums[3 * k + 1]);
      occupied += 1.0;
    }
  }
  cancel_n_rendezvous += 1.0;
  cancel_abs_before += diag[0];
  cancel_n_sites += diag[1];
  cancel_n_removed += diag[2];
  cancel_abs_after += abs_after;
  cancel_n_occupied += occupied;

  bank.resize(n_keep);
  return n_removed;
}

} // namespace

void transport()
{
  thread_fission_banks.assign(num_threads(), {});
  simulation::shared_secondary_bank_read.clear();
  simulation::shared_secondary_bank_write.clear();

  if (perturbations.empty()) {
    // Generation 0: external source particles
    transport_tracks(simulation::work_per_rank, false);
  } else {
    // Generation 0: noise source sampled in the static phase
    int64_t n_src = sample_noise_source();
    calculate_work(n_src);
    transport_tracks(n_src, true);
    track_offset += n_src;
    simulation::shared_secondary_bank_read.clear();
  }

  vector<double> batch_sum;
  int gen = 0;
  double sum_first = 0.0;
  while (true) {
    collect_fission_banks();

    // Weight cancellation at the fission-generation rendezvous. The bank is in
    // canonical order here, so results do not depend on the thread count.
    if (settings::noise_cancellation_mesh != C_NONE) {
      cancel_and_roulette(simulation::shared_secondary_bank_write, track_offset);
    }

    int64_t n_alive =
      synchronize_global_secondary_bank(simulation::shared_secondary_bank_write);

    // Bank statistics (global)
    double stats[4] = {0.0, 0.0, 0.0, 0.0};
    const auto& bank = simulation::shared_secondary_bank_write;
    for (int64_t i = 0; i < bank.size(); ++i) {
      double a = std::hypot(bank[i].wgt, bank[i].wgt_im);
      stats[0] += a;
      stats[1] += a * a;
      stats[2] += bank[i].wgt;
      stats[3] += bank[i].wgt_im;
    }
#ifdef OPENMC_MPI
    MPI_Allreduce(MPI_IN_PLACE, stats, 4, MPI_DOUBLE, MPI_SUM, mpi::intracomm);
#endif
    if (n_alive == 0)
      break;
    ++gen;

    if (gen_sum_abs_wgt.size() < gen) {
      gen_sum_abs_wgt.push_back(0.0);
      gen_sum_abs_wgt_sq.push_back(0.0);
      gen_sum_wgt_re.push_back(0.0);
      gen_sum_wgt_im.push_back(0.0);
      gen_n_sites.push_back(0.0);
      gen_n_batches.push_back(0);
    }
    gen_sum_abs_wgt[gen - 1] += stats[0];
    gen_sum_abs_wgt_sq[gen - 1] += stats[1];
    gen_sum_wgt_re[gen - 1] += stats[2];
    gen_sum_wgt_im[gen - 1] += stats[3];
    gen_n_sites[gen - 1] += n_alive;
    gen_n_batches[gen - 1] += 1;
    batch_sum.push_back(stats[0]);

    if (mpi::master) {
      write_message(7, " Noise fission generation {:<6} sites: {:<10} "
                       "sum|w|: {:.6e}",
        gen, n_alive, stats[0]);
    }

    // Guards
    if (gen > settings::noise_max_generations) {
      fatal_error(fmt::format("Noise mode: batch {} exceeded the maximum "
                              "number of fission generations ({}).",
        simulation::current_batch, settings::noise_max_generations));
    }
    // Growth guard: in a convergent game the expected sum of |w| decays by
    // K~ < 1 per generation, so a sustained rise far above the first fission
    // generation signals a divergent (K~ > 1) game. A consecutive-increase
    // test on a moving average is not used because it fires by chance in
    // near-critical convergent games (bank size ~ critical random walk).
    if (gen == 1)
      sum_first = stats[0];
    if (stats[0] > settings::noise_growth_limit * sum_first) {
      fatal_error(fmt::format(
        "Noise mode: the sum of |w| over the fission bank grew by a factor "
        "{:.3g} over the first fission generation (limit {}; batch {}, "
        "generation {}). The fission-generation game is divergent at this "
        "frequency without weight cancellation; increase eta or enable "
        "cancellation.",
        stats[0] / sum_first, settings::noise_growth_limit,
        simulation::current_batch, gen));
    }

    // Transport the generation
    calculate_work(n_alive);
    simulation::shared_secondary_bank_read =
      std::move(simulation::shared_secondary_bank_write);
    simulation::shared_secondary_bank_write = SharedArray<SourceSite>();
    transport_tracks(simulation::shared_secondary_bank_read.size(), true);
    track_offset += n_alive;
    simulation::shared_secondary_bank_read.clear();
  }

  batch_generations.push_back(gen);
  calculate_work(settings::n_particles);
  weight_scale = 1.0;
}

void write_statepoint(hid_t file_id)
{
  hid_t g = create_group(file_id, "noise");
  write_dataset(g, "frequency", settings::noise_omega / (2.0 * PI));
  write_dataset(g, "eta", settings::noise_eta);
  write_dataset(g, "keff", settings::noise_keff);
  write_dataset(g, "generations_per_batch", batch_generations);
  write_dataset(g, "gen_sum_abs_wgt", gen_sum_abs_wgt);
  write_dataset(g, "gen_sum_abs_wgt_sq", gen_sum_abs_wgt_sq);
  write_dataset(g, "gen_sum_wgt_re", gen_sum_wgt_re);
  write_dataset(g, "gen_sum_wgt_im", gen_sum_wgt_im);
  write_dataset(g, "gen_n_sites", gen_n_sites);
  write_dataset(g, "gen_n_batches", gen_n_batches);
  if (!perturbations.empty()) {
    write_dataset(g, "nskip", nskip);
    write_dataset(g, "static_inactive", static_inactive);
    write_dataset(g, "static_k", static_k_gen);
    write_dataset(g, "source_n_sites", batch_n_source);
    write_dataset(g, "source_weight_scale", batch_scale);
    // Per batch and perturbation: (re, im) sum of sampled source weights
    write_dataset(g, "source_sum", batch_source_sum);
    write_dataset(g, "n_perturbations", static_cast<int>(perturbations.size()));
  }
  if (settings::noise_cancellation_mesh != C_NONE) {
    write_dataset(g, "cancel_n_rendezvous", cancel_n_rendezvous);
    write_dataset(g, "cancel_n_sites", cancel_n_sites);
    write_dataset(g, "cancel_n_occupied", cancel_n_occupied);
    write_dataset(g, "cancel_abs_before", cancel_abs_before);
    write_dataset(g, "cancel_abs_after", cancel_abs_after);
    write_dataset(g, "cancel_n_removed", cancel_n_removed);
  }
  close_group(g);
}

} // namespace noise
} // namespace openmc
