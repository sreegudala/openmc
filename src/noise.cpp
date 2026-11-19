#include "openmc/noise.h"

#include <algorithm>
#include <cmath>

#include <fmt/core.h>

#include "openmc/bank.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
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
#include "openmc/vector.h"

namespace openmc {
namespace noise {

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

} // namespace

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
}

void initialize()
{
  reset();

  if (settings::event_based) {
    warning("Event-based transport is not available in noise mode; "
            "history-based transport will be used.");
    settings::event_based = false;
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
}

void roulette(Particle& p)
{
  if (p.wgt_abs() < settings::weight_cutoff) {
    russian_roulette(p, settings::weight_survive);
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
    if (w_abs < settings::weight_cutoff) {
      if (settings::weight_survive * prn(p.current_seed()) < w_abs) {
        double f = settings::weight_survive / w_abs;
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

void sample_reaction(Particle& p)
{
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

//! Gather thread fission banks into shared_secondary_bank_write in canonical
//! order (parent global index, then progeny index)
void collect_fission_banks()
{
  int64_t n = 0;
  for (const auto& b : thread_fission_banks)
    n += b.size();
  vector<SourceSite> all;
  all.reserve(n);
  for (auto& b : thread_fission_banks) {
    all.insert(all.end(), b.begin(), b.end());
    b.clear();
  }
  std::sort(all.begin(), all.end(), [](const SourceSite& a, const SourceSite& b) {
    return a.parent_id < b.parent_id ||
           (a.parent_id == b.parent_id && a.progeny_id < b.progeny_id);
  });
  simulation::shared_secondary_bank_write.resize(0);
  simulation::shared_secondary_bank_write.extend_uninitialized(n);
  for (int64_t i = 0; i < n; ++i)
    simulation::shared_secondary_bank_write[i] = all[i];
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
    if (settings::noise_cancellation_roulette &&
        w_abs < settings::weight_cutoff) {
      uint64_t seed =
        init_seed(id_base + offset + i + 1, STREAM_URR_PTABLE);
      if (w_abs > 0.0 && settings::weight_survive * prn(&seed) < w_abs) {
        double f = settings::weight_survive / w_abs;
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
  double local[4] = {abs_before, static_cast<double>(bank.size()),
    static_cast<double>(n_removed), abs_outside};
#ifdef OPENMC_MPI
  MPI_Allreduce(MPI_IN_PLACE, local, 4, MPI_DOUBLE, MPI_SUM, mpi::intracomm);
#endif
  double abs_after = local[3];
  double occupied = 0.0;
  for (int64_t k = 0; k < n_bins; ++k) {
    if (sums[3 * k + 2] > 0.0) {
      abs_after += std::hypot(sums[3 * k], sums[3 * k + 1]);
      occupied += 1.0;
    }
  }
  cancel_n_rendezvous += 1.0;
  cancel_abs_before += local[0];
  cancel_n_sites += local[1];
  cancel_n_removed += local[2];
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

  // Generation 0: external source particles
  transport_tracks(simulation::work_per_rank, false);

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
