# Stage 7 preparation log (brief `briefs/noise_stage7_brief.md`)

Date: 2026-11-20. Branch code at commit 255bb94f5 plus this commit; line
numbers refer to that tree.

## N1. Copy of the ce_mg_triga transposition

- `common/adjoint_lib.py` is a verbatim copy of:
  - `step3a_adjoint.py` lines 68–171 (`_nbins`, `_dense`, `_write_sparse`,
    `make_adjoint_library`);
  - `common.py` lines 415–450 (`FISSION_DATASETS`, `remove_fission`, which
    `make_adjoint_library` calls).
- Source: ce_mg_triga commit e018a837d4859a68e1ec6ab9417f09b667941727
  (2026-11-19). Only the imports are new; the header records the provenance
  and the "do not edit here" rule.
- The code matches the brief's description in every point:
  - nu-scatter is transposed and keeps the μ histogram;
  - Σa' = max(0, Σt − in-scatter);
  - the row multiplicity is in-scatter / (Σt − Σa');
  - Σt is unchanged;
  - all fission data are removed.

  One addition: it also checks the forward balance and refuses a library
  that is off by more than 1e-2.
- **Sync fixture.** The TRIGA repository does not ship a fixture yet.
  `common/fixtures/adjoint_lib/make_fixture.py` makes one by importing
  `make_adjoint_library` from the TRIGA checkout itself. The commit used is
  recorded in `fixtures/adjoint_lib/SOURCE`.
  - `forward.h5`: a random 4-group library with 8-bin histogram scattering,
    upscatter, multiplicity 1.3 on two transfers, and fission with chi + beta
    and two families.
  - `adjoint_step3a.h5`: its transpose. TRIGA's report: max multiplicity
    1.97, 3 groups with absorption zeroed.
- `common/test_adjoint_lib.py` checks that the copy reproduces the fixture
  exactly (every dataset and attribute). It **passes**.
- When the copy is refreshed, rerun `make_fixture.py --triga <checkout>`.
  If the TRIGA repository adds its own fixture, use that one instead.

## N2. Extension to the noise operator (`common/adjoint_noise_lib.py`)

**Fission form: matrix.** No C++ changes are needed.
- With no `chi` dataset, OpenMC's MG reader takes `prompt-nu-fission`
  [g_in, g_out] and `delayed-nu-fission` [k, g_in, g_out] as production
  matrices (`XsData::fission_matrix_no_beta_from_hdf5`, src/xsdata.cpp:376).
  It stores the family production per g_in and a normalised χ row per g_in.
- `Mgxs::sample_fission_energy` (src/mgxs.cpp:525) chooses prompt or family
  k from the g_in row sums and the outgoing group from that row.
- `noise::create_fission_sites` (src/noise.cpp:271) reads only the total
  production `nu_fission(g_in)` and `sample_fission_energy`. The transposed
  row sums are therefore exactly the adjoint production.
- Without delayed families, a single `nu-fission` matrix is written.

**How the forward data are read.** `production_matrices` rebuilds the
forward matrices the way OpenMC's readers do. It handles all six
forward forms:
- chi + nu-fission, without or with beta (beta per family or per family and
  group);
- chi-prompt/chi-delayed + prompt/delayed vectors;
- matrix forms with and without beta.

χ is normalised as OpenMC normalises it.

**Keys in the transposed library.**
- `prompt-nu-fission` and `delayed-nu-fission`: the transposes.
- `decay-rate`: copied.
- `inverse-velocity` and `total`: unchanged (kept by N1).
- `absorption`: N1's balance value. It is not a reaction rate and contains
  no fission, which is correct for the noise collision:
  - the continuing weight is (Σt − Σa')/Σ̃t × c_g = in-scatter/Σ̃t;
  - fission is a separate branch with expected nu_fission/(Σ̃t k) sites.
- `fission`, `kappa-fission`, `chi` and `beta` are absent. Scores other
  than flux are meaningless in adjoint runs.

**Reads of non-nu `fission` in the noise path: none.**
- The noise path uses `macro_xs().nu_fission`, `.total`, `.absorption` and
  `.noise`.
- `fission` and `kappa-fission` are read only for tally scores.

**Zero-production rows.** Rows of the transposed matrices whose production
is zero (groups with χ = 0) give a 0/0 = NaN χ row in OpenMC.
- Such rows are never sampled, because their production, and hence the
  number of sites, is zero.
- `XsData::combine` only scales them.
- They would poison a mixture of several XSdata in one material; our
  libraries are macroscopic, one XSdata per material.
- The N3b runs (thermal χ = 0) confirm this in practice.

**Checks** (`common/test_adjoint_noise_lib.py`, 9 cases, all **pass**). Each
check runs on three forward libraries: the fixture (chi + beta), the N3
material and the N3x material (prompt/delayed vectors).
- Transposing twice returns the original fission and delayed matrices,
  decay rates, nu-scatter, Σt and 1/v (rtol 1e-14).
- The infinite-medium k of the transposed library equals the forward k
  (rel 1e-12). Example, N3 material: 1.577161055718 for both.
- What OpenMC actually transports with the transposed library is exactly Aᵀ:
  - scatter including the rebalance and multiplicity, where absorption is
    zeroed;
  - prompt and delayed production;
  - 1/v and λ.

## N3a. Complex weight under multiplicity: audit

Every place that scales the weight. ✓ means both parts are scaled by the
same real factor (or the complex product is correct). ✗ means wrong for
complex weights.

| Place | Code | Status |
|---|---|---|
| Scatter multiplicity | src/physics_mg.cpp:90–101 `scatter()`: complex branch passes `mult = 1` to `sample_scatter`, then `p.scale_wgt(mult)` (line 100); real branch (wgt_im = 0) multiplies `wgt` in place | ✓ |
| Implicit capture with rebalanced Σa' | src/noise.cpp:463–464 `p.scale_wgt(sigma_s / sigma_t)`, σs = Σ̃t − Σω − Σa' ≥ 0 since Σa' ≤ Σt | ✓ |
| Survival biasing | Noise mode forces implicit capture in `noise::sample_reaction`; `physics_mg.cpp absorption()` (which also handles wgt_im, line 268–273) is not called | ✓ |
| Copy reaction | src/noise.cpp:456–459 `create_secondary(wr + wi/η, …, wi − wr/η)` = w(1 − i/η) | ✓ |
| Roulette (collision) | src/noise.cpp:264–269 → src/physics_common.cpp:12–35 `russian_roulette`: norm-based, `scale_wgt(weight_survive/|w|)` (line 29) | ✓ |
| Roulette (fission sites) | src/noise.cpp:326–332, both parts | ✓ |
| Delayed factor | src/noise.cpp:312–321, complex multiply | ✓ |
| Fission site weight | src/noise.cpp:293–294 copies both parts | ✓ |
| Cancellation and bank roulette | src/noise.cpp:703–765, averaging (752–753) and roulette (755–763) on both parts | ✓ |
| Kill | include/openmc/particle_data.h:665–669 zeroes both parts; every `wgt() = 0` on the branch is `kill()` | ✓ |
| Sampled source (static phase) | src/noise.cpp:355 `w = p.wgt()`: static weights are real (ω = 0, no copies, delayed factor 1); scatter term includes `mult` (line 419) | ✓ |
| Tally scoring | src/tallies/tally_scoring.cpp:2533–2562 uses wgt_im / wgt_im_last | ✓ |
| Splitting | Only in weight windows, see N3c | ✗ |
| `wgt_born` | src/noise.cpp:494 sets `wgt_born = Re w`; used only by `survival_normalization` and weight windows, neither of which is used by the noise game | (✗ under N3c) |

Conclusion: the noise game is complex-safe under multiplicity. The
transposed |w| game also has the same convergence as the forward one. Its
expected-|w| transfer matrix is D⁻¹Kᵀ against the forward D⁻¹K (D =
diag Σ̃t), and the two are similar matrices.

**Second moment (new finding).** Multiplicity changes the *second* moment.
`adjoint_noise_lib.second_moment_radius` is the spectral radius of the
per-collision second-moment operator Q = p f², with p the outgoing-group
probabilities and f the weight factor per collision.
- Forward libraries (f ≤ 1) always have a radius below 1.
- A transposed group with c = 2.4 and self-transfer probability 0.35 gives
  2.4² × 0.35 = 2.0. The variance is then infinite in an infinite medium.
- The first N3 candidate material did exactly this (kept as **N3x**, not
  gated): radius 0.77 forward, **2.05 transposed**. With 3 replicas × 5
  batches × 2000 particles at 100 Hz, fission off, the two estimators
  disagree (thermal track length 11.4 vs collision 8.6, reference 9.8) and
  |z| reaches 15.
- The gated N3 material is designed to have a radius of 0.33 at the same
  multiplicity, 2.425.
- This is the mechanism behind the TRIGA adjoints' weight-tail bias, and the
  reason Stage 9b needs N3c.

**Test (N3a).** `stage7/run_stage7_prep.py --only n3a`. The adjoint source
is a thermal unit source; the reference is `ref_inf.solve(..., adjoint=True)`
with fission off.
- 100 Hz (plateau) and 3 kHz (break) at η = 1.
- 10 kHz (above) at η = 3. At η = 1 the thermal copy loop has gain
  x/(1+x)·√2 = 1.10 > 1 (x = 3.6): the Stage 2 collision-based cut-off,
  identical for forward and adjoint. The smoke run at η = 1 diverged
  (means ~ −4000), and at η = 3 it passes.

Gate run: TACC (pending), 20 × 20000 × 20 as in Stage 3.

## N3b. Fission data in the N2 form

- **Form: matrix** (above). The prompt and delayed noise fission paths
  handle it with no change: `create_fission_sites` and `sample_source_sites`
  both call `sample_fission_energy`, and the delayed factor uses the sampled
  family's `decay-rate`.
- **Test (N3b).** `stage7/run_stage7_prep.py --only n3b`.
  - The N3 material with fission and the six TRIGA families (β = 715 pcm,
    Keepin λ); Λ = 37 µs.
  - noise keff = k_inf = 1.577161 (forward and transposed agree to 1e-15).
  - 1 Hz (plateau) and 30 Hz (break) at η = 1, 100 Hz (above) at η = 3,
    where K̃cb = 0.9966 (convergent; the same regime as Stage 3b at 100 Hz).
  - Smoke run (3 × 5 × 2000): 1 and 30 Hz pass. At 100 Hz the Re z reaches
    3.3, from a heavy-tailed case with only 3 replicas (Stage 3 showed small
    replicas understate the spread); the gate run decides.

Gate run: TACC (pending).

## N3c. Weight windows in noise mode: audit

**Do they apply?** Yes, but incorrectly. Nothing in noise mode disables them.
- After `noise::sample_reaction`, `collision_mg` (src/physics_mg.cpp:36–45)
  calls `apply_weight_window`.
- Surface crossings call `apply_weight_windows` (src/particle.cpp:391).

**What is wrong for complex weights.** All of it is in
`apply_weight_window`, src/weight_windows.cpp:971–1053.
1. Every decision uses the real part only: `double weight = p.wgt()`
   (line 991).
2. The weight cutoff `if (p.wgt() < weight_cutoff) kill` (line 994) kills
   **every particle with Re w ≤ cutoff, including all negative Re weights**.
   Negative Re weights are routine: copies carry w(1 − i/η), and above the
   break Re δφ itself changes sign.
3. Splitting calls `p.split(weight / n_split)` (line 1041) with the default
   `wgt_im = 0`, so split copies lose the imaginary part. The surviving
   particle gets `p.wgt() = weight / n_split` (line 1044) and keeps its full
   Im part. The phase is not preserved.
4. Roulette goes through the norm-based `russian_roulette` (line 1051, ✓),
   but `weight_survive` is computed from Re w.
5. `ww_factor`/`max_lb_ratio` (lines 1002–1010) compare Re w.
6. Window normalisation uses `wgt_born / wgt_ww_born` (line 988), and the
   noise loop sets `wgt_born = Re w` (src/noise.cpp:494). For a sampled
   source, Re w can be ~0 or negative.
7. Outside the windows' mesh, `apply_russian_roulette` uses the unscaled
   `settings::weight_cutoff` (src/physics_mg.cpp:41). For a sampled source
   the noise loop instead scales its thresholds by `weight_scale` (the mean
   source modulus, ~1e-3 to 1e-2), so most particles would be rouletted
   there.

**Interaction with cancellation.** Cancellation acts on the fission bank at
generation boundaries (src/noise.cpp:703); windows act within a generation.
- They do not conflict directly.
- But fission sites are born with the parent's weight and are not
  checked against windows until their first collision or surface. A
  windowed game should also check sites at birth. This is a design choice
  for the amendment.

**What is missing** (to be done only after a design-doc amendment, as the
brief requires):
- Use |w| in every comparison (cutoff, upper/lower bound, max_lb_ratio,
  survival weight).
- Split with `p.split(w_re / n, w_im / n)` and scale the survivor with
  `scale_wgt(1/n)`.
- Define `wgt_born` as |w| at birth in the noise loop.
- Scale window bounds by the noise `weight_scale` for sampled sources, or
  normalise the source to unit mean modulus.
- Decide whether fission sites are checked at birth.

**Test.** Not run: windows cannot be tested before this is fixed.

**Guard (2026-11-20, approved by the user).** `noise::initialize`
(src/noise.cpp) now stops with a fatal error when weight windows are on in
noise mode, so they can no longer silently give wrong answers. Regression
test: `tests/unit_tests/test_noise.py::test_noise_refuses_weight_windows`.

## Gate-run results

TACC LS6 job 3598276, commit dc5189570, 20 replicas × 20000 particles × 20
batches, seeds 1 + 7919 r. Log: `stage7_prep_log_n3a_n3b_n3x.md`.
**6/6 gated checks pass** (STD, both estimators). Pooled z over the gated
cases: n = 72, mean +0.25, std 1.06, no |z| ≥ 3.

| Case | max \|z\| (tl, col) | rel σ of mean (thermal) | amplitude check | notes |
|---|---|---|---|---|
| N3a 100 Hz | 1.86, 1.66 | 0.02 % | 0.15 % (3 bins evaluated) | phase 0.01 deg |
| N3a 3 kHz | 1.14, 1.73 | 0.9 % | 0.6 % (2 bins) | phase ≤ 0.47 deg |
| N3a 10 kHz (η 3) | 1.21, 1.43 | 4.4 % | not evaluated (σ too large) | |
| N3b 1 Hz | 2.09, 2.09 | 1.1 % | not evaluated | 725 gens/batch |
| N3b 30 Hz | 1.34, 1.32 | 1.8 % | not evaluated | 1718 gens/batch |
| N3b 100 Hz (η 3) | 0.63, 0.62 | 3.7 % | not evaluated | 5055 gens/batch (max 19358) |
| N3x 100 Hz (not gated) | 13.95, 8.58 | — | — | infinite variance, as predicted |

**N3a** passes, with high precision at the plateau and the break: the
copied transposition with multiplicity 2.4 transports the adjoint
correctly, complex weights included.

**N3b** passes STD but with **low power**.
- The replica distributions are clean: no outliers, largest deviation ~2σ
  of the replica spread.
- The relative σ of the mean is 1.1–3.7 %, about 4–10 times that of the
  forward Stage 3b at the same particle count. The design doc's 2 %
  amplitude tolerance was therefore never evaluated.
- The z-test alone excludes amplitude biases above about 3 %, 5 % and 11 %
  at 1, 30 and 100 Hz.
- Follow-up before 7a (which promotes N3b to a gate):
  - N3f, the forward run on the same material, tells whether the adjoint
    game is intrinsically noisier or the material is just expensive. An
    adjoint fast-group collision yields νΣf′/(Σ̃t k) ≈ 2.5 sites, and 100 Hz
    needs 2.3 times the generations of Stage 3b.
  - N3b at 4× particles with fresh seeds.

**N3x** behaves as the second-moment radius (2.05) predicts.
- 5–6 of 6 bins fail, with max |z| 14 (track length) and 8.6 (collision).
- The means are 8–10 % low in every group, the direction of the TRIGA
  weight-tail bias: rare high-weight histories carry the missing mean.
