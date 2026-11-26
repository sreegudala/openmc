# Status

Updated: 2026-11-20.
Stage definitions and gates: design doc §5. Overview: [`README.md`](README.md).

## Stages

| Stage | What it adds | Status | Record |
|---|---|---|---|
| 0 | Preflight: reference scripts, Abeille build, library check | **Passed** (with notes below) | `stage0/stage0_log.md` |
| 1 | Complex weights, complex tallies, Re–Im covariance | **Passed** (35/35, exact) | `stage1/stage1_log.md` |
| 2 | $i\omega/v$ copy reaction, no fission | **Passed** (43/43) | `stage2/stage2_log_2a_2b_f7.md` |
| 3 | Fission with $1/k$, delayed neutrons, generation loop, growth guard | **Passed** | `stage3/stage3_gate.md` |
| 4 | Space dependence (slabs vs $S_N$) | **Passed** (9/9); STD-C criterion introduced | `stage4/stage4_gate.md` |
| 5 | Weight cancellation | **Passed**; 5d gated on Im FOM (amendment) | `stage5/stage5_gate.md` |
| 6 | Sampled noise source from $\delta\Sigma$ perturbations | **Passed**: 6a/6b 22/22, 6c 6/6 vs Abeille (global amplitude resolved to about ±5 %, shape to 0.3 %) | `stage6/stage6_gate.md` |
| 7 prep | Brief N1–N3: copied TRIGA transposition, noise extension, audits | **Done**: N1/N2 unit tests pass; N3a/N3b 6/6, and N3b at 4× particles plus the N3f forward control 6/6 (job 3599802; adjoint FOM 4–11 % of forward on this material). N3c (weight windows) moved to 9b prep | `stage7/stage7_prep_log.md` |
| 7 | Adjoint noise and reciprocity (amended 2026-11-20, 2026-11-20b) | **Passed** (30/30 after the 1 Hz rerun; 7c scattering checked against S_N, see gate) | `stage7/stage7_gate.md` |
| 9b prep | N3c: weight windows with complex weights; N5: shared TRIGA library (amendment 2026-11-22, runs before Stage 8) | **N3c passed** (unit tests; deep penetration 6/6, FOM ×240 and ×8600; infinite-medium adjoint 6/6, unbiased, no gain). N5 waits for the regenerated library | `prep9b/n3c_gate.md` |
| 8 | Vibrating component (amendments 2026-11-22b, 2026-11-23) | **Passed**: 8a 6/6, 8b 4/4, 8c 16/18 (the two n = 1 thermal misses against Abeille dispositioned as an error in the reference: the adjoint contraction agrees with OpenMC and is 7.6σ and 5.5σ from Abeille) | `stage8/stage8_gate.md` |
| 9 | Many-group benchmark (C5G7) and TRIGA dry run | **9a passed** (12/15 with disposition: 3 single reflector cells at 3.6–3.8σ; relative noise within about 1 % and 0.1° of FEMFFUSION SP3, whose default noise operator had to be replaced by its full one). 9b: waiting for the N5 library | `stage9/stage9a_gate.md` |

## Next

- Abeille n = 1 without cancellation (job 3612742, running): append its
  result to `stage8/stage8_gate.md`.
- Decide whether the per-cell tolerance rule should allow for the number
  of cells (stage9a_gate.md, end) before 9b; needs an amendment.
- Stage 9b (TRIGA): waiting for the N5 library.

## Open items and known limitations

- **Variance of transposed libraries.** A multiplicity above 1 can make the
  second moment infinite (`adjoint_noise_lib.second_moment_radius` ≥ 1; N3x
  has 2.05). Check this for every adjoint test material; the TRIGA library
  needs windows (N3c, implemented 2026-11-22).
- **OpenMC MG multiplicity (found in ce_mg_triga 540ab07).** MGMC samples
  the outgoing group from the nu-scatter row and then multiplies by M,
  so a multiplicity that varies along a row is counted twice
  (`ce_mg_triga/openmc_issues/ISSUE_mg_multiplicity_row.md`, not posted).
  N1 makes M constant along every row, so the adjoint libraries here are
  not affected, and `effective_operator` models the same sampling. Check
  any forward library with (n,2n), such as the N5 library above 6 MeV.
- **N5, shared TRIGA library.** Regenerate the CASMO-70 library once in
  ce_mg_triga (`step1d_mgxs.py`) with inverse velocity, delayed
  nu-fission and χ_d per family, β and tallied λ_k. Both projects then use
  it, transposed by N1 + N2. When it arrives:
  - check forward and transposed MG k against CE;
  - check the low-ω response against G₀ with the IFP β_eff and Λ_eff.

- **MPI.** The local build has no MPI.
  - Stage 1d's 1-vs-4-rank check has not been run.
  - The sampled noise source refuses more than one MPI process (fatal
    error).
  - Runs use threads plus concurrent replicas instead.
- **FEMFFUSION**: built locally (`~/dev/external/FEMFFUSION`, commit
  1611f90, Ubuntu deal.II 9.5.1 / PETSc 3.19; its 99 built-in tests pass,
  including the noise and vibration ones). The TACC build script is
  `tacc/build_femffusion.sh` (deal.II module). Its cross-checks in 4b, 6c
  and 9a are informational; 9a's gate depends on it.
  - On LS6 (`dealii/9.4.0-real`, Intel 19.1, PETSc/SLEPc 3.16) it
    passes its test suite after six build-time workarounds, listed in
    `tacc/README.md`. Two are real FEMFFUSION bugs that GCC hides: 26
    `void` PETSc callbacks, and integer cell indices computed through
    `log10`/`pow`. The ROM unit tests are skipped: that deal.II has no
    LAPACK.
- **Stage 8 geometry: resolved.** The paper's x-vibrating pin at
  (−2.52, −2.52) and Abeille's y-vibrating pin at (+2.52, −2.52) are related
  by a 90° rotation of a fully symmetric assembly, so the fields are
  rotations of each other. 8c uses Abeille's input. The Eq. 26 form is
  arccos, as in Abeille. Stage 8's gate adds the adjoint contraction
  ⟨χ, δS⟩ as a second estimator (amendment 2026-11-22b).
- **Abeille local patch.** The reference build carries one fix (kcol
  restored before `record_generation`); see `tacc/build_abeille.sh`.
- **Plots.** The design doc asks for Bode plots, z histograms and
  $\sum|w|$ curves per stage. They have not been produced yet; the logs
  hold the numbers.
- **Heavy-tailed statistics.** Twice (3b, 5a), 20 small replicas
  understated the spread. Gate runs now use at least 20000 particles per
  replica, and misses are rerun with fresh seeds before being treated as
  bias.

## Design doc amendments

| Date | Stage | Change |
|---|---|---|
| 2026-11-18 | 4 | STD-C for correlated mesh tallies: integral z plus Hotelling $T^2$ on principal components replaces the z-histogram test |
| 2026-11-19 | 5 | 5d gates the Im FOM ratio; mesh averaging cannot reduce Re variance when weights are in phase |
| 2026-11-20 | 7 | Bilinear convention (transposed library at +ω, conjugation only in post-processing); adjoint = forward solver on the N1 + N2 transposed library (no code change); tests 7a–7f redefined per the brief |
| 2026-11-20b | 7 | 7e dropped (negative frequencies are rejected and not needed under the bilinear convention); 7c uses P0 + P1 for both perturbations (S_N: P0 alone is 25 % off for a scattering change) and a paired test; penetration statistics deferred to 9b |
| 2026-11-22 | 9b prep | Stage 9b preparation (N3c weight windows, N5) runs before Stage 8. Noise-mode windows: absolute (source weight scale), modulus comparisons, phase-preserving split, fission sites checked at birth, windows replace the noise roulette inside their mesh; tests N3c-u/a/b |
| 2026-11-22b | 8 | Pin orientation resolved by symmetry (Abeille's used); gate adds the adjoint contraction ⟨χ, δS⟩ (P0 + P1, from tallies) as a second estimator; statistics sized from a measured forward σ |
| 2026-11-23 | 8 | Vibration coefficient c_n derived; sampling rule (each side's own kernel and fission data, normalised by the collision material's Σt); `vibration` input |
| 2026-11-26 | 9a | Reference uses FEMFFUSION's full noise operator (`-matrixfree_type_A full_allocated`): the default non_diagonal SP_N noise operator omits the cross-moment Marshak terms at vacuum boundaries (δφ/φ₀ fell to 0.75 at the vacuum edge in SP3, 0.60 in SP5, not in SP1) |
| 2026-11-25 | 9a | Reference: FEMFFUSION `ref1_FE2` SP3 of examples/2D_C5G7_NN as shipped (UO2 kinetics on the perturbed MOX mix reproduced); relative noise per pin cell and group; outside the perturbed assembly, amplitude 5 % and phase 2° where OpenMC is precise enough, ≥ half the cells evaluated, both cancellation meshes (51, 102); OpenMC k from a separate eigenvalue run, σ_k ≤ 3 pcm |

## Log

- 2026-11-26: **Stage 9a passed** (12/15 with disposition).
  - OpenMC k = 1.186489 ± 2.2 pcm (job 3615390); noise at cancellation 51
    and 102 (jobs 3615392/3).
  - Against FEMFFUSION's shipped noise operator every group failed, in the
    reflector at the vacuum boundary: that operator omits the cross-moment
    Marshak terms (amendment 2026-11-26). With its full operator: relative
    noise median −0.8 % / −1.0 %, core sums within 1 % and 0.1°, phases
    within 1.5°; 3 single reflector cells over 5 % at 3.6–3.8σ, dispositioned.
- 2026-11-25: Stage 9a set up (amendment 2026-11-25). `stage9/c5g7.py`
  reads FEMFFUSION's C5G7 noise example in place; OpenMC k = 1.18674 ±
  0.00024 (C5G7 transport reference 1.18655). FEMFFUSION reference run
  locally: ref1_FE2 SP3 k 1.18271; relative noise outside the perturbed
  assembly within 1.8 %/0.06° of ref0_FE1, 2.5 % of SP1, 0.5 % of SP5.
  A 30-batch OpenMC probe is within about 5 % and 3° of FEMFFUSION,
  before k is pinned down (at 1 Hz 1 pcm of k is about 0.1° and 0.2 %).
- 2026-11-25: FEMFFUSION built on LS6 (`tacc/build_femffusion.sh`): ALL
  TEST PASSED (ROM tests skipped). Workarounds in `tacc/README.md`;
  `tacc/diag_femffusion.sh` diagnoses start-up crashes.
- 2026-11-25: FEMFFUSION built locally (Ubuntu packages; 99/99 built-in
  tests pass). OpenMC is unaffected (a fresh configure still finds serial
  HDF5). TACC build script added.
- 2026-11-25: **Stage 8 gate passed.**
  - 8c gate runs (job 3612744, 960 batches): forward against contraction
    passes on all 6 detector/mesh pairs; the meshes are consistent; n = 2
    passes against Abeille.
  - n = 1 thermal against Abeille fails at both meshes, dispositioned as
    the reference (`stage8/diag/`).
- 2026-11-24: **8c near-pin discrepancy resolved on the OpenMC side**
  (`stage8/diag/README.md`).
  - n = 1 half-pin detectors: the adjoint contraction (no cancellation,
    independent of the forward vibration sampling) gives +0.740j ± 0.003
    above and −0.632j ± 0.003 below.
  - OpenMC forward agrees with it at no cancellation, 170 and 340 (within
    0.8σ). Abeille is 3.1× and 2.8× larger (7.6σ and 5.5σ off).
  - Abeille agrees with S_N in 1-D (slab, and a vibrating plate with two
    opposite interfaces), so its 2-D error is specific to this case; its
    cancellation is the remaining suspect.
  - Next: Abeille n = 1 with the shipped setting (cancellation off).
- 2026-11-24: **8a and 8b passed 10/10** (job 3609592, 100000 × 400 static
  particles, cancellation 1 cm).
  - 8a n = 1: detector within 0.05 % of S_N (σ 0.9 %); forward against
    contraction within 0.16 %; contraction within 0.21 % of S_N (σ 0.18 %).
  - 8a n = 2: forward σ 8.9 %, contraction σ 0.26 %, consistent with S_N.
  - 8b: response / ε at 0.05 and 0.1 cm is 1.005 ± 0.014 and
    0.990 ± 0.012 of the 0.2 cm value.
  - 8c diagnostic (job 3609594): OpenMC at cancellation 170 and 340
    agree with each other. Against Abeille the near-pin dipole has the
    same shape but about 3× the amplitude in Abeille (thermal upper half
    of the pin cell: 2.29 ± 0.20 against 0.95 ± 0.27, about 4σ).
  - Abeille on the 8a slab (local) agrees with S_N (0.973), and on a 1-D
    vibrating plate (two opposite interfaces, local) its dipole swing
    matches S_N (0.12 against 0.125). OpenMC on the same plate (10
    replicas, 0.125 cm cancellation) gives 0.1292 ± 0.0072 against S_N
    0.1262 and Abeille 0.1197, so all three agree in 1-D. The 8c factor of
    about 3 is specific to 2-D (lattice, finite layer width, or
    cancellation on a 170 × 170 mesh). Next: adjoint contraction for
    half-pin detectors (no cancellation anywhere), and an OpenMC forward
    run without cancellation.
- 2026-11-24: Stage 8 probes (jobs 3607700–3607705).
  - Slab: both harmonics agree with S_N.
  - 8c: the forward estimate is very noisy, for Abeille too; the adjoint
    contraction is about 700× more efficient for a detector.
  - Abeille n = 1 done (600 batches); n = 2 crashed.
- 2026-11-23: Stage 8 started.
  - Amendment 2026-11-23: c_n derived (arccos form; Abeille's conventions
    confirmed in its source, Δ = negative − positive). Each side's own
    kernel and fission data are sampled at collisions in the layer, which
    works where the material present has no fission.
  - Implemented the vibration source with unit tests: the source integral
    is −iπ/2 for n = 1 and 0 for n = 2.
  - `ref_sn1d` now accepts several materials (it compared materials
    elementwise).
  - Drivers for 8a/8b/probe and 8c are smoke-tested.
- 2026-11-23: **N3c passed** (jobs 3606525, 3606526).
  - Deep-penetration water slab: windowed detector within 0.11 % / 0.18 %
    of S_N; on vs off consistent; FOM ×240 at 1 kHz and ×8600 at 10 kHz.
  - Infinite-medium adjoint with group windows: unbiased (6/6), but the
    gap to forward is not closed. N3 gains 1–2× (adjoint/forward 0.07–0.25
    with windows against 0.05–0.30 without); CORTEX is 10–25× worse (rare
    thermal → fast transfer; windows cannot bias a reaction channel).
  - 10 kHz field check (logged): thermal bins near the source are starved
    by detector-optimised windows.
- 2026-11-22: N3c implemented (design doc amendment 2026-11-22).
  - Code: `apply_weight_window` (modulus, phase-preserving split, absolute
    windows); collision roulette hand-over; fission sites checked at
    birth.
  - Found: `settings.cpp` had refused windows since the first noise
    commit, so the N3c audit's "they apply" was wrong (correction noted
    in the prep log). Windows also forced OpenMC's shared secondary bank,
    which noise mode refuses.
  - Fixed during smoke tests: a segfault from the random-number stream
    not being set at a fission site's birth check.
  - Non-noise runs are bit-identical to the pre-change build (MG window
    case, 1 thread).
  - Unit tests: 27 noise tests pass, including complex-phase windows
    against the analytic answer and generators refused.
- 2026-11-22: N1 copy refreshed to ce_mg_triga 540ab07 (`transpose_library`
  plus its HDF5 wrapper, both verbatim). Sync tests: the whole-file fixture
  regenerated byte-identical to the e018a83 one, and the pure function
  matches TRIGA's own `sync_arrays.npz`. 13 unit tests pass.
- 2026-11-22: **Stage 7 gate passed.** Rerun of 7b/7d at 1 Hz (job
  3602980, 4× particles, fresh seeds): forward −0.29 % from S_N (z −1.0),
  reciprocity −0.43 % ± 0.5 %. The first-run +1.5 % was 20-replica
  spread, as in 3b and 5a.
- 2026-11-21: Stage 7 gate runs (commit c12f32cd8). 7a 10/10, pooled z
  mean +0.38, std 0.99; N3 amplitudes within 0.55 %, CORTEX within
  about 1 σ (1 %). Slab job 19/20:
  - 7b reciprocity passes at 1, 30 and 100 Hz (|fwd/adj| − 1 = +2.0 %,
    +1.1 %, +0.1 %, all within 2σ);
  - 7c absorber: fwd/contraction within 0.02–1.3 %;
  - low-ω check, 7f and 7d at 30 and 100 Hz pass;
  - adjoint FOM 0.3–0.6 of forward in the slab;
  - miss: 7d forward 1 Hz (see Next).

  N3b/N3f follow-up 6/6.
- 2026-11-20: Stage 7 driver written and smoke-tested locally (all paths
  run). New: `common/pt_lib.py`, a verbatim copy of ce_mg_triga's P0/P1
  contraction (540ab07) with a sync fixture, and
  `adjoint_noise_lib.pt_terms_complex`. S_N on the 4b slab: reciprocity
  exact to 1e-15; break 29 Hz (Λ_eff 29 µs, β_eff 535 pcm). Amendment
  2026-11-20b. Harness: an exact 0 = 0 with no variance now gives z = 0
  (it was NaN, i.e. a fail; no earlier gate hit this case).
- 2026-11-20: **Stage 6 gate passed.** 6c (CORTEX at 1 Hz) vs Abeille
  (job 3596959, 20 replicas, 219–230 batches each at the 21 h cap):
  two-sample STD-C passes in both groups (Hotelling p 0.63, max |z| 1.4).
  The OpenMC modulus is 1.1–1.5 % below Abeille's, 0.6 combined σ; the
  2 % check passes but only resolves about ±5 %, because Abeille's σ is
  2.2 %. The spatial pattern agrees to 0.3 %, the phase to 0.14 deg.
- 2026-11-20: Stage 7 prep gate run (TACC 3598276): N3a and N3b pass 6/6.
  N3a is precise: 0.15 % amplitude at 100 Hz. N3b's σ is 1–4 %, too large
  to test the 2 % tolerance. N3x (infinite variance) fails as predicted:
  means 8–10 % low, |z| up to 14. Abeille 6c still running (21 h cap).
- 2026-11-20: noise mode refuses weight windows (fatal error) until N3c
  is designed.
- 2026-11-20: Stage 7 brief N1 and N2 done. The TRIGA transposition was
  copied verbatim with a sync test, and the noise extension uses the matrix
  fission form. Unit tests pass: transposing twice is the identity, k is
  equal, and OpenMC transports exactly Aᵀ. N3 audits are written up: the
  noise game is complex-safe under multiplicity, weight windows are not, and
  a transposed library can have infinite variance. The Stage 7 amendment is
  in the design doc. The N3a/N3b gate run goes to TACC.
- 2026-11-19: Stage 6a/6b results committed (22/22). 6c OpenMC side done;
  Abeille built on TACC (node-local build) and running.
- 2026-11-19: Stage 5 gate passed.
- 2026-11-18: Stage 4 gate passed (STD-C amendment).
- Stage 3 gate passed; Stages 0–2 passed (initial commit 9e4eb3096).
