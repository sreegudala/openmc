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
| 7 prep | Brief N1–N3: copied TRIGA transposition, noise extension, audits | **Done**: N1/N2 unit tests pass; N3a/N3b 6/6, and N3b at 4× particles plus the N3f forward control 6/6 (job 3599802; adjoint FOM 4–11 % of forward on this material). N3c: weight windows are not complex-safe and noise mode refuses them | `stage7/stage7_prep_log.md` |
| 7 | Adjoint noise and reciprocity (amended 2026-11-20, 2026-11-20b) | **Passed** (30/30 after the 1 Hz rerun; 7c scattering checked against S_N, see gate) | `stage7/stage7_gate.md` |
| 8 | Vibrating component | Not started | |
| 9 | Many-group benchmark (C5G7) and TRIGA dry run | Not started (blocked, see below) | |

## Next: Stage 8 (vibrating component)

- First settle the open geometry question below: the vibrating-pin
  orientation differs between the paper summary (x) and Abeille's input
  (y).
- The vibration source has the same cancelling structure as the 7c
  scattering case: two large sampled terms of opposite sign. Expect a
  noisy forward estimate and plan the statistics for it.

## Open items and known limitations

- **N3c, weight windows in noise mode.** They apply but are wrong for
  complex weights: Re-only comparisons, the cutoff kills every Re w < 0, and
  split copies drop Im (`stage7/stage7_prep_log.md`). Noise mode now
  refuses them with a fatal error. A fix needs a design-doc amendment
  first. Stage 9b needs windows to reach the
  out-of-core detectors.
- **Variance of transposed libraries.** A multiplicity above 1 can make the
  second moment infinite (`adjoint_noise_lib.second_moment_radius` ≥ 1; N3x
  has 2.05). Check this for every adjoint test material; the TRIGA library
  will need windows (N3c).
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
- **FEMFFUSION** is not built (needs deal.II, PETSc and SLEPc). Its
  cross-checks in 4b, 6c and 9a are informational; 9a's gate depends on it.
- **Stage 8 geometry.** The vibrating-pin orientation differs between the
  paper summary (x) and Abeille's input (y); settle this before Stage 8.
  The Eq. 26 form is resolved: arccos, as in Abeille.
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

## Log

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
