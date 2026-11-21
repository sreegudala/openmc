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
| 7 | Adjoint noise and reciprocity (amended 2026-11-20, 2026-11-20b) | **In progress**: 7a 10/10 (job 3600163); slab job 19/20 (job 3600164). The miss is 7d forward at 1 Hz (+1.5 %, z 3.5), being rerun with fresh seeds at 4× particles | `stage7/stage7_log_7a.md`, `stage7/stage7_log_7b_lf_7c_7f.md` |
| 8 | Vibrating component | Not started | |
| 9 | Many-group benchmark (C5G7) and TRIGA dry run | Not started (blocked, see below) | |

## Next: Stage 7 rerun, then the gate record

- **7d forward 1 Hz.** The forward field is uniformly 1.5 % above S_N
  (integral z +3.5, Hotelling p 0.010, shape fine). The adjoint at 1 Hz
  agrees with S_N (−0.5 ± 1.0 %), and the reciprocity test passes
  (t 1.88).
  - Not a heavy tail: the replicas lie within 0.97–1.05.
  - Not the reference: S_N k is converged to 3e-7, and 1.5 % would need
    a k error of 8e-5 (0.19 % per 1e-5 at 1 Hz).
  - Rerun with fresh seeds at 4× particles:
    `run_stage7.py --only 7b --freqs-7b 1 --particles 80000 --base-seed
    100001 --tag p80k`.
- **7c scattering case.** The forward sampled source has 40–150 % σ: the
  copy (−δΣ_t ψ₀) and the sampled scatter (+δΣ_s) cancel on average, so
  the paired test passes without power. The contraction matches S_N
  within 1σ (0.8–1.3 %) at all three frequencies, with P1 = 27 % of the
  answer; a wrong P1 sign or factor would be 40–65σ off. The gate record
  will state this. Cancelling sources of this kind return in Stage 8
  (vibration).
- Then `stage7/stage7_gate.md` and Stage 8.

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
- **Refresh the N1 copy.** ce_mg_triga 540ab07 (2026-11-20) split
  `transpose_library` out of `make_adjoint_library` as a pure function
  for this branch (`CHANGELOG_transpose.md`, its own sync fixture). The
  arithmetic is unchanged, so `common/adjoint_lib.py` (copied at e018a83)
  still matches. Refreshing it means moving the HDF5 I/O into
  `adjoint_noise_lib.py`.
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
