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
| 6 | Sampled noise source from $\delta\Sigma$ perturbations | **In progress**: 6a/6b pass (22/22), 6c waiting for Abeille | `stage6/stage6_log_6a_6b_nskip.md` |
| 7 prep | Brief N1–N3: copied TRIGA transposition, noise extension, audits | **In progress**: N1/N2 done (unit tests pass); N3a/N3b gate run pending on TACC; N3c audit: weight windows not complex-safe | `stage7/stage7_prep_log.md` |
| 7 | Adjoint noise and reciprocity (amended 2026-11-20) | Not started | |
| 8 | Vibrating component | Not started | |
| 9 | Many-group benchmark (C5G7) and TRIGA dry run | Not started (blocked, see below) | |

## Current work: Stage 6c (CORTEX assembly at 1 Hz vs Abeille)

- **OpenMC:** done. TACC job 3596962, 20 replicas × 600 batches × 20000
  particles, cancellation mesh 68×68. Output in `stage6/stage6c_openmc.npz`
  (not committed, 18 MB).
  - The static flux shape matches Abeille's to a median 0.06% (worst cell
    0.9%).
  - The peak noise response is in the perturbed pin, so the orientation is
    right.
  - Per-cell σ of the mean is 1.2%, about the same as the assembly
    integral's. The field is dominated by one global (point-kinetic) shape.
- **Abeille:** running. TACC job 3596959 started 2026-11-19 15:16 CST, 20
  replicas × 600 batches × 1e5 particles, about 10–11 h per replica,
  21 h cap. The Abeille static normalisation is committed
  (`stage6/stage6c_abeille_static.npz`).
- **Next:** rsync `stage6c_abeille.npz` → `python stage6/run_stage6c.py analyze`
  → `stage6/stage6_gate.md` → commit.
  - The gate is a two-sample STD-C per group on cells above 10% of max,
    plus a 2% modulus check on the assembly integral and the perturbed pin.
  - Each code is about 1.2% uncertain, so the 2% check is only about 1.5σ
    wide and a miss may need more batches before calling it a bias.

## Next: Stage 7 preparation (brief `briefs/noise_stage7_brief.md`)

- TACC job: `stage7/run_stage7_prep.py` (N3a and N3b against
  `ref_inf.py`'s adjoint, plus the informational N3x).
- Approval needed: should noise mode refuse weight windows (fatal error)
  until N3c is implemented?

## Open items and known limitations

- **N3c, weight windows in noise mode.** They apply but are wrong for
  complex weights: Re-only comparisons, the cutoff kills every Re w < 0, and
  split copies drop Im (`stage7/stage7_prep_log.md`). A fix needs a
  design-doc amendment first. Stage 9b needs windows to reach the
  out-of-core detectors.
- **Variance of transposed libraries.** A multiplicity above 1 can make the
  second moment infinite (`adjoint_noise_lib.second_moment_radius` ≥ 1; N3x
  has 2.05). Check this for every adjoint test material; the TRIGA library
  will need windows (N3c).
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

## Log

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
