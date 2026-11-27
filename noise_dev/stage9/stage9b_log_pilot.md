# Stage 9b pilot (local, 2026-11-26)

Smoke runs of `run_stage9b.py` on the local machine (14 threads), before
the LS6 gate runs. Sizes are far below the gate's; the numbers show what
works and what does not, not results.

## Forward (mesh A, 1 Hz, 4 replicas × 4 batches × 5000, `--wwmap`)

- **The pipeline works.** δF/F₀ = −0.0496 + 0.0060i (N5's G₀ fit predicts
  −0.0470 at this frequency). Each batch runs about 730 fission generations.
- **Linearity.** One seed at ε = 0.01, 0.1 and 1 gives responses in the
  ratio of ε to 2e-13 (`stage9b_log_linearity.md`). The strength is
  therefore a scale factor, not a test.
- **The static d₀ at the far detectors is barely sampled.** It is zero in
  some replicas at det3–5 and arc5, while the noise δd is sampled at every
  detector (σ 4–30 % at this size). The two-mesh gate therefore compares
  δd, not δd/d₀.
- **Static k** per static generation scatters by 1.8 % at 5000 particles.
  The replica mean is 0.9985 ± 0.003 against `noise.keff` 1.00137. At the
  gate size it reaches about 20 pcm, so it stays logged.

## Solver bug found: weight windows at birth, several threads

- **Symptom.** The first adjoint smoke run segfaulted in `prn`, called from
  `russian_roulette`, from `apply_weight_window`, from
  `initialize_particle_track` (gdb).
- **Cause.** Noise mode's absolute windows can roulette a source particle at
  birth, before its first step sets the random-number stream.
  `ParticleData::stream_` is not initialised, so on a fresh thread-local
  `Particle` the seed pointer is garbage. Non-noise windows put every source
  particle at its window centre and never draw there.
- **Fix.** Set `STREAM_TRACKING` before the birth check (commit 3d32d7862),
  as the fission-site branch of `noise::transport_tracks` already did.
- **Test.** `test_noise_weight_windows_birth_roulette_threads`, 4 threads,
  fails without the fix. All 35 noise unit tests and 36 weight-window unit
  tests pass.

## Adjoint, arc0, 1 Hz: the first attempt fails

The adjoint used 4 replicas × 4 batches × 5000, without cancellation.
Contraction against the forward δd at arc0, which is −5.60e-4 + 7.08e-5i
(Im/Re −0.126, matching G₀'s phase):

| windows | contraction | Im/Re | fwd / contraction | generations per batch |
|---|---|---|---|---|
| CADIS from the forward MC φ₀ (8 groups) | −3.25e-4 + 6.5e-6i | −0.020 | 1.74 | 201 |
| flat 0.5–2 | −8.4e-5 | −0.0003 | 6.7 | 52 |
| none | −4.1e-7 | 0 | 1400 | 5 |

The contraction itself is converged in its discretisation. The P1 share is
5e-4, and 2 cm voxels change it by 4e-4, 4 cm by 3e-3.

## Diagnosis

- **The library is right.** With the delayed factors at 1 Hz, the
  transposed library's infinite-medium noise operator equals the transpose
  of the forward one to 4e-15 in three fuels, water and graphite. The
  exception is fuel_10147_2358, at 1.3e-5: its balance gap, absorbed by N1.
  The prompt/delayed split per family and the λ_k are identical, and
  deterministic reciprocity holds to 1e-12.
- **Infinite-medium test.** fuel_10146_2337 alone in a reflective box,
  `noise.keff` = k∞ = 1.52943, at 1 Hz, with 2000 particles × 10 batches per
  replica. Volume-integrated flux per source particle against the exact
  solve, 16 replicas:

  | run | Re | Im | exact |
  |---|---|---|---|
  | forward, no windows | 4913 ± 153 | −519 ± 32 | 4670 − 438i |
  | forward, roulette-only window | 4626 ± 117 | −474 ± 65 | |
  | forward, flat window | 4656 ± 171 | −365 ± 32 | |
  | adjoint, no windows | (4 replicas) 11.3, 3 % of exact | ~0 | 390 − 37i |
  | adjoint, flat window | 291 ± 88 | −3.5 ± 1.3 | |
  | adjoint, no windows + cancellation | 12.5 ± 2 | ~0 | |
  | adjoint, flat + cancellation | 432 ± 111 | −29 ± 13 | |
  | adjoint, group importance windows (70 groups) | 451 ± 66 | −41 ± 11 | |
  | adjoint, group windows (70) + cancellation | 388 ± 51 | −35 ± 8 | |
  | adjoint, group windows (CASMO-8) + cancellation | 387 ± 95 | −33 ± 13 | |

  Group importance is ψ = |A_f⁻¹·1|, lower bound C/ψ_g, ratio 5.
  Cancellation is one mesh cell, per group.
- **What the table shows.**
  - The forward game is unbiased with or without windows. Its spread is
    within about 2.5σ.
  - The adjoint without windows is the infinite-second-moment case of the
    transposed library (radius 1.9 in fuel; N5). Typical sample means are
    far too low, and the replica σ means nothing. N5's transposed k without
    windows (0.586) was the same thing.
  - Windows without importance in energy (flat) leave the adjoint heavy-
    tailed: its Im is wrong by 30× at a quoted 1.3σ.
  - Energy-resolved importance windows give an unbiased adjoint, and 8
    groups are enough (FOM 0.60 against 0.40 for 70 groups).
  - Cancellation further stabilises it. The adjoint FOM is about 1/11 of the
    forward's here, as in Stage 7.
- **In the TRIGA run** the forward MC map is unscored or noisy exactly where
  the adjoint particles start (the pool). The arc0 failure above combines
  that with no cancellation.

## Change

- **Importance:** ce_mg_triga's random-ray forward flux
  (`results/step3/forward_flux.npz`, read in place). It is CASMO-8, on the
  same 5 cm mesh as the windows, and covers the pool to r = 81 cm.
- **Adjoint cancellation:** mesh A, as the forward runs.

## Adjoint, arc0, 1 Hz, with the change

Same size: 4 replicas × 4 batches × 5000.

| windows | cancellation | contraction | σ | Im/Re | fwd / contraction | generations per batch |
|---|---|---|---|---|---|---|
| CADIS, random-ray map | A | −4.54e-4 + 2.7e-5i | 19 % | −0.059 | 1.24 | 240 |

- **The magnitude agrees** within about 1.3σ.
- **Im is still low** by about 5σ on 4 replicas: per replica 1.1–4.4e-5,
  against the forward's 7.1e-5. The infinite-medium runs with the same
  ingredients were unbiased at 16 replicas, with skewed small-sample means.
  Whether this is skew or a residual bias needs the gate-size run.
- **Cost.** 438 s for 4 replicas, half the forward's per history. At the
  gate size the contraction's σ is expected near 2–3 % per detector, and the
  forward's near 0.5 % at arc0.

## LS6 pilot, 1 Hz (jobs 3617424–3617428, commit 6a358d348)

20 replicas × 25 batches × 10⁴ per run. Forward mesh A; adjoints with the
random-ray CADIS windows and cancellation mesh A.

**Forward** (wall time about 47 min per replica, 732 generations per batch):

| detector | δd | σ(δd) | σ(d₀) | δd/d₀ |
|---|---|---|---|---|
| det0 | −4.05e-5 + 4.9e-6i | 0.9 % | 4.4 % | −0.0455 + 0.0055i |
| det1 | −1.14e-5 + 1.4e-6i | 1.1 % | 7.1 % | −0.0417 + 0.0050i |
| det2 | −2.98e-6 + 3.5e-7i | 1.3 % | 16 % | −0.0485 + 0.0057i |
| det3 | −8.85e-7 + 1.1e-7i | 2.1 % | 27 % | −0.0608 + 0.0072i |
| det4 | −3.07e-7 + 3.7e-8i | 3.2 % | 46 % | −0.0667 + 0.0080i |
| det5 | −1.09e-7 + 1.3e-8i | 3.8 % | 58 % | −0.0515 + 0.0061i |
| arc0 | −5.15e-4 + 6.2e-5i | 0.9 % | 1.7 % | −0.0449 + 0.0054i |
| arc1 | −1.58e-4 + 1.9e-5i | 1.0 % | 2.0 % | −0.0444 + 0.0053i |
| arc2 | −4.66e-5 + 5.6e-6i | 0.9 % | 4.7 % | −0.0411 + 0.0049i |
| arc3 | −1.59e-5 + 1.9e-6i | 1.1 % | 8.0 % | −0.0442 + 0.0054i |
| arc4 | −6.24e-6 + 7.5e-7i | 1.6 % | 11 % | −0.0425 + 0.0051i |
| arc5 | −2.80e-6 + 3.4e-7i | 2.1 % | 21 % | −0.0475 + 0.0057i |

- **δF/F₀** = −0.0460 + 0.0053i (σ 0.9 %). N5's 1 Hz value, on a regular
  25 × 25 × 8 mesh, is −0.0467 + 0.0054i, within about 1σ.
- **The noise δd** is sampled to 1–4 % everywhere. The static d₀ is not
  (16–58 % at the far point detectors), which confirms δd for the mesh
  gate. δd/d₀ is −0.041 to −0.049 wherever d₀ is usable: the response is
  global, as point kinetics predicts.
- **Static k** 1.00149 ± 0.00026: +12 ± 26 pcm against `noise.keff`, too
  coarse to test N5's 17 pcm.

**Reciprocity:**

| | forward δd | contraction | \|f/c\| − 1 | phase | Im/Re fwd, contr | paired t (Re, Im) | Hotelling p | FOM adj/fwd |
|---|---|---|---|---|---|---|---|---|
| arc0 | −5.15e-4 (0.9 %) | −4.83e-4 (3.4 %) | +6.6 % | −0.16° | −0.121, −0.118 | −2.09, +1.27 | 0.15 | 0.11 |
| arc3 | −1.59e-5 (1.1 %) | −1.32e-5 (9.0 %) | +20 % | −0.62° | −0.121, −0.110 | −2.54, +0.87 | 0.07 | 0.06 |
| arc5 | | | | | | | | |

- **The phase is now right.** The adjoint's Im/Re matches the forward's,
  and the local pilot's Im deficit was small-sample skew.
- **The adjoint game converges like the forward.** The Σ|w| generation
  ratio is 0.9932 for arc0 against 0.9931 forward, with chains up to 1118
  generations.
- **Both contractions are low,** at 2.1σ and 2.5σ in Re, which combine to
  p ≈ 0.01. Possible causes:
  - small-sample skew, since arc3's per-replica Im has one value 10× the
    rest;
  - smoothing by cancellation in the adjoint (bin averaging moves weight
    within a bin; only the forward has an A-against-B check);
  - a systematic not yet identified.
- **arc5 aborted** at batch 6 on the growth guard. Its first fission
  generation held about one site, because few adjoint particles get from
  81 cm into the core, and the normal rise to 106 sites exceeded the
  guard's 100×. The guard is now 10⁶ for adjoint runs (`--growth-limit`),
  and the gated generation ratio judges boundedness.

## LS6 pilot, second round (1 Hz)

- forward mesh B (needed for the gate anyway);
- adjoint arc0 with cancellation mesh B (tests adjoint cancellation
  smoothing);
- adjoint arc0 and arc3 with a new seed and 4× the batches (skew or bias);
- adjoint arc5 rerun with the raised guard.
