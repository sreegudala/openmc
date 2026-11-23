# Stage 7 gate record

Gate (design doc Stage 7, amendments 2026-11-20 and 2026-11-20b):
- **7a:** STD against the `ref_inf` adjoint.
- **7b, 7f:** two-sample test, forward against adjoint.
- **7c:** paired test, forward against the contraction.
- **7d:** STD-C of the forward and adjoint fields against `ref_sn1d`.
- **Low-ω check (7b):** two-sample STD-C of the noise adjoint against a static adjoint run, and each against S_N.
- **7e:** dropped by amendment 2026-11-20b.

The adjoint is the unchanged noise solver run on the transposed library (`common/adjoint_noise_lib.py`, N1 + N2). The bilinear convention applies: same +ω, no conjugation in the run.

All runs on TACC LS6, 20 concurrent replicas × 12 threads.

## Runs

| Run | Job | Commit | Settings | Result |
|---|---|---|---|---|
| 7a | 3600163 | c12f32cd8 | 20 × 20; N3 80000, CORTEX 1e6 particles; seeds 1 + 7919 r | **10/10** (`stage7_log_7a.md`) |
| 7b, 7d, low-ω, 7c, 7f | 3600164 | c12f32cd8 | slab 20000 × 20 × 20; 7c 60000 static particles, nskip 20, 300 inactive; adjoint seeds offset by 1000003 | 19/20 (`stage7_log_7b_lf_7c_7f.md`); miss: 7d forward 1 Hz |
| 7b, 7d at 1 Hz | 3602980 | 208c5db03 | 80000 × 20 × 20, seeds 100001 + 7919 r | **3/3** (`stage7_log_7b_p80k.md`) |
| N3b, N3f (prep) | 3599802 | c12f32cd8 | 80000 × 20 × 20, seeds 100001 + 7919 r | 6/6 (`stage7_prep_log.md`) |

## Results

**7a, infinite medium, full operator, delayed neutrons on.**

| Material | f (Hz) | max amplitude error | σ of mean | max phase error |
|---|---|---|---|---|
| N3 (transposed multiplicity 2.4) | 0.1 / 1 / 10 / 30 / 100 | 0.52 / 0.20 / 0.55 / 0.51 / 0.32 % | 0.5–0.7 % | ≤ 0.67 deg |
| CORTEX fuel | 0.1 / 1 / 10 / 30 / 100 | 0.09 / 0.60 / 0.11 / 0.27 / 1.09 % | 0.9–1.3 % | ≤ 0.12 deg |

- Pooled z over 100 entries: mean +0.38, std 0.99; none at |z| ≥ 3.
- The transposed k equals the forward k to 1e-10 for both materials.

**The adjoint is the more expensive game here.**
- N3 adjoint FOM is 4–11 % of the forward control N3f.
- In the transposed CORTEX fuel only the fast group produces. A thermal adjoint particle reaches it with probability 0.5 %, which is why it needs 1e6 particles per batch.
- Slab FOM adjoint/forward: 0.3–0.8 for 7b, 0.6 for 7f.

**7b, reciprocity with an external source** (4b slab, fast source in [−10, −5] cm, thermal detector in [20, 25] cm; S_N reciprocity exact to 1e-15).

| f | forward / adjoint − 1 | phase diff | Welch t (Re, Im) | run |
|---|---|---|---|---|
| 1 Hz | −0.43 % (σ 0.29 %, 0.42 %) | −0.02 deg | −0.87, −0.09 | p80k (gate) |
| 1 Hz | +2.0 % (σ 0.44 %, 1.0 %) | −0.01 deg | +1.88, −0.57 | first run |
| 30 Hz (break) | +1.1 % (σ 0.63 %, 1.5 %) | −0.03 deg | +0.76, −0.62 | |
| 100 Hz | +0.07 % (σ 0.70 %, 1.45 %) | +0.52 deg | +0.75, +0.08 | |

**7d, S_N fields.** The forward and adjoint 5 cm mesh fields (both groups) pass STD-C at 30 and 100 Hz, and at 1 Hz in the p80k run.
- p80k at 1 Hz: forward max amplitude error 0.30 %, adjoint 0.14 %.

*Disposition of the first-run 7d forward 1 Hz miss.*
- **The miss.** The forward field was 1.5 % above S_N in every bin (integral z +3.5, Hotelling p 0.010, shape normal). The adjoint at 1 Hz agreed with S_N (−0.5 ± 1.0 %).
- **Not the reference.** S_N k is converged to 3e-7, and a 1.5 % shift would need a k error of 8e-5 (0.19 % per 1e-5 at 1 Hz).
- **Not a few outliers.** The replicas lay within 0.97–1.05 of their mean, with no outliers.
- **The rerun.** 4× particles and fresh seeds give −0.29 % (z −1.0, σ 0.29 %); a real +1.5 % would have appeared at about 5σ.
- **Conclusion.** This is the third case (after 3b and 5a) of 20 small replicas understating the spread of these weight games.

**Low-ω check** (40 cm prompt-only slab, k = 0.849, noise keff = 1, 0.1 Hz, η 1000). The noise adjoint against a plain OpenMC fixed-source adjoint on the same transposed library:
- **Noise against static, Re:** max |Re ratio − 1| = 0.43 %, Hotelling p 0.42.
- **Each against S_N:** both pass. The noise run's Im is checked against S_N, where |Im/Re| = 1.8e-4.

**7c, reciprocity with the sampled source** (1 cm cell P = [30, 31] cm; contraction `pt_lib.pt_terms`, copied from ce_mg_triga 540ab07, P0 + P1, applied to complex χ by linearity).

| Case | f | forward / contraction − 1 | forward σ | contraction / S_N − 1 (σ) | P1 share |
|---|---|---|---|---|---|
| absorber | 1 / 30 / 100 Hz | −0.02 / +0.35 / +1.3 % | 0.8–0.9 % | −0.5 / −1.3 / +0.03 % (0.9–1.3 %) | 0.1 % |
| scattering | 1 / 30 / 100 Hz | no power | 40–155 % | −0.3 / −1.4 / −0.07 % (0.8–1.3 %) | 27 % |

- **Absorber.** The paired test passes with σ under 1 %. Forward and contraction also agree with S_N.
- **Scattering.** The paired test passes, but **without power**.
  - The forward sampled source is a copy (−δΣ_t ψ₀) and a sampled scatter (+δΣ_s), which cancel on average because δΣ_a = 0. The forward σ is therefore 40–155 %.
  - The contraction itself agrees with S_N within 1σ at every frequency. P1 is 27 % of the answer, so a wrong P1 sign or factor would be 40–65σ off.
  - This checks the contraction (sign, factor 3, OpenMC's Y1m normalisation, complex χ) against S_N, not against forward Monte Carlo. Recorded as such.

**7f, histogram library** (32 bins, μ̄ 0.3; OpenMC k 0.985346 ± 0.000089).
- 1 Hz: forward/adjoint − 1 = +0.23 %.
- 30 Hz: +1.0 %.
- Welch |t| ≤ 0.95, Hotelling p ≥ 0.61.

## Verdict

Stage 7 **passed**: 7a 10/10, 7b 3/3, 7d 6/6, low-ω 3/3, 7c 6/6 (the scattering case by contraction against S_N, see above), 7f 2/2.

The adjoint noise solution from the transposed library reproduces:
- the exact adjoint in infinite media, with multiplicity 2.4;
- the forward detector response by reciprocity in a two-group slab, with isotropic and histogram scattering;
- the S_N adjoint field;
- the static adjoint at low ω.

**Carried forward:**
- **Cancelling sampled sources.** These are very noisy (7c scattering). Stage 8's vibration source has the same structure.
- **Adjoint efficiency.** The adjoint can be 10–25 times less efficient than the forward game (N3, CORTEX). The TRIGA application will need weight windows (N3c).
