# Stage 4 gate record

Gate (design doc Stage 4, with the STD-C amendment of 2026-11-18 for
correlated mesh tallies): STD-C on mesh cells where |dphi| > 10 % of its
maximum for 4a (0.1, 1, 10, 100 Hz) and 4b (1, 10 Hz); 4c amplitude within
2 % and phase within 1 deg; OpenMC k consistent with the S_N k.

## Gate run

TACC LS6 job 3592829, commit a2e0a4b80, 20 replicas x 20000 particles x 20
batches, seeds 1 + 7919 r. Log: `stage4_log_4a_4b_4c_ksens.md`.
**9/9 gate checks passed.**

| Case | max bin \|z\| | max amp err | max phase err | integral z (Re, Im) | Hotelling k, p |
|---|---|---|---|---|---|
| 4a 0.1 Hz | 1.13 | 0.56 % | 0.15 deg | -0.97, -0.71 | 5, 0.66 |
| 4a 1 Hz | 1.52 | 0.92 % | 0.05 deg | -0.87, +0.00 | 4, 0.19 |
| 4a 10 Hz | 1.53 | 0.58 % | 0.21 deg | +0.41, -1.32 | 5, 0.42 |
| 4a 100 Hz (eta 3) | 2.86 | 1.34 % | 0.27 deg | +0.31, +1.07 | 10, 0.10 |
| 4b 1 Hz | 1.25 | 0.52 % | 0.05 deg | -1.18, +1.13 | 2, 0.35 |
| 4b 10 Hz | 1.14 | 0.23 % | 0.13 deg | +0.19, -1.10 | 2, 0.50 |

- 4c point-kinetic limit: G0_dyn = 152.73 - 8.59i vs G0 = 153.39 - 8.60i,
  amplitude 0.43 %, phase 0.01 deg.
- k: 4a OpenMC 0.999998 +- 0.000019 vs S_N 1.00000406 (z = -0.30); 4b
  OpenMC 0.998087 +- 0.000101 vs S_N 0.99819815 (z = -1.10).
- k sensitivity (4a, 0.1 Hz, k_SN +- 2 sigma_k, 5 replicas each): mean
  |dphi|/ref 0.9960 both ways, phase shift -0.14 / +0.06 deg.
- Pooled raw z over all 580 entries: mean -0.30, std 0.96.

## Follow-ups run alongside (commit 778cb3b53, logged with the old STD
## labels; re-scored here with STD-C)

- 4a 100 Hz, 4x particles (80000), fresh seeds 100001 + 7919 r (job
  3592817, `stage4_log_4a_p80k_100Hz.md`): max |z| 1.46, amplitude 0.35 %,
  phase 0.22 deg, integral z +0.82 / +0.94, Hotelling k = 9, p = 0.79. The
  smooth +2 % (Re) trend toward the right boundary in the gate run falls to
  +0.5 % with half the sigma, so it was a correlated fluctuation, not a
  bias. The S_N reference changes by < 2e-4 on refining to 960 cells or
  S128.
- 4b rerun with the corrected reference normalisation (job 3592818,
  `stage4_log_4b_fixed.md`): identical numbers to the gate run.

## Fixes made during this stage

- 4a at 100 Hz uses eta = 3: eta = 1 gives K~cb = 1.0053 and the growth
  guard fired (TACC job 3589519).
- 4b reference: removed a wrong factor 5 (OpenMC/ref was 0.200 in every
  bin).
- Gate statistic: STD-C (design doc amendment), because mesh bins have
  replica correlation 0.9-1.0 and the raw z-histogram test fails for
  unbiased results.

Stage 4: PASSED.
