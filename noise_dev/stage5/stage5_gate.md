# Stage 5 gate record

Gate (design doc, with amendments of 2026-11-18 and 2026-11-19): 5a, 5c
STD; 5b, 5e on-vs-off comparison at h = 1 and 0.5 cm with no trend under an
h/2 origin shift (paired STD-C; h = 2 cm recorded only); 5d median Im FOM
ratio on/off above 1 at the 95 % bootstrap lower bound at 0.1 and 1 Hz, Re
recorded; sum|w| plateaus (every batch terminates, growth guard never
fires).

All runs on TACC LS6, 20 concurrent replicas x 12 threads.

## Runs

| Run | Job | Commit | Settings | Result |
|---|---|---|---|---|
| 5a, 5c | 3592995 | 1d41001dc | 20000 x 20 x 20, seeds 1 + 7919 r | 5/6; 5a 100 Hz Re z = -3.13 (`stage5_log_5a_5c.md`) |
| 5b, 5d, 5e | 3592996 | 1d41001dc | slab 20000 x 20 x 20 | 5b 6/6, 5e 2/2; 5d (original gate) 1 Hz Re 0.79 (`stage5_log_5b_5d_5e.md`) |
| 5a, 5c | 3594120 | c7777124d | 80000 x 20 x 20, seeds 100001 + 7919 r | **6/6** (`stage5_log_5a_5c_p80k.md`) |
| 5d | 3594121 | c7777124d | 20 replicas, seeds 500001 + 7919 r | intervals too wide to decide (`stage5_log_5d_rep.md`) |
| 5d | (rep200) | 484a26d0d | 200 replicas, seeds 700001 + 7919 r | **gate run** (`stage5_log_5d_rep200.md`) |
| 5d tuning | (tune5d) | bd812a01a | 60 x 20 batches, pooled batch variances | `tune5d_log.md` |

## Results

**5a, 5c (gate run, 80000 particles).** Every amplitude within 0.2 % of
the reference, phases within 0.11 deg, max |z| 1.55. With cancellation the
measured generation ratio equals the physical K (0.1 Hz: 0.99474 vs
0.99473; 100 Hz: 0.99232 vs 0.99241) instead of K~cb. 3c at 100 Hz with
eta = 1, divergent without cancellation, now terminates (mean 1022
generations per batch) and meets STD; the 1e4 Hz point of 3a at eta = 1
meets STD with 5e5 particles.

*Disposition of the first-run 5a misses.* The first run had 100 Hz Re
z = -3.13 and amplitudes ~0.86 % low at 0.1 and 1 Hz (|z| 2.3-2.5). The
code was reviewed (bin totals conserved, roulette unbiased with its own
stream, unique ids); in an infinite medium mesh averaging is exactly
unbiased. The replication with 4x particles and fresh seeds shows
deviations of at most 0.19 %; a real 0.86 % bias would have appeared there
at about 5 sigma. This repeats the Stage 3 finding that 20 replicas can
understate the spread of these heavy-tailed weight games.

**5b, 5e.** Paired STD-C passes at every mesh, including h = 2 cm (max |t|
1.31, Hotelling p 0.49-0.74 for 5b; 5e h = 1: p 0.94, h = 0.5: p 0.41).
No trend with the h/2 shift. Max |on/off - 1| over gated bins: 5b 0.3-1.0 %,
5e 0.9-1.8 %, within statistics.

**5d (gate run, 200 replicas), median FOM ratio on/off, 95 % bootstrap:**

| f | Re | Im | Im gate |
|---|---|---|---|
| 0.1 Hz | 1.39 [1.07, 1.80] | 1.84 [1.43, 2.42] | pass |
| 1 Hz | 0.82 [0.63, 1.05] | 2.62 [1.97, 3.51] | pass |

**5d tuning study** (ratios from pooled batch variances; replica/batch
variance cross-check 0.65-1.44):

| Setting | 1 Hz Re | 1 Hz Im | 0.1 Hz Re | 0.1 Hz Im |
|---|---|---|---|---|
| h = 0.5 | 0.90 [0.79, 1.02] | 2.66 [2.31, 3.00] | 1.08 [0.98, 1.20] | 1.76 [1.54, 2.04] |
| h = 1 | 0.95 [0.85, 1.08] | 2.91 [2.57, 3.25] | 1.12 [0.98, 1.28] | 1.71 [1.49, 1.96] |
| h = 2 | 0.88 [0.78, 0.99] | 2.76 [2.38, 3.15] | 0.99 [0.88, 1.10] | 1.60 [1.39, 1.84] |
| h = 4 | 0.92 [0.83, 1.03] | 2.69 [2.36, 3.07] | 1.02 [0.91, 1.16] | 1.69 [1.49, 1.94] |
| h = 1, no bank roulette | 0.92 [0.82, 1.02] | 2.84 [2.50, 3.21] | 1.03 [0.91, 1.18] | 1.61 [1.41, 1.83] |
| h = 1, 4x particles | 0.98 [0.86, 1.12] | 2.92 [2.56, 3.36] | | |

Diagnostics: sum|w| after/before averaging 1.000 (1 Hz), 0.999 (0.1 Hz);
sites per occupied bin 71 (h = 0.5) to 520 (h = 4); bank roulette removes
< 2e-7 of sites. In this problem the in-phase weights all share a sign, so
averaging removes only the phase dispersion carried mainly by delayed
neutrons; the Re variance is population noise that averaging does not
touch. Hence the 2026-11-19 amendment to gate 5d on Im. Cancellation costs
0-8 % in time.

**sum|w|.** Every batch of every run terminated by extinction of the bank;
the growth guard never fired (largest logged generations per batch: 2603).

## Not done

The design doc's plots (sum|w| and |sum w| against generation, t
histograms) were not produced; the data needed for them are in the
`stage5_vals_*.npz` and `tune5d_vals.npz` files (generation sums for 5a/5c,
replica values for all runs).

New options added at this stage: `settings.noise['cancellation_roulette']`
and cancellation diagnostics in the statepoint `noise` group.

Stage 5: PASSED.
