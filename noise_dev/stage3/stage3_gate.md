# Stage 3 gate record

Gate: STD for 3a, 3b, 3d; generation ratio within 1 % of K~cb in convergent
cases; 3c guard fires with ratio > 1 and within 1 % of K~cb.

## Runs

| Run | Where | Particles x batches x replicas | Seeds | Result |
|---|---|---|---|---|
| 3a, 3b, 3c, 3d | local, 12 threads | 5000 x 20 x 20 | 1 + 7919 r | 37/38; 3b 1 Hz Re z = 3.05 (`stage3_log_3a_3b_3c_3d.md`) |
| 3b (all 7 frequencies) | TACC LS6, 20 x 12 threads | 20000 x 20 x 20 | 100001 + 7919 r | 14/14 (`stage3_log_3b_p20k.md`) |

## Disposition of the 3b 1 Hz miss

Local: 1377.0 - 157.2i vs reference 1404.2 - 168.1i (-1.9 % in Re,
z = 3.05). Replication with fresh seeds and 4x particles: 1408.4 - 167.5i
(+0.29 %, max |z| = 0.51, amplitude 0.29 %, phase 0.05 deg). A real -1.9 %
bias would have appeared there at about z = -3; it is excluded. The full 3b
replication passes every frequency with pooled z mean -0.23, std 0.99
(n = 28); the local 5000-particle run had pooled std 1.45 (n = 116, two
nearly identical estimators per case). Small replicas evidently understate
the spread of these heavy-tailed weight games: later gates use >= 20000
particles per replica.

The generation-ratio check (measured sum|w| ratio vs K~cb) agrees to
< 1.1e-4 relative in every convergent case of both runs.

## Other observations

- 3c: guard fires at a growth factor of 100; measured ratio 1.00510 vs
  K~cb 1.00526 (rel. diff 1.6e-4).
- The TACC replication's 0.1 and 0.3 Hz results are bit-identical between a
  4-thread run and a 12-thread run with the same seeds (jobs 3589518 and
  3591800): results do not depend on the thread count.

Stage 3: PASSED.
