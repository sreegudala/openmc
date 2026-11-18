# Stage 4 log

```
Stage 4a: S_N (480 cells, S64) k = 1.00000406
Stage 4b: two-group slab
  homogenized fuel/water k_inf = 1.03790
  critical width ~ 102.55 cm; using 100 cm
  S_N k = 0.99819815
  OpenMC k = 0.998087 +- 0.000101; z = -1.10
  4b 1 Hz: PASS n=72 bad=0 z_mean=-0.02 z_std=1.16 max|z|=1.25; amp: 36 evaluated, max rel err 5.19e-03; phase: 36 evaluated, max 0.051 deg
    integral z: Re -1.18 Im +1.13; shape z: n=72 mean=-0.04 std=1.13 max|z|=2.06; replica corr (Re) neighbour 1.00, >=10 bins apart 0.99
  4b 10 Hz: FAIL n=72 bad=0 z_mean=-0.45 z_std=0.65 max|z|=1.14; amp: 36 evaluated, max rel err 2.25e-03; phase: 36 evaluated, max 0.133 deg
    integral z: Re +0.19 Im -1.10; shape z: n=72 mean=-0.04 std=0.78 max|z|=2.02; replica corr (Re) neighbour 1.00, >=10 bins apart 1.00
Pooled z: n=144 mean=-0.24 std=0.96

Stage 4 (4b): 2/3 gate checks passed
  FAILED: 4b 10 Hz
```
