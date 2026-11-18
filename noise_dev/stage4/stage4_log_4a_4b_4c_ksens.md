# Stage 4 log

```
Stage 4a: S_N (480 cells, S64) k = 1.00000406
  OpenMC eigenvalue k = 0.999998 +- 0.000019; (k_MC - k_SN)/sigma = -0.30
  4a 0.1 Hz: eta = 1.0, infinite-medium K~cb = 0.99432
  4a 0.1 Hz: FAIL n=112 bad=0 z_mean=-0.79 z_std=0.19 max|z|=1.13; amp: 56 evaluated, max rel err 5.58e-03; phase: 56 evaluated, max 0.152 deg
    integral z: Re -0.97 Im -0.71; shape z: n=112 mean=-0.09 std=0.84 max|z|=1.68; replica corr (Re) neighbour 1.00, >=10 bins apart 0.95
    diffusion vs S_N (informational): max |rel diff| 0.009; T=3570s
  4a 1 Hz: eta = 1.0, infinite-medium K~cb = 0.99370
  4a 1 Hz: FAIL n=112 bad=0 z_mean=-0.43 z_std=0.56 max|z|=1.52; amp: 56 evaluated, max rel err 9.23e-03; phase: 56 evaluated, max 0.054 deg
    integral z: Re -0.87 Im +0.00; shape z: n=112 mean=-0.09 std=1.25 max|z|=2.18; replica corr (Re) neighbour 1.00, >=10 bins apart 0.92
    diffusion vs S_N (informational): max |rel diff| 0.009; T=3263s
  4a 10 Hz: eta = 1.0, infinite-medium K~cb = 0.99467
  4a 10 Hz: FAIL n=112 bad=0 z_mean=-0.43 z_std=0.83 max|z|=1.53; amp: 56 evaluated, max rel err 5.78e-03; phase: 56 evaluated, max 0.213 deg
    integral z: Re +0.41 Im -1.32; shape z: n=112 mean=+0.04 std=1.10 max|z|=2.16; replica corr (Re) neighbour 1.00, >=10 bins apart 0.90
    diffusion vs S_N (informational): max |rel diff| 0.010; T=3846s
  4a 100 Hz: eta = 3.0, infinite-medium K~cb = 0.99809
  4a 100 Hz: FAIL n=100 bad=0 z_mean=+0.44 z_std=1.43 max|z|=2.86; amp: 32 evaluated, max rel err 1.34e-02; phase: 24 evaluated, max 0.268 deg
    integral z: Re +0.31 Im +1.07; shape z: n=100 mean=+0.10 std=1.66 max|z|=3.91; replica corr (Re) neighbour 0.99, >=10 bins apart 0.53
    diffusion vs S_N (informational): max |rel diff| 0.022; T=10872s
  k sensitivity 0.1 Hz, k = k_SN -2 sigma_k (0.999966): mean |dphi| ratio to reference 0.9960, mean phase shift -0.14 deg
  k sensitivity 0.1 Hz, k = k_SN +2 sigma_k (1.000042): mean |dphi| ratio to reference 0.9960, mean phase shift +0.06 deg
Test 4c: G0_dyn = 152.730-8.592j (|.|=152.97, -3.22 deg) vs G0 = 153.388-8.596j (|.|=153.63, -3.21 deg): amp diff 0.428%, phase diff -0.01 deg PASS
Stage 4b: two-group slab
  homogenized fuel/water k_inf = 1.03790
  critical width ~ 102.55 cm; using 100 cm
  S_N k = 0.99819815
  OpenMC k = 0.998087 +- 0.000101; z = -1.10
  4b 1 Hz: FAIL n=72 bad=72 z_mean=-368.91 z_std=592.75 max|z|=1037.74; amp: 36 evaluated, max rel err 8.01e-01; phase: 36 evaluated, max 0.051 deg
    integral z: Re -961.92 Im +219.26; shape z: n=72 mean=-0.04 std=1.13 max|z|=2.06; replica corr (Re) neighbour 1.00, >=10 bins apart 0.99
  4b 10 Hz: FAIL n=72 bad=72 z_mean=-151.95 z_std=711.74 max|z|=932.90; amp: 36 evaluated, max rel err 8.00e-01; phase: 36 evaluated, max 0.133 deg
    integral z: Re -863.41 Im +555.16; shape z: n=72 mean=-0.04 std=0.78 max|z|=2.02; replica corr (Re) neighbour 1.00, >=10 bins apart 1.00
Pooled z: n=580 mean=-64.90 std=347.53

Stage 4 (4a,4b,4c,ksens): 3/9 gate checks passed
  FAILED: 4a 0.1 Hz
  FAILED: 4a 1 Hz
  FAILED: 4a 10 Hz
  FAILED: 4a 100 Hz
  FAILED: 4b 1 Hz
  FAILED: 4b 10 Hz
```
