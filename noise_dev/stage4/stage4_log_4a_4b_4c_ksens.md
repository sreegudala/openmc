# Stage 4 log

```
Stage 4a: S_N (480 cells, S64) k = 1.00000406
  OpenMC eigenvalue k = 0.999998 +- 0.000019; (k_MC - k_SN)/sigma = -0.30
  4a 0.1 Hz: eta = 1.0, infinite-medium K~cb = 0.99432
  4a 0.1 Hz: PASS (STD-C) bins [ok] n=112 bad=0 max|z|=1.13; amp 56 evaluated, max rel err 5.58e-03; phase 56 evaluated, max 0.152 deg | integral [ok] z Re=-0.97 Im=-0.71 | Hotelling [ok] k=5 T2=4.2 p=0.658 | logged: raw hist mean=-0.79 std=0.19, shape hist mean=-0.09 std=0.84, neighbour corr 1.00
    diffusion vs S_N (informational): max |rel diff| 0.009; T=3170s
  4a 1 Hz: eta = 1.0, infinite-medium K~cb = 0.99370
  4a 1 Hz: PASS (STD-C) bins [ok] n=112 bad=0 max|z|=1.52; amp 56 evaluated, max rel err 9.23e-03; phase 56 evaluated, max 0.054 deg | integral [ok] z Re=-0.87 Im=+0.00 | Hotelling [ok] k=4 T2=8.3 p=0.189 | logged: raw hist mean=-0.43 std=0.56, shape hist mean=-0.09 std=1.25, neighbour corr 1.00
    diffusion vs S_N (informational): max |rel diff| 0.009; T=3090s
  4a 10 Hz: eta = 1.0, infinite-medium K~cb = 0.99467
  4a 10 Hz: PASS (STD-C) bins [ok] n=112 bad=0 max|z|=1.53; amp 56 evaluated, max rel err 5.78e-03; phase 56 evaluated, max 0.213 deg | integral [ok] z Re=+0.41 Im=-1.32 | Hotelling [ok] k=5 T2=6.8 p=0.416 | logged: raw hist mean=-0.43 std=0.83, shape hist mean=+0.04 std=1.10, neighbour corr 1.00
    diffusion vs S_N (informational): max |rel diff| 0.010; T=3690s
  4a 100 Hz: eta = 3.0, infinite-medium K~cb = 0.99809
  4a 100 Hz: PASS (STD-C) bins [ok] n=100 bad=0 max|z|=2.86; amp 32 evaluated, max rel err 1.34e-02; phase 24 evaluated, max 0.268 deg | integral [ok] z Re=+0.31 Im=+1.07 | Hotelling [ok] k=10 T2=44.1 p=0.100 | logged: raw hist mean=+0.44 std=1.43, shape hist mean=+0.10 std=1.66, neighbour corr 0.99
    diffusion vs S_N (informational): max |rel diff| 0.022; T=10364s
  k sensitivity 0.1 Hz, k = k_SN -2 sigma_k (0.999966): mean |dphi| ratio to reference 0.9960, mean phase shift -0.14 deg
  k sensitivity 0.1 Hz, k = k_SN +2 sigma_k (1.000042): mean |dphi| ratio to reference 0.9960, mean phase shift +0.06 deg
Test 4c: G0_dyn = 152.730-8.592j (|.|=152.97, -3.22 deg) vs G0 = 153.388-8.596j (|.|=153.63, -3.21 deg): amp diff 0.428%, phase diff -0.01 deg PASS
Stage 4b: two-group slab
  homogenized fuel/water k_inf = 1.03790
  critical width ~ 102.55 cm; using 100 cm
  S_N k = 0.99819815
  OpenMC k = 0.998087 +- 0.000101; z = -1.10
  4b 1 Hz: PASS (STD-C) bins [ok] n=72 bad=0 max|z|=1.25; amp 36 evaluated, max rel err 5.19e-03; phase 36 evaluated, max 0.051 deg | integral [ok] z Re=-1.18 Im=+1.13 | Hotelling [ok] k=2 T2=2.3 p=0.353 | logged: raw hist mean=-0.02 std=1.16, shape hist mean=-0.04 std=1.13, neighbour corr 1.00
  4b 10 Hz: PASS (STD-C) bins [ok] n=72 bad=0 max|z|=1.14; amp 36 evaluated, max rel err 2.25e-03; phase 36 evaluated, max 0.133 deg | integral [ok] z Re=+0.19 Im=-1.10 | Hotelling [ok] k=2 T2=1.5 p=0.502 | logged: raw hist mean=-0.45 std=0.65, shape hist mean=-0.04 std=0.78, neighbour corr 1.00
Pooled z: n=580 mean=-0.30 std=0.96

Stage 4 (4a,4b,4c,ksens): 9/9 gate checks passed
```
