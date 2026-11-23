# Stage 7 log

```
Slab (4b mixture, 100 cm, S_N 400 cells S64): k = 0.99819815; source (-10.0, -5.0) fast, detector (20.0, 25.0) thermal, perturbed cell (30.0, 31.0)
Test 7b: reciprocity, external source (and 7d: S_N fields)
  S_N reciprocity 1 Hz: forward +1.16449e+02-5.78602e+00j, adjoint +1.16449e+02-5.78602e+00j
  7b 1 Hz: PASS two-sample Welch t Re=-0.87 Im=-0.09, Hotelling p=0.597; forward +1.16108e+02-5.81088e+00j adjoint +1.16610e+02-5.79976e+00j; |fwd/adj|-1 = -0.0043, phase diff -0.018 deg; sigma rel fwd 0.0029 adj 0.0042; FOM fwd 4.1 adj 3.3 (adj/fwd 0.81); fwd vs S_N z Re=-1.04 Im=-0.42; adj vs S_N z Re=+0.34 Im=-0.13; S_N +1.16449e+02-5.78602e+00j
    T fwd 29812s adj 17563s; fwd generations/batch 1379; adj generations/batch 1118
  7d forward 1 Hz: PASS (STD-C) bins [ok] n=72 bad=0 max|z|=1.06; amp 36 evaluated, max rel err 3.01e-03; phase 36 evaluated, max 0.024 deg | integral [ok] z Re=-1.01 Im=-0.47 | Hotelling [ok] k=2 T2=2.3 p=0.351 | logged: raw hist mean=-0.74 std=0.29, shape hist mean=+0.03 std=0.89, neighbour corr 1.00
  7d adjoint 1 Hz: PASS (STD-C) bins [ok] n=80 bad=0 max|z|=0.37; amp 40 evaluated, max rel err 1.43e-03; phase 40 evaluated, max 0.015 deg | integral [ok] z Re=+0.31 Im=-0.11 | Hotelling [ok] k=2 T2=0.1 p=0.958 | logged: raw hist mean=+0.09 std=0.23, shape hist mean=-0.02 std=0.87, neighbour corr 1.00

Stage 7 (7b): 3/3 gate checks passed
```
