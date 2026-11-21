# Stage 7 log

```
Slab (4b mixture, 100 cm, S_N 400 cells S64): k = 0.99819815; source (-10.0, -5.0) fast, detector (20.0, 25.0) thermal, perturbed cell (30.0, 31.0)
Test 7b: reciprocity, external source (and 7d: S_N fields)
  S_N reciprocity 1 Hz: forward +1.16449e+02-5.78602e+00j, adjoint +1.16449e+02-5.78602e+00j
  7b 1 Hz: PASS two-sample Welch t Re=+1.88 Im=-0.57, Hotelling p=0.189; forward +1.18203e+02-5.95109e+00j adjoint +1.15849e+02-5.81078e+00j; |fwd/adj|-1 = +0.0203, phase diff -0.011 deg; sigma rel fwd 0.0044 adj 0.0101; FOM fwd 6.94 adj 2.26 (adj/fwd 0.33); fwd vs S_N z Re=+3.52 Im=-1.29; adj vs S_N z Re=-0.52 Im=-0.12; S_N +1.16449e+02-5.78602e+00j
    T fwd 7603s adj 4369s; fwd generations/batch 1134; adj generations/batch 857
  7d forward 1 Hz: FAIL (STD-C) bins [FAIL] n=72 bad=36 max|z|=3.65; amp 36 evaluated, max rel err 1.52e-02; phase 36 evaluated, max 0.042 deg | integral [FAIL] z Re=+3.54 Im=-1.24 | Hotelling [ok] k=2 T2=12.7 p=0.010 | logged: raw hist mean=+1.13 std=2.39, shape hist mean=+0.03 std=1.31, neighbour corr 1.00
  7d adjoint 1 Hz: PASS (STD-C) bins [ok] n=80 bad=0 max|z|=0.55; amp 0 evaluated, max rel err nan; phase 40 evaluated, max 0.033 deg | integral [ok] z Re=-0.51 Im=-0.14 | Hotelling [ok] k=2 T2=0.4 p=0.829 | logged: raw hist mean=-0.31 std=0.20, shape hist mean=+0.06 std=0.79, neighbour corr 1.00
  S_N reciprocity 30 Hz: forward +5.42079e+01-6.02855e+01j, adjoint +5.42079e+01-6.02855e+01j
  7b 30 Hz: PASS two-sample Welch t Re=+0.76 Im=-0.62, Hotelling p=0.730; forward +5.39728e+01-6.02129e+01j adjoint +5.34213e+01-5.95276e+01j; |fwd/adj|-1 = +0.0110, phase diff -0.033 deg; sigma rel fwd 0.0063 adj 0.0153; FOM fwd 1.96 adj 0.568 (adj/fwd 0.29); fwd vs S_N z Re=-0.69 Im=+0.19; adj vs S_N z Re=-1.22 Im=+0.73; S_N +5.42079e+01-6.02855e+01j
    T fwd 12805s adj 7496s; fwd generations/batch 1698; adj generations/batch 1272
  7d forward 30 Hz: PASS (STD-C) bins [ok] n=72 bad=0 max|z|=0.85; amp 36 evaluated, max rel err 3.07e-03; phase 36 evaluated, max 0.131 deg | integral [ok] z Re=-0.76 Im=+0.21 | Hotelling [ok] k=3 T2=0.6 p=0.912 | logged: raw hist mean=-0.28 std=0.48, shape hist mean=+0.03 std=0.48, neighbour corr 1.00
  7d adjoint 30 Hz: PASS (STD-C) bins [ok] n=80 bad=0 max|z|=1.43; amp 0 evaluated, max rel err nan; phase 0 evaluated, max nan deg | integral [ok] z Re=-1.26 Im=+0.80 | Hotelling [ok] k=2 T2=1.7 p=0.460 | logged: raw hist mean=-0.23 std=1.06, shape hist mean=-0.13 std=1.10, neighbour corr 1.00
  S_N reciprocity 100 Hz: forward +5.20292e+00-3.18111e+01j, adjoint +5.20292e+00-3.18111e+01j
  7b 100 Hz: PASS two-sample Welch t Re=+0.75 Im=+0.08, Hotelling p=0.739; forward +5.42736e+00-3.16551e+01j adjoint +5.13439e+00-3.16823e+01j; |fwd/adj|-1 = +0.0007, phase diff +0.524 deg; sigma rel fwd 0.0070 adj 0.0145; FOM fwd 1.24 adj 0.463 (adj/fwd 0.37); fwd vs S_N z Re=+1.39 Im=+1.01; adj vs S_N z Re=-0.19 Im=+0.42; S_N +5.20292e+00-3.18111e+01j
    T fwd 16673s adj 10242s; fwd generations/batch 2110; adj generations/batch 1588
  7d forward 100 Hz: PASS (STD-C) bins [ok] n=72 bad=0 max|z|=2.41; amp 36 evaluated, max rel err 8.51e-03; phase 36 evaluated, max 0.619 deg | integral [ok] z Re=+1.34 Im=+1.65 | Hotelling [ok] k=5 T2=16.3 p=0.071 | logged: raw hist mean=+1.45 std=0.44, shape hist mean=+0.02 std=1.07, neighbour corr 1.00
  7d adjoint 100 Hz: PASS (STD-C) bins [ok] n=78 bad=0 max|z|=0.44; amp 0 evaluated, max rel err nan; phase 0 evaluated, max nan deg | integral [ok] z Re=-0.11 Im=+0.24 | Hotelling [ok] k=5 T2=2.6 p=0.839 | logged: raw hist mean=+0.04 std=0.25, shape hist mean=-0.02 std=0.66, neighbour corr 1.00
Test 7c: reciprocity, sampled source in P (contraction pt_lib.pt_terms, P0 + P1)
  7c abs 1 Hz: PASS paired t Re=+0.01 Im=-0.09, Hotelling p=0.996; forward delta d -6.20602e-03+2.74907e-04j contraction -6.20728e-03+2.76250e-04j (P0 -6.21086e-03+2.76398e-04j, P1 +3.57576e-06-1.48770e-07j); |fwd/contr|-1 = -0.0002, phase diff +0.012 deg; sigma rel fwd 0.0084
    vs S_N (delta d / static thermal flux in P, -7.63411e-02+3.35932e-03j): forward z Re=+0.56 Im=+0.05, contraction z Re=+0.53 Im=+0.19; T fwd 4671s; generations/batch 951
  7c abs 30 Hz: PASS paired t Re=-0.15 Im=+0.30, Hotelling p=0.956; forward delta d -3.27043e-03+2.85183e-03j contraction -3.26507e-03+2.83519e-03j (P0 -3.26702e-03+2.83669e-03j, P1 +1.94985e-06-1.49917e-06j); |fwd/contr|-1 = +0.0035, phase diff -0.120 deg; sigma rel fwd 0.0088
    vs S_N (delta d / static thermal flux in P, -4.04048e-02+3.52607e-02j): forward z Re=+1.26 Im=-1.05, contraction z Re=+1.22 Im=-0.96; T fwd 6962s; generations/batch 1409
  7c abs 100 Hz: PASS paired t Re=-0.85 Im=+1.19, Hotelling p=0.435; forward delta d -9.94447e-04+1.64446e-03j contraction -9.76678e-04+1.62546e-03j (P0 -9.77423e-04+1.62631e-03j, P1 +7.44327e-07-8.46369e-07j); |fwd/contr|-1 = +0.0134, phase diff +0.162 deg; sigma rel fwd 0.0081
    vs S_N (delta d / static thermal flux in P, -1.19008e-02+1.99029e-02j): forward z Re=-1.68 Im=+2.12, contraction z Re=-0.23 Im=-0.09; T fwd 8535s; generations/batch 1718
  7c scat 1 Hz: PASS paired t Re=-1.47 Im=+1.15, Hotelling p=0.279; forward delta d +6.72164e-04+2.38129e-04j contraction +1.82591e-03-8.09342e-05j (P0 +1.33619e-03-5.96355e-05j, P1 +4.89724e-04-2.12986e-05j); |fwd/contr|-1 = -0.6098, phase diff +22.046 deg; sigma rel fwd 1.1707
    vs S_N (delta d / static thermal flux in P, +2.24304e-02-9.84509e-04j): forward z Re=-1.47 Im=+1.15, contraction z Re=-0.42 Im=-0.18; T fwd 11623s; generations/batch 1150
  7c scat 30 Hz: PASS paired t Re=-0.24 Im=-2.19, Hotelling p=0.125; forward delta d +8.00494e-04-2.37928e-03j contraction +9.62054e-04-8.31981e-04j (P0 +7.00501e-04-6.16113e-04j, P1 +2.61553e-04-2.15868e-04j); |fwd/contr|-1 = +0.9737, phase diff -30.552 deg; sigma rel fwd 0.3902
    vs S_N (delta d / static thermal flux in P, +1.18998e-02-1.03476e-02j): forward z Re=-0.25 Im=-2.17, contraction z Re=-1.41 Im=+0.99; T fwd 19109s; generations/batch 1840
  7c scat 100 Hz: PASS paired t Re=-1.30 Im=+0.71, Hotelling p=0.367; forward delta d -2.25739e-04-2.46311e-04j contraction +2.91166e-04-4.80577e-04j (P0 +2.03454e-04-3.59295e-04j, P1 +8.77117e-05-1.21281e-04j); |fwd/contr|-1 = -0.4054, phase diff -73.715 deg; sigma rel fwd 1.5478
    vs S_N (delta d / static thermal flux in P, +3.54077e-03-5.88317e-03j): forward z Re=-1.30 Im=+0.71, contraction z Re=+0.26 Im=+0.24; T fwd 24521s; generations/batch 2242
Test 7b low frequency: 40 cm slab, prompt only, S_N k = 0.849394; noise keff = 1 (fixed-source multiplication), 0.1 Hz, eta 1000; thermal adjoint source in (5.0, 10.0)
  S_N: max |Im/Re| of the 0.1 Hz adjoint 1.77e-04; max |noise/static - 1| (Re) 2.39e-08
  lf noise vs static (Re): PASS two-sample STD-C n=30 bad=0 max|t|=0.94; integral t Re=-0.72 Im=+0.00; Hotelling k=5 p=0.421; Re noise/static - 1: max |.| 0.0043
  lf noise vs S_N: PASS (STD-C) bins [ok] n=30 bad=0 max|z|=1.82; amp 15 evaluated, max rel err 1.60e-03; phase 15 evaluated, max 0.000 deg | integral [ok] z Re=+0.47 Im=-1.36 | Hotelling [ok] k=8 T2=8.4 p=0.712 | logged: raw hist mean=-0.37 std=0.88, shape hist mean=+0.00 std=1.10, neighbour corr 0.97
  lf static vs S_N: PASS (STD-C) bins [ok] n=15 bad=0 max|z|=1.70; amp 15 evaluated, max rel err 5.42e-03; phase 15 evaluated, max 0.000 deg | integral [ok] z Re=+1.26 Im=+0.00 | Hotelling [ok] k=5 T2=12.6 p=0.138 | logged: raw hist mean=+1.23 std=0.29, shape hist mean=+0.19 std=1.27, neighbour corr 0.98
    T noise 258s static 125s; noise generations/batch 49
Test 7f: 7b with the 32-bin histogram library (mu_bar 0.3); OpenMC k = 0.985346 +- 0.000089
  7f 1 Hz: PASS two-sample Welch t Re=+0.18 Im=+0.41, Hotelling p=0.891; forward +1.12962e+02-5.53983e+00j adjoint +1.12698e+02-5.66602e+00j; |fwd/adj|-1 = +0.0023, phase diff +0.071 deg; sigma rel fwd 0.0067 adj 0.0115; FOM fwd 3.07 adj 1.82 (adj/fwd 0.59)
    T fwd 7165s adj 4166s; fwd generations/batch 1095; adj generations/batch 833
  7f 30 Hz: PASS two-sample Welch t Re=+0.59 Im=-0.95, Hotelling p=0.610; forward +5.45385e+01-5.78816e+01j adjoint +5.40246e+01-5.72651e+01j; |fwd/adj|-1 = +0.0102, phase diff -0.036 deg; sigma rel fwd 0.0068 adj 0.0119; FOM fwd 1.79 adj 1.01 (adj/fwd 0.57)
    T fwd 11912s adj 6921s; fwd generations/batch 1683; adj generations/batch 1283

Stage 7 (7b,lf,7c,7f): 19/20 gate checks passed
  FAILED: 7d forward 1 Hz
```
