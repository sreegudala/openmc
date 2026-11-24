# Stage 8 log

```
Slab: 4b mixture with water gap [10.0, 14.0] cm, interface x0 = 10.0 vibrating along x (fuel negative side); S_N k = 0.97609437; detector (20.0, 25.0) thermal; cancellation 1 cm
  8a n=1 field dphi/phi0 vs S_N: PASS (STD-C) bins [ok] n=400 bad=3 max|z|=4.70; amp 0 evaluated, max rel err nan; phase 200 evaluated, max 0.155 deg | integral [ok] z Re=-0.06 Im=+0.05 | Hotelling [ok] k=2 T2=0.1 p=0.936 | logged: raw hist mean=-0.00 std=0.38, shape hist mean=-0.22 std=5.50, neighbour corr 1.00; rel sigma per bin median 0.0095 max 0.0105
  8a n=1 detector dd/phi0(D) vs S_N: PASS n=2 bad=0 z_mean=-0.01 z_std=0.10 max|z|=0.09; amp: 0 evaluated, max rel err nan; phase: 1 evaluated, max 0.006 deg; MC -3.85569e-02-8.29390e-01j S_N -3.84885e-02-8.29816e-01j; rel sigma 0.0091; T 37325s; generations/batch 425
  8a n=1 second estimator: PASS paired t Re=+0.04 Im=-0.17, Hotelling p=0.951; forward dd -1.64314e-02-3.53456e-01j contraction -1.64478e-02-3.52883e-01j (P0 +2.33533e-02+4.91501e-01j, P1 -8.81077e-04-1.64781e-02j, fission -3.89201e-02-8.27906e-01j); |fwd/contr|-1 = +0.0016; rel sigma fwd 0.0091 contr 0.0018; adjoint T 54557s, generations/batch 821
    contraction / phi0(D) vs S_N (logged): n=2 bad=0 z_mean=+0.52 z_std=1.07 max|z|=1.27; amp: 1 evaluated, max rel err 2.13e-03; phase: 1 evaluated, max 0.013 deg
  8a n=2 field dphi/phi0 vs S_N: PASS (STD-C) bins [ok] n=400 bad=2 max|z|=4.62; amp 0 evaluated, max rel err nan; phase 0 evaluated, max nan deg | integral [ok] z Re=+1.65 Im=-1.31 | Hotelling [ok] k=2 T2=2.7 p=0.298 | logged: raw hist mean=+0.16 std=1.49, shape hist mean=+0.20 std=1.33, neighbour corr 1.00; rel sigma per bin median 0.0939 max 0.1077
  8a n=2 detector dd/phi0(D) vs S_N: PASS n=2 bad=0 z_mean=+0.20 z_std=2.14 max|z|=1.72; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; MC -1.89539e-02+1.29523e-03j S_N -2.18197e-02+1.59359e-03j; rel sigma 0.0887; T 32728s; generations/batch 259
  8a n=2 second estimator: PASS paired t Re=+1.67 Im=-1.23, Hotelling p=0.292; forward dd -8.07668e-03+5.51809e-04j contraction -9.25760e-03+6.69302e-04j (P0 -2.92339e-03+2.27006e-04j, P1 +1.24060e-04-9.71032e-06j, fission -6.45826e-03+4.52006e-04j); |fwd/contr|-1 = -0.1278; rel sigma fwd 0.0886 contr 0.0026; adjoint T 54677s, generations/batch 822
    contraction / phi0(D) vs S_N (logged): n=2 bad=0 z_mean=-0.13 z_std=2.70 max|z|=2.04; amp: 1 evaluated, max rel err 4.52e-03; phase: 1 evaluated, max 0.042 deg
Test 8b: n = 1 response / eps
  8b eps=0.05 field dphi/phi0 vs S_N: PASS (STD-C) bins [ok] n=400 bad=2 max|z|=5.66; amp 0 evaluated, max rel err nan; phase 200 evaluated, max 0.199 deg | integral [ok] z Re=-1.01 Im=-0.41 | Hotelling [ok] k=2 T2=1.6 p=0.492 | logged: raw hist mean=-0.71 std=0.51, shape hist mean=-0.11 std=3.98, neighbour corr 1.00; rel sigma per bin median 0.0147 max 0.0165
  8b eps=0.05 detector dd/phi0(D) vs S_N: PASS n=2 bad=0 z_mean=-0.63 z_std=0.46 max|z|=0.95; amp: 0 evaluated, max rel err nan; phase: 1 evaluated, max 0.100 deg; MC -1.00266e-02-2.08340e-01j S_N -9.62332e-03-2.07476e-01j; rel sigma 0.0139; T 30809s; generations/batch 241
  8b eps=0.05: detector vs S_N PASS (n=2 bad=0 z_mean=-0.63 z_std=0.46 max|z|=0.95; amp: 0 evaluated, max rel err nan; phase: 1 evaluated, max 0.100 deg); (dd/phi0)/eps vs eps=0.2: PASS Welch t Re=-0.83 Im=-0.29; ratio +1.00486e+00-1.64222e-03j; S_N ratio +1.00010e+00-9.78793e-07j
  8b eps=0.1 field dphi/phi0 vs S_N: PASS (STD-C) bins [ok] n=400 bad=4 max|z|=3.92; amp 0 evaluated, max rel err nan; phase 200 evaluated, max 0.252 deg | integral [ok] z Re=+0.97 Im=+0.89 | Hotelling [ok] k=2 T2=1.0 p=0.642 | logged: raw hist mean=+0.93 std=0.37, shape hist mean=-0.02 std=4.39, neighbour corr 1.00; rel sigma per bin median 0.0131 max 0.0146
  8b eps=0.1 detector dd/phi0(D) vs S_N: PASS n=2 bad=0 z_mean=+0.93 z_std=0.05 max|z|=0.97; amp: 0 evaluated, max rel err nan; phase: 1 evaluated, max 0.030 deg; MC -1.88192e-02-4.10374e-01j S_N -1.92449e-02-4.14915e-01j; rel sigma 0.0123; T 32501s; generations/batch 318
  8b eps=0.1: detector vs S_N PASS (n=2 bad=0 z_mean=+0.93 z_std=0.05 max|z|=0.97; amp: 0 evaluated, max rel err nan; phase: 1 evaluated, max 0.030 deg); (dd/phi0)/eps vs eps=0.2: PASS Welch t Re=+0.78 Im=+0.69; ratio +9.89550e-01+6.21644e-04j; S_N ratio +1.00002e+00-7.54637e-07j

Stage 8 (8a,8b): 10/10 gate checks passed
```
