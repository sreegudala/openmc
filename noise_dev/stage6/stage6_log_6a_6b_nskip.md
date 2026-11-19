# Stage 6 log

```
Test 6a: infinite medium (case 3b), uniform dSigma_a/Sigma_a = 1e-3
  6a 0.1 Hz (eta 1.0): PASS n=2 bad=0 z_mean=+0.03 z_std=0.57 max|z|=0.43; amp: 1 evaluated, max rel err 1.62e-03; phase: 1 evaluated, max 0.021 deg; dphi/phi0 = -0.167053+0.062054j ref -0.166806+0.061891j (-G0 eps -0.166806+0.061891j)
  6a 0.1 Hz source integral: mean -9.989490e-04+0.000000e+00j vs -dSigma_t phi_0 -9.988541e-04; z Re -0.33, Im +0.00 -> PASS
    static k 0.99994 +- 0.00022 (per generation, correlated); source sites/batch 36267; generations/batch 1507; T=11670s
  6a 0.3 Hz (eta 1.0): PASS n=2 bad=0 z_mean=+0.44 z_std=0.28 max|z|=0.64; amp: 1 evaluated, max rel err 1.85e-03; phase: 1 evaluated, max 0.046 deg; dphi/phi0 = -0.148719+0.030862j ref -0.149019+0.030800j (-G0 eps -0.149019+0.030800j)
  6a 0.3 Hz source integral: mean -9.999891e-04+0.000000e+00j vs -dSigma_t phi_0 -1.000065e-03; z Re +0.26, Im +0.00 -> PASS
    static k 0.99995 +- 0.00023 (per generation, correlated); source sites/batch 36305; generations/batch 1175; T=8731s
  6a 1 Hz (eta 1.0): PASS n=2 bad=0 z_mean=-0.75 z_std=0.91 max|z|=1.40; amp: 1 evaluated, max rel err 1.24e-04; phase: 1 evaluated, max 0.068 deg; dphi/phi0 = -0.140462+0.016649j ref -0.140425+0.016815j (-G0 eps -0.140425+0.016815j)
  6a 1 Hz source integral: mean -1.000162e-03+0.000000e+00j vs -dSigma_t phi_0 -9.997549e-04; z Re -1.25, Im +0.00 -> PASS
    static k 1.00017 +- 0.00023 (per generation, correlated); source sites/batch 36311; generations/batch 1039; T=7689s
  6a 3 Hz (eta 1.0): PASS n=2 bad=0 z_mean=-0.04 z_std=0.59 max|z|=0.46; amp: 1 evaluated, max rel err 1.22e-03; phase: 1 evaluated, max 0.024 deg; dphi/phi0 = -0.136823+0.020815j ref -0.136982+0.020899j (-G0 eps -0.136982+0.020899j)
  6a 3 Hz source integral: mean -9.999949e-04+0.000000e+00j vs -dSigma_t phi_0 -9.999752e-04; z Re -0.10, Im +0.00 -> PASS
    static k 0.99967 +- 0.00023 (per generation, correlated); source sites/batch 36305; generations/batch 1022; T=7861s
  6a 10 Hz (eta 1.0): PASS n=2 bad=0 z_mean=+0.95 z_std=1.03 max|z|=1.68; amp: 1 evaluated, max rel err 3.11e-03; phase: 1 evaluated, max 0.107 deg; dphi/phi0 = -0.118028+0.050404j ref -0.118490+0.050340j (-G0 eps -0.118490+0.050340j)
  6a 10 Hz source integral: mean -1.000176e-03+0.000000e+00j vs -dSigma_t phi_0 -1.000149e-03; z Re -0.11, Im +0.00 -> PASS
    static k 0.99986 +- 0.00023 (per generation, correlated); source sites/batch 36312; generations/batch 1138; T=8455s
  6a 30 Hz (eta 1.0): PASS n=2 bad=0 z_mean=+0.42 z_std=0.65 max|z|=0.88; amp: 1 evaluated, max rel err 1.75e-03; phase: 1 evaluated, max 0.113 deg; dphi/phi0 = -0.054225+0.068181j ref -0.054455+0.068193j (-G0 eps -0.054455+0.068193j)
  6a 30 Hz source integral: mean -1.000128e-03+0.000000e+00j vs -dSigma_t phi_0 -1.000214e-03; z Re +0.24, Im +0.00 -> PASS
    static k 1.00006 +- 0.00023 (per generation, correlated); source sites/batch 36310; generations/batch 1776; T=14493s
  6a 100 Hz (eta 3.0): PASS n=2 bad=0 z_mean=+1.70 z_std=0.37 max|z|=1.96; amp: 1 evaluated, max rel err 6.61e-03; phase: 1 evaluated, max 0.374 deg; dphi/phi0 = -0.007461+0.032000j ref -0.007619+0.031741j (-G0 eps -0.007619+0.031741j)
  6a 100 Hz source integral: mean -1.000509e-03+0.000000e+00j vs -dSigma_t phi_0 -1.000036e-03; z Re -1.54, Im +0.00 -> PASS
    static k 1.00000 +- 0.00023 (per generation, correlated); source sites/batch 36324; generations/batch 2413; T=21357s
Test 6b: slab 4a, dSigma_a/Sigma_a = 0.01 in [14, 15] (dSigma_t/Sigma_t = 0.001); S_N k = 1.00000406
  6b 0.1 Hz (eta 1.0, nskip 20): PASS (STD-C) bins [ok] n=120 bad=0 max|z|=1.41; amp 33 evaluated, max rel err 5.86e-03; phase 60 evaluated, max 0.047 deg | integral [ok] z Re=+0.27 Im=+0.17 | Hotelling [ok] k=9 T2=10.3 p=0.726 | logged: raw hist mean=+0.21 std=0.40, shape hist mean=-0.17 std=0.91, neighbour corr 0.99
    dphi/phi0 at [14,15]: -4.30595e-02+3.77494e-03j ref -4.31645e-02+3.76317e-03j; at [0,1]: -2.46930e-02+3.46576e-03j ref -2.48264e-02+3.46751e-03j
  6b 0.1 Hz source integral: mean -1.870395e-04+0.000000e+00j vs -dSigma_t phi_0 -1.869041e-04; z Re -2.07, Im +0.00 -> PASS
    static k 0.99997 (S_N 1.00000); source sites/batch 20336; generations/batch 1034; T=4221s
  6b 1 Hz (eta 1.0, nskip 20): PASS (STD-C) bins [ok] n=120 bad=0 max|z|=1.73; amp 58 evaluated, max rel err 9.94e-03; phase 60 evaluated, max 0.029 deg | integral [ok] z Re=+1.33 Im=-0.29 | Hotelling [ok] k=11 T2=7.0 p=0.968 | logged: raw hist mean=+0.45 std=0.81, shape hist mean=-0.19 std=0.77, neighbour corr 0.98
    dphi/phi0 at [14,15]: -4.29552e-02+1.61623e-03j ref -4.30861e-02+1.62939e-03j; at [0,1]: -2.46360e-02+1.49493e-03j ref -2.47778e-02+1.50030e-03j
  6b 1 Hz source integral: mean -1.866295e-04+0.000000e+00j vs -dSigma_t phi_0 -1.866473e-04; z Re +0.20, Im +0.00 -> PASS
    static k 1.00004 (S_N 1.00000); source sites/batch 20298; generations/batch 958; T=3947s
  6b 10 Hz (eta 1.0, nskip 20): PASS (STD-C) bins [ok] n=120 bad=0 max|z|=1.38; amp 60 evaluated, max rel err 5.59e-03; phase 60 evaluated, max 0.107 deg | integral [ok] z Re=+0.49 Im=-0.61 | Hotelling [ok] k=11 T2=10.6 p=0.889 | logged: raw hist mean=-0.06 std=0.65, shape hist mean=-0.06 std=0.70, neighbour corr 0.99
    dphi/phi0 at [14,15]: -3.87125e-02+1.08042e-02j ref -3.87719e-02+1.08414e-02j; at [0,1]: -2.04802e-02+9.79265e-03j ref -2.05404e-02+9.84253e-03j
  6b 10 Hz source integral: mean -1.873737e-04+0.000000e+00j vs -dSigma_t phi_0 -1.874434e-04; z Re +0.63, Im +0.00 -> PASS
    static k 0.99998 (S_N 1.00000); source sites/batch 20377; generations/batch 1087; T=4388s
  6b 100 Hz (eta 3.0, nskip 20): PASS (STD-C) bins [ok] n=120 bad=0 max|z|=2.49; amp 25 evaluated, max rel err 6.96e-03; phase 21 evaluated, max 0.472 deg | integral [ok] z Re=-0.49 Im=-1.98 | Hotelling [ok] k=12 T2=20.3 p=0.714 | logged: raw hist mean=-0.86 std=1.00, shape hist mean=-0.27 std=1.75, neighbour corr 0.99
    dphi/phi0 at [14,15]: -1.31113e-02+1.15969e-02j ref -1.30359e-02+1.16132e-02j; at [0,1]: 7.22477e-04+5.26801e-03j ref 6.82565e-04+5.34440e-03j
  6b 100 Hz source integral: mean -1.866838e-04+0.000000e+00j vs -dSigma_t phi_0 -1.866168e-04; z Re -0.81, Im +0.00 -> PASS
    static k 0.99995 (S_N 1.00000); source sites/batch 20303; generations/batch 2457; T=10547s
nskip sensitivity (logged, not gated): 6b at 1 Hz, nskip 3
  6b_nskip3 1 Hz (eta 1.0, nskip 3): PASS (STD-C) bins [ok] n=120 bad=0 max|z|=1.64; amp 33 evaluated, max rel err 4.35e-03; phase 60 evaluated, max 0.118 deg | integral [ok] z Re=+0.31 Im=-0.71 | Hotelling [ok] k=9 T2=4.7 p=0.958 | logged: raw hist mean=-0.20 std=0.66, shape hist mean=+0.09 std=1.35, neighbour corr 0.99
    dphi/phi0 at [14,15]: -4.31314e-02+1.63743e-03j ref -4.30861e-02+1.62939e-03j; at [0,1]: -2.47809e-02+1.48374e-03j ref -2.47778e-02+1.50030e-03j
  6b_nskip3 1 Hz source integral: mean -1.861827e-04+0.000000e+00j vs -dSigma_t phi_0 -1.861062e-04; z Re -1.33, Im +0.00 -> PASS
    static k 0.99995 (S_N 1.00000); source sites/batch 20249; generations/batch 953; T=3493s

Stage 6 (6a,6b,nskip): 22/22 gate checks passed
```
