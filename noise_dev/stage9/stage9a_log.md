# Stage 9a log

```
Stage 9a: C5G7 noise, OpenMC against FEMFFUSION ref1_FE2_SP3
  k: OpenMC 1.186489 +- 0.000022 (C5G7 transport reference 1.18655); FEMFFUSION ref1_FE2_SP3 1.18271, ref0_FE1_SP3 1.18539, ref1_FE2_SP1 1.18335, ref1_FE2_SP5 1.18220, ref1_FE2_SP3_nondiag 1.18271, ref1_FE2_SP5_nondiag 1.18220
  k sigma 2.21 pcm -> PASS (3 pcm)
  FEMFFUSION ref0_FE1_SP3 vs ref1_FE2_SP3 (outside the perturbed assembly): fuel max |amp ratio - 1| 0.0181, max |dphase| 0.06 deg; reflector max |amp ratio - 1| 0.0096, max |dphase| 0.02 deg
  FEMFFUSION ref1_FE2_SP1 vs ref1_FE2_SP3 (outside the perturbed assembly): fuel max |amp ratio - 1| 0.0246, max |dphase| 0.06 deg; reflector max |amp ratio - 1| 0.1272, max |dphase| 0.03 deg
  FEMFFUSION ref1_FE2_SP5 vs ref1_FE2_SP3 (outside the perturbed assembly): fuel max |amp ratio - 1| 0.0053, max |dphase| 0.03 deg; reflector max |amp ratio - 1| 0.0054, max |dphase| 0.02 deg
  FEMFFUSION ref1_FE2_SP3_nondiag vs ref1_FE2_SP3 (outside the perturbed assembly): fuel max |amp ratio - 1| 0.0015, max |dphase| 0.00 deg; reflector max |amp ratio - 1| 0.5431, max |dphase| 0.01 deg
  FEMFFUSION ref1_FE2_SP5_nondiag vs ref1_FE2_SP3 (outside the perturbed assembly): fuel max |amp ratio - 1| 0.0042, max |dphase| 0.02 deg; reflector max |amp ratio - 1| 0.7136, max |dphase| 0.02 deg

  OpenMC run 'c51': 20 replicas, cancellation 51, keff 1.186489, T=51086s
    g1: PASS: evaluated amp 2045/2312, phase 2310/2312; max |amp ratio - 1| 0.0431, max |dphase| 0.36 deg; all cells: amp ratio - 1 median -0.0079 [5-95 % -0.0213, +0.0089], dphase median -0.09 deg; z hist (logged) mean -0.25 std 0.93
    g2: PASS: evaluated amp 1993/2312, phase 2303/2312; max |amp ratio - 1| 0.0401, max |dphase| 0.70 deg; all cells: amp ratio - 1 median -0.0080 [5-95 % -0.0206, +0.0114], dphase median -0.09 deg; z hist (logged) mean -0.24 std 0.92
    g3: FAIL: evaluated amp 1866/2312, phase 2299/2312; max |amp ratio - 1| 0.0566, max |dphase| 0.61 deg; all cells: amp ratio - 1 median -0.0076 [5-95 % -0.0209, +0.0218], dphase median -0.08 deg; z hist (logged) mean -0.15 std 0.94
    g4: FAIL: evaluated amp 1752/2312, phase 2269/2312; max |amp ratio - 1| 0.0537, max |dphase| 0.61 deg; all cells: amp ratio - 1 median -0.0074 [5-95 % -0.0221, +0.0339], dphase median -0.09 deg; z hist (logged) mean -0.12 std 0.94
    g5: PASS: evaluated amp 1710/2312, phase 2264/2312; max |amp ratio - 1| 0.0432, max |dphase| 0.94 deg; all cells: amp ratio - 1 median -0.0074 [5-95 % -0.0216, +0.0506], dphase median -0.09 deg; z hist (logged) mean -0.07 std 0.99
    g6: PASS: evaluated amp 2100/2312, phase 2303/2312; max |amp ratio - 1| 0.0398, max |dphase| 0.61 deg; all cells: amp ratio - 1 median -0.0078 [5-95 % -0.0173, +0.0100], dphase median -0.09 deg; z hist (logged) mean -0.22 std 0.93
    g7: PASS: evaluated amp 2166/2312, phase 2307/2312; max |amp ratio - 1| 0.0447, max |dphase| 0.47 deg; all cells: amp ratio - 1 median -0.0078 [5-95 % -0.0167, +0.0056], dphase median -0.09 deg; z hist (logged) mean -0.24 std 0.92
    against ref1_FE2_SP3_nondiag (logged): g1 FAIL (max amp 0.370), g2 FAIL (max amp 0.225), g3 FAIL (max amp 0.140), g4 FAIL (max amp 0.120), g5 FAIL (max amp 0.095), g6 FAIL (max amp 0.301), g7 FAIL (max amp 0.281)
    g1 perturbed pin: OpenMC 0.01245 at  +86.24 deg, FEMFFUSION 0.01152 at  +85.97 deg
    g4 perturbed pin: OpenMC 0.01385 at  +86.63 deg, FEMFFUSION 0.01334 at  +86.53 deg
    g7 perturbed pin: OpenMC 0.02879 at  +88.38 deg, FEMFFUSION 0.02878 at  +88.39 deg
    g1 perturbed assembly (logged): amp ratio - 1 [-0.0126, +0.0805], dphase [-0.10, +0.26] deg
    g7 perturbed assembly (logged): amp ratio - 1 [-0.0386, +0.0024], dphase [-0.16, -0.01] deg
    g1 core sum dphi/phi0: OpenMC 0.00966 at  +85.15 deg (sigma 0.0069), FEMFFUSION 0.00973 at  +85.23 deg; ratio - 1 -0.0075, dphase -0.09 deg
    g2 core sum dphi/phi0: OpenMC 0.00961 at  +85.13 deg (sigma 0.0069), FEMFFUSION 0.00968 at  +85.21 deg; ratio - 1 -0.0076, dphase -0.09 deg
    g3 core sum dphi/phi0: OpenMC 0.00961 at  +85.13 deg (sigma 0.0069), FEMFFUSION 0.00968 at  +85.21 deg; ratio - 1 -0.0075, dphase -0.08 deg
    g4 core sum dphi/phi0: OpenMC 0.00962 at  +85.13 deg (sigma 0.0069), FEMFFUSION 0.00969 at  +85.22 deg; ratio - 1 -0.0074, dphase -0.09 deg
    g5 core sum dphi/phi0: OpenMC 0.00962 at  +85.13 deg (sigma 0.0069), FEMFFUSION 0.00970 at  +85.22 deg; ratio - 1 -0.0076, dphase -0.09 deg
    g6 core sum dphi/phi0: OpenMC 0.00964 at  +85.13 deg (sigma 0.0069), FEMFFUSION 0.00971 at  +85.22 deg; ratio - 1 -0.0078, dphase -0.09 deg
    g7 core sum dphi/phi0: OpenMC 0.00967 at  +85.14 deg (sigma 0.0069), FEMFFUSION 0.00975 at  +85.22 deg; ratio - 1 -0.0080, dphase -0.09 deg
    g1 static shape (fuel, sum-normalised): |OpenMC/FEMFFUSION - 1| median 0.0076, max 0.0803
    g7 static shape (fuel, sum-normalised): |OpenMC/FEMFFUSION - 1| median 0.0194, max 0.0744

  OpenMC run 'c102': 20 replicas, cancellation 102, keff 1.186489, T=50203s
    g1: PASS: evaluated amp 2017/2312, phase 2309/2312; max |amp ratio - 1| 0.0441, max |dphase| 0.64 deg; all cells: amp ratio - 1 median -0.0096 [5-95 % -0.0192, +0.0108], dphase median +0.05 deg; z hist (logged) mean -0.96 std 0.69
    g2: FAIL: evaluated amp 1972/2312, phase 2306/2312; max |amp ratio - 1| 0.0533, max |dphase| 0.76 deg; all cells: amp ratio - 1 median -0.0097 [5-95 % -0.0203, +0.0152], dphase median +0.05 deg; z hist (logged) mean -0.95 std 0.73
    g3: PASS: evaluated amp 1876/2312, phase 2289/2312; max |amp ratio - 1| 0.0476, max |dphase| 0.64 deg; all cells: amp ratio - 1 median -0.0096 [5-95 % -0.0226, +0.0233], dphase median +0.05 deg; z hist (logged) mean -0.90 std 0.79
    g4: PASS: evaluated amp 1752/2312, phase 2271/2312; max |amp ratio - 1| 0.0455, max |dphase| 1.06 deg; all cells: amp ratio - 1 median -0.0093 [5-95 % -0.0229, +0.0400], dphase median +0.05 deg; z hist (logged) mean -0.82 std 0.84
    g5: PASS: evaluated amp 1707/2312, phase 2262/2312; max |amp ratio - 1| 0.0404, max |dphase| 1.47 deg; all cells: amp ratio - 1 median -0.0093 [5-95 % -0.0239, +0.0450], dphase median +0.05 deg; z hist (logged) mean -0.81 std 0.86
    g6: PASS: evaluated amp 2080/2312, phase 2300/2312; max |amp ratio - 1| 0.0384, max |dphase| 0.65 deg; all cells: amp ratio - 1 median -0.0096 [5-95 % -0.0189, +0.0078], dphase median +0.05 deg; z hist (logged) mean -1.03 std 0.67
    g7: PASS: evaluated amp 2134/2312, phase 2306/2312; max |amp ratio - 1| 0.0348, max |dphase| 0.46 deg; all cells: amp ratio - 1 median -0.0097 [5-95 % -0.0174, +0.0035], dphase median +0.05 deg; z hist (logged) mean -1.06 std 0.63
    against ref1_FE2_SP3_nondiag (logged): g1 FAIL (max amp 0.358), g2 FAIL (max amp 0.252), g3 FAIL (max amp 0.186), g4 FAIL (max amp 0.094), g5 FAIL (max amp 0.060), g6 FAIL (max amp 0.292), g7 FAIL (max amp 0.267)
    g1 perturbed pin: OpenMC 0.01245 at  +86.35 deg, FEMFFUSION 0.01152 at  +85.97 deg
    g4 perturbed pin: OpenMC 0.01389 at  +86.73 deg, FEMFFUSION 0.01334 at  +86.53 deg
    g7 perturbed pin: OpenMC 0.02882 at  +88.43 deg, FEMFFUSION 0.02878 at  +88.39 deg
    g1 perturbed assembly (logged): amp ratio - 1 [-0.0147, +0.0807], dphase [+0.03, +0.37] deg
    g7 perturbed assembly (logged): amp ratio - 1 [-0.0363, +0.0014], dphase [-0.05, +0.10] deg
    g1 core sum dphi/phi0: OpenMC 0.00963 at  +85.28 deg (sigma 0.0076), FEMFFUSION 0.00973 at  +85.23 deg; ratio - 1 -0.0099, dphase +0.05 deg
    g2 core sum dphi/phi0: OpenMC 0.00958 at  +85.26 deg (sigma 0.0076), FEMFFUSION 0.00968 at  +85.21 deg; ratio - 1 -0.0099, dphase +0.05 deg
    g3 core sum dphi/phi0: OpenMC 0.00959 at  +85.27 deg (sigma 0.0076), FEMFFUSION 0.00968 at  +85.21 deg; ratio - 1 -0.0099, dphase +0.05 deg
    g4 core sum dphi/phi0: OpenMC 0.00960 at  +85.27 deg (sigma 0.0076), FEMFFUSION 0.00969 at  +85.22 deg; ratio - 1 -0.0099, dphase +0.05 deg
    g5 core sum dphi/phi0: OpenMC 0.00960 at  +85.27 deg (sigma 0.0076), FEMFFUSION 0.00970 at  +85.22 deg; ratio - 1 -0.0099, dphase +0.05 deg
    g6 core sum dphi/phi0: OpenMC 0.00962 at  +85.27 deg (sigma 0.0076), FEMFFUSION 0.00971 at  +85.22 deg; ratio - 1 -0.0100, dphase +0.05 deg
    g7 core sum dphi/phi0: OpenMC 0.00965 at  +85.27 deg (sigma 0.0076), FEMFFUSION 0.00975 at  +85.22 deg; ratio - 1 -0.0100, dphase +0.05 deg
    g1 static shape (fuel, sum-normalised): |OpenMC/FEMFFUSION - 1| median 0.0077, max 0.0873
    g7 static shape (fuel, sum-normalised): |OpenMC/FEMFFUSION - 1| median 0.0192, max 0.0726
  'c51' vs 'c102' g1 (logged): consistent, max|z| 3.60, Hotelling p 1.000
  'c51' vs 'c102' g7 (logged): consistent, max|z| 3.22, Hotelling p 1.000

Stage 9a: 12/15 gate checks passed
  FAILED: 9a 'c51' group 3
  FAILED: 9a 'c51' group 4
  FAILED: 9a 'c102' group 2
```
