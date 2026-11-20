# Stage 6c log

```
Stage 6c: OpenMC 20 replicas, Abeille 20 replicas (batches [223, 226, 230, 224, 228, 226, 227, 225, 230, 226, 225, 225, 224, 223, 219, 229, 227, 229, 225, 228]); Abeille static normalisation rel sigma 8.25e-04; Abeille static k 0.99912 +- 0.00006
  static phi_0 / phi_0(centre, thermal): max |OpenMC / Abeille - 1| = 9.074e-03 (fast and thermal, OpenMC replica mean)
  fast: PASS two-sample STD-C: n=38088 bad=0 max|z|=1.16, z hist mean -0.76 std 0.19; integral z Re=-0.94 Im=-0.58; Hotelling k=2 T2=0.9 p=0.635; per-cell modulus OpenMC/Abeille - 1: median -0.014, 5-95 % [-0.016, -0.013]
    fast assembly integral: OpenMC 3.55043e+02+6.72526e+03j Abeille 3.76535e+02+6.82317e+03j; |modulus ratio - 1| = 0.0145 (sigma OpenMC 0.0116, Abeille 0.0223), phase diff +0.14 deg -> PASS (2 %)
    fast perturbed pin: OpenMC 3.04210e-01+6.01331e+00j Abeille 3.22167e-01+6.09983e+00j; |modulus ratio - 1| = 0.0143 (sigma OpenMC 0.0110, Abeille 0.0213), phase diff +0.13 deg -> PASS (2 %)
  thermal: PASS two-sample STD-C: n=38088 bad=0 max|z|=1.37, z hist mean -0.76 std 0.20; integral z Re=-0.94 Im=-0.59; Hotelling k=2 T2=0.9 p=0.634; per-cell modulus OpenMC/Abeille - 1: median -0.015, 5-95 % [-0.018, -0.012]
    thermal assembly integral: OpenMC 1.05703e+02+2.00146e+03j Abeille 1.12090e+02+2.03095e+03j; |modulus ratio - 1| = 0.0146 (sigma OpenMC 0.0116, Abeille 0.0222), phase diff +0.14 deg -> PASS (2 %)
    thermal perturbed pin: OpenMC 8.44434e-02+1.93236e+00j Abeille 8.94371e-02+1.95395e+00j; |modulus ratio - 1| = 0.0111 (sigma OpenMC 0.0096, Abeille 0.0184), phase diff +0.12 deg -> PASS (2 %)

Stage 6c: 6/6 gate checks passed
```
