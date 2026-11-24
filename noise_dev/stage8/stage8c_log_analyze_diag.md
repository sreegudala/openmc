# Stage 8c log

```
Stage 8c n=2
  Abeille batches [600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600]
  c=170 fast vs Abeille: PASS n=38070 bad=0 max|t|=2.27; integral t Re=+0.62 Im=+0.43; Hotelling k=2 p=0.670; modulus OpenMC/Abeille - 1 median -0.386
  c=170 thermal vs Abeille: PASS n=38004 bad=8 max|t|=4.98; integral t Re=+0.61 Im=+0.43; Hotelling k=2 p=0.672; modulus OpenMC/Abeille - 1 median -0.382
  c=340 fast vs Abeille: PASS n=38070 bad=0 max|t|=2.01; integral t Re=+0.86 Im=-0.60; Hotelling k=2 p=0.679; modulus OpenMC/Abeille - 1 median -0.767
  c=340 thermal vs Abeille: PASS n=38004 bad=7 max|t|=3.76; integral t Re=+0.85 Im=-0.60; Hotelling k=2 p=0.684; modulus OpenMC/Abeille - 1 median -0.764
  fast cancellation 170 vs 340: PASS integral t Re=-0.44 Im=+0.91; Hotelling p=0.666
  thermal cancellation 170 vs 340: PASS integral t Re=-0.44 Im=+0.91; Hotelling p=0.666
  c=170 detector vs contraction: PASS paired t Re=-0.22 Im=-0.11; forward -6.45082e-04+2.93811e-05j contraction -5.61268e-04+3.58857e-05j (P0 -1.17330e-04+2.78655e-07j, P1 -1.00557e-05-1.67504e-07j, fission -4.33883e-04+3.57746e-05j)
    pin detector (thermal), forward det tally: -3.57161e-01+1.60932e-02j +- 2.16e-01
    pin detector (thermal), forward map sum: -3.57161e-01+1.60932e-02j +- 2.16e-01
    pin detector (thermal), contraction: -3.10812e-01+1.98733e-02j +- 2.88e-03
    pin detector (thermal), Abeille map sum: -6.84185e-01-8.74766e-03j +- 1.79e-01
  c=340 detector vs contraction: PASS paired t Re=+0.33 Im=-0.99; forward -3.33357e-04-1.23827e-04j contraction -5.62371e-04+3.59301e-05j (P0 -1.17555e-04+3.12960e-07j, P1 -1.08650e-05-1.31656e-07j, fission -4.33951e-04+3.57488e-05j)
    pin detector (thermal), forward det tally: -1.84867e-01-6.85447e-02j +- 3.92e-01
    pin detector (thermal), forward map sum: -1.84867e-01-6.85447e-02j +- 3.92e-01
    pin detector (thermal), contraction: -3.11393e-01+1.98957e-02j +- 3.21e-03
    pin detector (thermal), Abeille map sum: -6.84185e-01-8.74766e-03j +- 1.79e-01

Stage 8c: 8/8 gate checks passed
```
