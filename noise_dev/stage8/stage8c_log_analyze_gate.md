# Stage 8c log

```
Stage 8c n=1
  Abeille batches [600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600]
  c=170 fast vs Abeille: PASS n=14750 bad=19 max|t|=6.26; integral t Re=-0.07 Im=-1.37; Hotelling k=2 p=0.329; modulus OpenMC/Abeille - 1 median -0.423
  c=170 thermal vs Abeille: FAIL n=462 bad=76 max|t|=14.41; integral t Re=-0.05 Im=-1.86; Hotelling k=4 p=0.000; modulus OpenMC/Abeille - 1 median -0.874
  c=340 fast vs Abeille: PASS n=14750 bad=11 max|t|=4.56; integral t Re=-1.18 Im=-0.41; Hotelling k=2 p=0.513; modulus OpenMC/Abeille - 1 median -0.519
  c=340 thermal vs Abeille: FAIL n=462 bad=45 max|t|=11.20; integral t Re=-1.18 Im=-0.82; Hotelling k=3 p=0.000; modulus OpenMC/Abeille - 1 median -0.656
  fast cancellation 170 vs 340: PASS integral t Re=+1.28 Im=-0.87; Hotelling p=0.193
  thermal cancellation 170 vs 340: PASS integral t Re=+1.31 Im=-0.90; Hotelling p=0.178
  pin_top c=170 forward vs contraction: PASS paired t Re=+0.61 Im=-1.12; P0 +6.13168e-04+6.88477e-01j, P1 +2.21490e-04+3.87509e-02j, fission +5.88745e-04+9.98361e-03j
    pin_top (thermal, normalised), forward map sum: +1.85891e-02+5.85415e-01j +- 1.37e-01
    pin_top (thermal, normalised), contraction: +1.42340e-03+7.37211e-01j +- 2.66e-03
    pin_top (thermal, normalised), Abeille map sum: +2.35735e-02+2.28802e+00j +- 2.03e-01
  pin_top c=340 forward vs contraction: PASS paired t Re=-1.21 Im=+0.23; P0 +6.69704e-04+6.90377e-01j, P1 +2.25233e-04+3.86372e-02j, fission +5.53606e-04+8.93078e-03j
    pin_top (thermal, normalised), forward map sum: -3.90851e-02+7.85481e-01j +- 2.10e-01
    pin_top (thermal, normalised), contraction: +1.44854e-03+7.37945e-01j +- 2.20e-03
    pin_top (thermal, normalised), Abeille map sum: +2.35735e-02+2.28802e+00j +- 2.03e-01
  pin_bottom c=170 forward vs contraction: PASS paired t Re=+0.60 Im=-1.13; P0 +5.22029e-04-6.56663e-01j, P1 +1.90687e-04-1.77219e-02j, fission +1.23471e-03+3.70536e-02j
    pin_bottom (thermal, normalised), forward map sum: +1.89811e-02-7.93010e-01j +- 1.40e-01
    pin_bottom (thermal, normalised), contraction: +1.94742e-03-6.37331e-01j +- 2.23e-03
    pin_bottom (thermal, normalised), Abeille map sum: +2.45283e-02-1.78004e+00j +- 2.10e-01
  pin_bottom c=340 forward vs contraction: PASS paired t Re=-1.23 Im=+0.33; P0 +6.10890e-04-6.55445e-01j, P1 +1.87854e-04-1.78459e-02j, fission +1.14716e-03+3.59616e-02j
    pin_bottom (thermal, normalised), forward map sum: -3.86035e-02-5.67450e-01j +- 2.12e-01
    pin_bottom (thermal, normalised), contraction: +1.94591e-03-6.37330e-01j +- 2.54e-03
    pin_bottom (thermal, normalised), Abeille map sum: +2.45283e-02-1.78004e+00j +- 2.10e-01
Stage 8c n=2
  Abeille batches [600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600, 600]
  c=170 fast vs Abeille: PASS n=38070 bad=0 max|t|=2.58; integral t Re=+0.13 Im=+0.91; Hotelling k=2 p=0.596; modulus OpenMC/Abeille - 1 median -0.051
  c=170 thermal vs Abeille: PASS n=38004 bad=19 max|t|=6.12; integral t Re=+0.11 Im=+0.91; Hotelling k=2 p=0.599; modulus OpenMC/Abeille - 1 median -0.045
  c=340 fast vs Abeille: PASS n=38070 bad=0 max|t|=2.89; integral t Re=+0.58 Im=+1.52; Hotelling k=2 p=0.183; modulus OpenMC/Abeille - 1 median -0.261
  c=340 thermal vs Abeille: PASS n=38004 bad=20 max|t|=6.21; integral t Re=+0.58 Im=+1.52; Hotelling k=2 p=0.183; modulus OpenMC/Abeille - 1 median -0.260
  fast cancellation 170 vs 340: PASS integral t Re=-0.60 Im=-1.05; Hotelling p=0.255
  thermal cancellation 170 vs 340: PASS integral t Re=-0.61 Im=-1.05; Hotelling p=0.607
  pin c=170 forward vs contraction: PASS paired t Re=-1.84 Im=+0.87; P0 -6.75752e-02+3.42694e-04j, P1 -5.75106e-03-7.91469e-05j, fission -2.36720e-01+1.94970e-02j
    pin (thermal, normalised), forward map sum: -4.92556e-01+3.93725e-02j +- 1.01e-01
    pin (thermal, normalised), contraction: -3.10046e-01+1.97605e-02j +- 2.43e-03
    pin (thermal, normalised), Abeille map sum: -6.84185e-01-8.74766e-03j +- 1.79e-01
  pin c=340 forward vs contraction: PASS paired t Re=-0.65 Im=+2.03; P0 -6.89068e-02+4.17491e-04j, P1 -5.65409e-03-7.97607e-05j, fission -2.34948e-01+1.93856e-02j
    pin (thermal, normalised), forward map sum: -3.94193e-01+7.54864e-02j +- 1.33e-01
    pin (thermal, normalised), contraction: -3.09509e-01+1.97234e-02j +- 2.50e-03
    pin (thermal, normalised), Abeille map sum: -6.84185e-01-8.74766e-03j +- 1.79e-01

Stage 8c: 16/18 gate checks passed
  FAILED: 8c n=1 c=170 thermal vs Abeille
  FAILED: 8c n=1 c=340 thermal vs Abeille
```
