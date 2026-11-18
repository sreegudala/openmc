# Stage 5 log (TACC job 3592996, commit 1d41001dc)

```
Test 5b: bias of mesh averaging, slab 4a at 1 Hz
  h=2.0 shift=0.0h: PASS paired: n=120 bad=0 max|t|=1.10; integral t Re=+0.33 Im=+0.28; Hotelling k=4 T2=2.4 p=0.733 | logged: paired t hist mean=+0.32 std=0.30, unpaired t hist mean=+0.31 std=0.30, on/off replica corr -0.11, max |on/off - 1| 9.84e-03
    (h=2 cm recorded as bias scale, not gated)
  h=2.0 shift=0.5h: PASS paired: n=120 bad=0 max|t|=1.31; integral t Re=+1.09 Im=+0.01; Hotelling k=5 T2=4.7 p=0.606 | logged: paired t hist mean=+0.53 std=0.60, unpaired t hist mean=+0.36 std=0.44, on/off replica corr 0.53, max |on/off - 1| 9.21e-03
    (h=2 cm recorded as bias scale, not gated)
  h=1.0 shift=0.0h: PASS paired: n=120 bad=0 max|t|=0.78; integral t Re=-0.33 Im=+0.53; Hotelling k=4 T2=2.7 p=0.693 | logged: paired t hist mean=+0.12 std=0.43, unpaired t hist mean=+0.13 std=0.40, on/off replica corr 0.29, max |on/off - 1| 2.96e-03
  h=1.0 shift=0.5h: PASS paired: n=120 bad=0 max|t|=1.01; integral t Re=+0.63 Im=+0.37; Hotelling k=4 T2=4.2 p=0.494 | logged: paired t hist mean=+0.50 std=0.27, unpaired t hist mean=+0.45 std=0.24, on/off replica corr 0.18, max |on/off - 1| 7.45e-03
  h=0.5 shift=0.0h: PASS paired: n=120 bad=0 max|t|=0.73; integral t Re=+0.51 Im=-0.24; Hotelling k=4 T2=2.3 p=0.740 | logged: paired t hist mean=+0.12 std=0.41, unpaired t hist mean=+0.12 std=0.38, on/off replica corr 0.09, max |on/off - 1| 5.53e-03
  h=0.5 shift=0.5h: PASS paired: n=120 bad=0 max|t|=0.93; integral t Re=-0.47 Im=+0.72; Hotelling k=4 T2=2.8 p=0.671 | logged: paired t hist mean=+0.12 std=0.60, unpaired t hist mean=+0.19 std=0.62, on/off replica corr 0.26, max |on/off - 1| 4.21e-03
Test 5d: FOM with and without cancellation, slab 4a
  0.1 Hz: median FOM ratio on/off Re=2.05 Im=3.66 (T off 3413s, on 3288s)
  1.0 Hz: median FOM ratio on/off Re=0.79 Im=4.15 (T off 3034s, on 3233s)
Test 5e: histogram scattering library (32 bins, mu_bar=0.3)
  hist h=1.0: PASS paired: n=120 bad=0 max|t|=1.30; integral t Re=-1.07 Im=+0.66; Hotelling k=5 T2=1.5 p=0.938 | logged: paired t hist mean=-0.15 std=0.85, unpaired t hist mean=-0.10 std=0.76, on/off replica corr 0.26, max |on/off - 1| 9.02e-03
  hist h=0.5: PASS paired: n=120 bad=0 max|t|=2.52; integral t Re=-2.43 Im=+1.07; Hotelling k=5 T2=6.9 p=0.407 | logged: paired t hist mean=-0.60 std=1.72, unpaired t hist mean=-0.52 std=1.49, on/off replica corr 0.24, max |on/off - 1| 1.81e-02

Stage 5 (5b,5d,5e): 7/8 gate checks passed
  FAILED: 5d 1.0 Hz FOM ratio > 1
```
