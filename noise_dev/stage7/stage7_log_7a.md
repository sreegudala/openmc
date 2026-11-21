# Stage 7 log

```
Test 7a: infinite medium, adjoint on the N2 library vs ref_inf adjoint (delayed neutrons on)
  n3: k_inf 1.5771610557 (transposed library 1.5771610557); max multiplicity 2.425; second-moment radius 0.333
  7a_n3_0.1Hz tl: PASS n=6 bad=0 z_mean=-0.44 z_std=1.32 max|z|=1.65; amp: 3 evaluated, max rel err 5.17e-03; phase: 3 evaluated, max 0.245 deg; mean=[240.849905-91.049257j 218.127488-82.461189j 215.769291-81.57759j ] ref=[240.011109-89.563034j 217.365521-81.115612j 215.005245-80.240325j] corr=[-0.35 -0.34 -0.34] FOM Re=[7.15 7.14 7.08] Im=[1.54 1.54 1.52]
  7a_n3_0.1Hz col: PASS n=6 bad=0 z_mean=-0.44 z_std=1.32 max|z|=1.64; amp: 3 evaluated, max rel err 5.16e-03; phase: 3 evaluated, max 0.245 deg; mean=[240.852042-91.046721j 218.126316-82.458581j 215.767102-81.577424j] ref=[240.011109-89.563034j 217.365521-81.115612j 215.005245-80.240325j] corr=[-0.34 -0.34 -0.34] FOM Re=[7.15 7.15 7.07] Im=[1.54 1.54 1.52]
  7a_n3_0.1Hz generations/batch mean=1327 max=2963  T=6603s
  7a_n3_1Hz tl: PASS n=6 bad=0 z_mean=-0.36 z_std=0.75 max|z|=1.05; amp: 3 evaluated, max rel err 1.95e-03; phase: 3 evaluated, max 0.086 deg; mean=[201.869821-27.751051j 182.817948-25.131999j 180.830173-24.866275j] ref=[201.517312-27.399885j 182.502406-24.815171j 180.518359-24.546557j] corr=[-0.43 -0.43 -0.43] FOM Re=[8.35 8.34 8.27] Im=[1.56 1.55 1.53]
  7a_n3_1Hz col: PASS n=6 bad=0 z_mean=-0.36 z_std=0.75 max|z|=1.04; amp: 3 evaluated, max rel err 1.95e-03; phase: 3 evaluated, max 0.085 deg; mean=[201.867522-27.750705j 182.819119-25.132798j 180.832886-24.863751j] ref=[201.517312-27.399885j 182.502406-24.815171j 180.518359-24.546557j] corr=[-0.43 -0.43 -0.42] FOM Re=[8.35 8.37 8.25] Im=[1.55 1.56 1.54]
  7a_n3_1Hz generations/batch mean=953 max=1778  T=4359s
  7a_n3_10Hz tl: PASS n=6 bad=0 z_mean=+0.27 z_std=1.04 max|z|=1.22; amp: 3 evaluated, max rel err 5.47e-03; phase: 3 evaluated, max 0.064 deg; mean=[151.173399-88.105933j 136.907316-79.789062j 135.420026-78.928652j] ref=[150.255497-87.79376j  136.077265-79.510114j 134.597806-78.645732j] corr=[-0.13 -0.13 -0.13] FOM Re=[7.4  7.37 7.2 ] Im=[6.83 6.84 6.75]
  7a_n3_10Hz col: PASS n=6 bad=0 z_mean=+0.27 z_std=1.05 max|z|=1.23; amp: 3 evaluated, max rel err 5.49e-03; phase: 3 evaluated, max 0.063 deg; mean=[151.177793-88.109674j 136.910234-79.796222j 135.422878-78.927523j] ref=[150.255497-87.79376j  136.077265-79.510114j 134.597806-78.645732j] corr=[-0.13 -0.13 -0.13] FOM Re=[7.37 7.38 7.26] Im=[6.83 6.83 6.73]
  7a_n3_10Hz generations/batch mean=1102 max=2664  T=5473s
  7a_n3_30Hz tl: PASS n=6 bad=0 z_mean=+1.28 z_std=0.43 max|z|=1.69; amp: 3 evaluated, max rel err 5.14e-03; phase: 3 evaluated, max 0.462 deg; mean=[50.880627-86.458838j 46.082166-78.308757j 45.581211-77.43384j ] ref=[50.444949-87.296825j 45.684257-79.060026j 45.187634-78.200453j] corr=[-0.  0.  0.] FOM Re=[0.75 0.75 0.74] Im=[2.12 2.14 2.07]
  7a_n3_30Hz col: PASS n=6 bad=0 z_mean=+1.27 z_std=0.43 max|z|=1.68; amp: 3 evaluated, max rel err 5.10e-03; phase: 3 evaluated, max 0.456 deg; mean=[50.883718-86.460988j 46.083352-78.304865j 45.570621-77.443659j] ref=[50.444949-87.296825j 45.684257-79.060026j 45.187634-78.200453j] corr=[ 0. -0. -0.] FOM Re=[0.75 0.75 0.75] Im=[2.13 2.13 2.1 ]
  7a_n3_30Hz generations/batch mean=2413 max=5871  T=14018s
  7a_n3_100Hz tl: PASS n=6 bad=0 z_mean=+1.75 z_std=0.57 max|z|=2.33; amp: 3 evaluated, max rel err 3.20e-03; phase: 3 evaluated, max 0.670 deg; mean=[6.262684-33.795016j 5.683653-30.614347j 5.610116-30.262976j] ref=[5.89605 -33.961086j 5.338732-30.756722j 5.280784-30.422318j] corr=[-0.   -0.01 -0.01] FOM Re=[0.13 0.13 0.12] Im=[5.43 5.38 5.32]
  7a_n3_100Hz col: PASS n=6 bad=0 z_mean=+1.75 z_std=0.57 max|z|=2.30; amp: 3 evaluated, max rel err 3.13e-03; phase: 3 evaluated, max 0.661 deg; mean=[6.269162-33.798219j 5.677661-30.610676j 5.611139-30.264987j] ref=[5.89605 -33.961086j 5.338732-30.756722j 5.280784-30.422318j] corr=[-0.01 -0.01  0.  ] FOM Re=[0.13 0.13 0.13] Im=[5.39 5.4  5.28]
  7a_n3_100Hz generations/batch mean=1803 max=4146  T=11608s
  cortex: k_inf 0.6178836947 (transposed library 0.6178836947); max multiplicity 1.000; second-moment radius 0.864
  7a_cortex_0.1Hz tl: PASS n=4 bad=0 z_mean=+0.14 z_std=0.05 max|z|=0.18; amp: 0 evaluated, max rel err nan; phase: 2 evaluated, max 0.056 deg; mean=[114.339929-14.647377j  38.494041 -5.191979j] ref=[114.227005-14.738543j  38.45367  -5.224479j] corr=[-0.17 -0.17] FOM Re=[0.65 0.58] Im=[0.05 0.05]
  7a_cortex_0.1Hz col: PASS n=4 bad=0 z_mean=+0.14 z_std=0.05 max|z|=0.18; amp: 0 evaluated, max rel err nan; phase: 2 evaluated, max 0.055 deg; mean=[114.340227-14.647475j  38.493711 -5.191959j] ref=[114.227005-14.738543j  38.45367  -5.224479j] corr=[-0.17 -0.17] FOM Re=[0.65 0.58] Im=[0.05 0.05]
  7a_cortex_0.1Hz generations/batch mean=985 max=2644  T=15990s
  7a_cortex_1Hz tl: PASS n=4 bad=0 z_mean=+0.04 z_std=0.70 max|z|=0.65; amp: 0 evaluated, max rel err nan; phase: 2 evaluated, max 0.040 deg; mean=[114.872541-2.016196j  38.682803-0.711887j] ref=[114.225139-1.928826j  38.45302 -0.680719j] corr=[-0.19 -0.19] FOM Re=[0.93 0.84] Im=[0.01 0.01]
  7a_cortex_1Hz col: PASS n=4 bad=0 z_mean=+0.04 z_std=0.70 max|z|=0.65; amp: 0 evaluated, max rel err nan; phase: 2 evaluated, max 0.040 deg; mean=[114.872877-2.016163j  38.683321-0.711888j] ref=[114.225139-1.928826j  38.45302 -0.680719j] corr=[-0.19 -0.19] FOM Re=[0.94 0.83] Im=[0.01 0.01]
  7a_cortex_1Hz generations/batch mean=884 max=2250  T=14199s
  7a_cortex_10Hz tl: PASS n=4 bad=0 z_mean=+0.20 z_std=0.08 max|z|=0.27; amp: 0 evaluated, max rel err nan; phase: 2 evaluated, max 0.042 deg; mean=[114.159703-4.661352j  38.431643-1.621892j] ref=[114.038933-4.734969j  38.388215-1.648127j] corr=[-0.08 -0.08] FOM Re=[0.98 0.88] Im=[0.02 0.02]
  7a_cortex_10Hz col: PASS n=4 bad=0 z_mean=+0.20 z_std=0.08 max|z|=0.27; amp: 0 evaluated, max rel err nan; phase: 2 evaluated, max 0.042 deg; mean=[114.160445-4.661583j  38.431935-1.621874j] ref=[114.038933-4.734969j  38.388215-1.648127j] corr=[-0.08 -0.08] FOM Re=[0.98 0.88] Im=[0.02 0.02]
  7a_cortex_10Hz generations/batch mean=887 max=2031  T=14189s
  7a_cortex_30Hz tl: PASS n=4 bad=0 z_mean=+0.25 z_std=0.61 max|z|=0.77; amp: 0 evaluated, max rel err nan; phase: 2 evaluated, max 0.094 deg; mean=[112.290056-13.420528j  37.777932 -4.667324j] ref=[112.55664 -13.630753j  37.872335 -4.742024j] corr=[-0.53 -0.53] FOM Re=[0.96 0.85] Im=[0.17 0.16]
  7a_cortex_30Hz col: PASS n=4 bad=0 z_mean=+0.25 z_std=0.61 max|z|=0.77; amp: 0 evaluated, max rel err nan; phase: 2 evaluated, max 0.094 deg; mean=[112.290323-13.420652j  37.777341 -4.667393j] ref=[112.55664 -13.630753j  37.872335 -4.742024j] corr=[-0.53 -0.53] FOM Re=[0.96 0.85] Im=[0.17 0.16]
  7a_cortex_30Hz generations/batch mean=896 max=1895  T=14520s
  7a_cortex_100Hz tl: PASS n=4 bad=0 z_mean=+0.39 z_std=0.79 max|z|=1.08; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[99.179957-39.601449j 33.226289-13.773258j] ref=[98.098427-39.358462j 32.840465-13.690668j] corr=[-0.91 -0.91] FOM Re=[0.68 0.59] Im=[0.16 0.16]
  7a_cortex_100Hz col: PASS n=4 bad=0 z_mean=+0.39 z_std=0.79 max|z|=1.08; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[99.179384-39.601366j 33.22641 -13.773384j] ref=[98.098427-39.358462j 32.840465-13.690668j] corr=[-0.91 -0.91] FOM Re=[0.68 0.59] Im=[0.16 0.16]
  7a_cortex_100Hz generations/batch mean=884 max=2554  T=14350s
  7a pooled z: n=100 mean=+0.38 std=0.99 |z|>=3: 0

Stage 7 (7a): 10/10 gate checks passed
```
