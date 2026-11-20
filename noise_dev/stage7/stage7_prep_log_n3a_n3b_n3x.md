# Stage 7 preparation log (N3a, N3b)

```
Test N3a: fission-free transposed library (copied step3a), report {'materials': 1, 'max_balance_error': 0.0, 'max_multiplicity': 2.4250000000000003, 'groups_absorption_zeroed': 1}
  N3a second-moment radius: forward 0.161, transposed 0.333
  n3a_100Hz tl: PASS n=6 bad=0 z_mean=-0.08 z_std=1.01 max|z|=1.86; amp: 3 evaluated, max rel err 1.54e-03; phase: 3 evaluated, max 0.010 deg; mean=[1.463383-0.06124j  1.317614-0.054989j 1.28608 -0.053913j] ref=[1.463709-0.061297j 1.317336-0.055208j 1.288062-0.053977j] corr=[ 0.15 -0.57 -0.35] FOM Re=[52908.34 32357.56 11859.13] Im=[650.09 560.69 342.91]
  n3a_100Hz col: PASS n=6 bad=0 z_mean=+0.01 z_std=1.01 max|z|=1.66; amp: 3 evaluated, max rel err 9.66e-04; phase: 3 evaluated, max 0.010 deg; mean=[1.463507-0.061162j 1.317305-0.055105j 1.286828-0.053701j] ref=[1.463709-0.061297j 1.317336-0.055208j 1.288062-0.053977j] corr=[-0.43 -0.43 -0.53] FOM Re=[146750.71 311724.61  24320.88] Im=[ 901.5  1134.87  309.55]
  n3a_100Hz generations/batch mean=0 max=0  T=123s
  n3a_3000Hz tl: PASS n=6 bad=0 z_mean=-0.64 z_std=0.56 max|z|=1.14; amp: 2 evaluated, max rel err 6.00e-03; phase: 1 evaluated, max 0.201 deg; mean=[0.567044-0.721188j 0.512002-0.649159j 0.492374-0.631408j] ref=[0.568686-0.714455j 0.511211-0.643492j 0.499911-0.629144j] corr=[0.44 0.22 0.76] FOM Re=[280.9  572.31 100.27] Im=[260.96 344.15 314.2 ]
  n3a_3000Hz col: PASS n=6 bad=0 z_mean=-0.90 z_std=0.57 max|z|=1.73; amp: 2 evaluated, max rel err 5.91e-03; phase: 2 evaluated, max 0.471 deg; mean=[0.566902-0.72254j  0.508893-0.651503j 0.497212-0.633189j] ref=[0.568686-0.714455j 0.511211-0.643492j 0.499911-0.629144j] corr=[0.73 0.57 0.86] FOM Re=[499.16 440.53 107.52] Im=[368.38 395.58 156.17]
  n3a_3000Hz generations/batch mean=0 max=0  T=50s
  n3a_10000Hz tl: PASS n=6 bad=0 z_mean=+0.49 z_std=0.92 max|z|=1.21; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[0.075731-0.314705j 0.072708-0.298962j 0.076848-0.289768j] ref=[0.079099-0.331246j 0.070252-0.298342j 0.068784-0.29169j ] corr=[-0.77 -0.25  0.21] FOM Re=[6.87 7.76 1.71] Im=[  6.68 186.91 264.36]
  n3a_10000Hz col: PASS n=6 bad=0 z_mean=+0.49 z_std=1.15 max|z|=1.43; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[0.075911-0.321755j 0.06628 -0.290326j 0.079458-0.288479j] ref=[0.079099-0.331246j 0.070252-0.298342j 0.068784-0.29169j ] corr=[-0.86 -0.9   0.79] FOM Re=[7.4  3.2  1.45] Im=[ 18.01  22.33 151.98]
  n3a_10000Hz generations/batch mean=0 max=0  T=78s
Test N3b: transposed library with fission and delayed data, report {'materials': 1, 'max_balance_error': 0.0, 'max_multiplicity': 2.4250000000000003, 'groups_absorption_zeroed': 1, 'fission_form': 'matrix', 'fissionable': ['m']}
  k_inf forward (ref_inf) 1.577161055718, transposed library 1.577161055718
  N3b second-moment radius: forward 0.161, transposed 0.333
  n3b_1Hz tl: PASS n=6 bad=0 z_mean=+0.45 z_std=1.79 max|z|=2.09; amp: 0 evaluated, max rel err nan; phase: 3 evaluated, max 0.292 deg; mean=[198.939001-26.017684j 180.153323-23.565582j 178.206615-23.30675j ] ref=[201.517312-27.399885j 182.502406-24.815171j 180.518359-24.546557j] corr=[-0.37 -0.38 -0.37] FOM Re=[6.55 6.54 6.41] Im=[1.22 1.22 1.21]
  n3b_1Hz col: PASS n=6 bad=0 z_mean=+0.45 z_std=1.79 max|z|=2.09; amp: 0 evaluated, max rel err nan; phase: 3 evaluated, max 0.293 deg; mean=[198.940817-26.016127j 180.167987-23.563331j 178.187787-23.303295j] ref=[201.517312-27.399885j 182.502406-24.815171j 180.518359-24.546557j] corr=[-0.38 -0.38 -0.37] FOM Re=[6.53 6.54 6.44] Im=[1.22 1.22 1.2 ]
  n3b_1Hz generations/batch mean=725 max=1956  T=1270s
  n3b_30Hz tl: PASS n=6 bad=0 z_mean=+0.85 z_std=0.51 max|z|=1.34; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[52.047992-86.773986j 47.137357-78.578169j 46.667231-77.692302j] ref=[50.444949-87.296825j 45.684257-79.060026j 45.187634-78.200453j] corr=[0.45 0.44 0.44] FOM Re=[0.48 0.48 0.47] Im=[1.   1.   0.98]
  n3b_30Hz col: PASS n=6 bad=0 z_mean=+0.85 z_std=0.51 max|z|=1.32; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[52.048657-86.769391j 47.14275 -78.584516j 46.636701-77.705399j] ref=[50.444949-87.296825j 45.684257-79.060026j 45.187634-78.200453j] corr=[0.44 0.44 0.44] FOM Re=[0.48 0.48 0.48] Im=[1.   1.   0.98]
  n3b_30Hz generations/batch mean=1718 max=4295  T=3761s
  n3b_100Hz tl: PASS n=6 bad=0 z_mean=+0.53 z_std=0.07 max|z|=0.63; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[6.362519-33.482803j 5.745954-30.320568j 5.622625-29.983484j] ref=[5.89605 -33.961086j 5.338732-30.756722j 5.280784-30.422318j] corr=[0.03 0.04 0.04] FOM Re=[0. 0. 0.] Im=[0.06 0.06 0.06]
  n3b_100Hz col: PASS n=6 bad=0 z_mean=+0.54 z_std=0.07 max|z|=0.62; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[6.351784-33.478378j 5.75733 -30.320103j 5.668368-29.989925j] ref=[5.89605 -33.961086j 5.338732-30.756722j 5.280784-30.422318j] corr=[0.03 0.03 0.01] FOM Re=[0. 0. 0.] Im=[0.06 0.06 0.06]
  n3b_100Hz generations/batch mean=5055 max=19358  T=18086s
Test N3x (informational, not gated): infinite-variance transpose, report {'materials': 1, 'max_balance_error': 0.0, 'max_multiplicity': 2.4, 'groups_absorption_zeroed': 1}
  N3x second-moment radius: forward 0.769, transposed 2.047
  n3x_100Hz tl: FAIL n=6 bad=5 z_mean=+3.69 z_std=8.86 max|z|=13.95; amp: 0 evaluated, max rel err nan; phase: 3 evaluated, max 0.727 deg; mean=[8.896032-0.137352j 7.930408-0.126096j 7.826327-0.135351j] ref=[9.840695-0.276878j 9.083706-0.256018j 8.841475-0.249137j] corr=[-0.2  -0.04  0.15] FOM Re=[ 5.01 13.11  9.14] Im=[1.66 1.39 0.71]
  n3x_100Hz col: FAIL n=6 bad=6 z_mean=+2.24 z_std=6.60 max|z|=8.58; amp: 0 evaluated, max rel err nan; phase: 3 evaluated, max 0.651 deg; mean=[8.758006-0.146799j 8.17162 -0.141681j 7.927608-0.140689j] ref=[9.840695-0.276878j 9.083706-0.256018j 8.841475-0.249137j] corr=[-0.14 -0.12 -0.02] FOM Re=[ 6.46  6.43 15.86] Im=[0.81 1.   0.86]
  n3x_100Hz generations/batch mean=0 max=0  T=113s
Pooled z: n=72 mean=+0.25 std=1.06 |z|>=3: 0

Stage 7 prep (n3a,n3b,n3x): 6/6 checks passed
```
