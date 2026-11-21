# Stage 7 preparation log (N3a, N3b)

```
Test N3b: transposed library with fission and delayed data, report {'materials': 1, 'max_balance_error': 0.0, 'max_multiplicity': 2.4250000000000003, 'groups_absorption_zeroed': 1, 'fission_form': 'matrix', 'fissionable': ['m']}
  k_inf forward (ref_inf) 1.577161055718, transposed library 1.577161055718
  N3b second-moment radius: forward 0.161, transposed 0.333
  n3b_1Hz tl: PASS n=6 bad=0 z_mean=-0.09 z_std=1.05 max|z|=1.05; amp: 3 evaluated, max rel err 3.96e-03; phase: 3 evaluated, max 0.152 deg; mean=[202.230447-28.041045j 183.157945-25.399042j 181.165098-25.123765j] ref=[201.517312-27.399885j 182.502406-24.815171j 180.518359-24.546557j] corr=[-0.3 -0.3 -0.3] FOM Re=[13.3  13.27 13.08] Im=[0.47 0.47 0.47]
  n3b_1Hz col: PASS n=6 bad=0 z_mean=-0.09 z_std=1.05 max|z|=1.05; amp: 3 evaluated, max rel err 3.99e-03; phase: 3 evaluated, max 0.152 deg; mean=[202.239035-28.041547j 183.155894-25.396243j 181.172164-25.124928j] ref=[201.517312-27.399885j 182.502406-24.815171j 180.518359-24.546557j] corr=[-0.3  -0.3  -0.29] FOM Re=[13.27 13.29 13.08] Im=[0.47 0.47 0.46]
  n3b_1Hz generations/batch mean=931 max=1763  T=4458s
  n3b_30Hz tl: PASS n=6 bad=0 z_mean=-1.24 z_std=0.76 max|z|=1.95; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[50.065038-88.83648j  45.339049-80.458477j 44.85064 -79.573431j] ref=[50.444949-87.296825j 45.684257-79.060026j 45.187634-78.200453j] corr=[0.61 0.62 0.61] FOM Re=[0.36 0.36 0.36] Im=[0.87 0.88 0.85]
  n3b_30Hz col: PASS n=6 bad=0 z_mean=-1.24 z_std=0.76 max|z|=1.94; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[50.060727-88.834087j 45.334939-80.449711j 44.859047-79.584019j] ref=[50.444949-87.296825j 45.684257-79.060026j 45.187634-78.200453j] corr=[0.61 0.61 0.61] FOM Re=[0.36 0.36 0.36] Im=[0.88 0.88 0.85]
  n3b_30Hz generations/batch mean=2440 max=5283  T=14335s
  n3b_100Hz tl: PASS n=6 bad=0 z_mean=-0.98 z_std=0.88 max|z|=1.89; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[5.824648-34.715342j 5.2439  -31.444478j 5.233205-31.15452j ] ref=[5.89605 -33.961086j 5.338732-30.756722j 5.280784-30.422318j] corr=[0.06 0.06 0.1 ] FOM Re=[0. 0. 0.] Im=[0.09 0.09 0.09]
  n3b_100Hz col: PASS n=6 bad=0 z_mean=-0.99 z_std=0.87 max|z|=1.79; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg; mean=[5.803709-34.730338j 5.249692-31.457972j 5.219479-31.123315j] ref=[5.89605 -33.961086j 5.338732-30.756722j 5.280784-30.422318j] corr=[0.07 0.07 0.07] FOM Re=[0. 0. 0.] Im=[0.09 0.1  0.09]
  n3b_100Hz generations/batch mean=7672 max=23231  T=67895s
Test N3f: forward control on the N3 material (forward library)
  n3f_1Hz tl: PASS n=6 bad=0 z_mean=+0.44 z_std=0.55 max|z|=0.96; amp: 3 evaluated, max rel err 1.04e-04; phase: 3 evaluated, max 0.018 deg; mean=[201.50575 -27.335444j  64.454126 -8.778411j 601.140709-81.890053j] ref=[201.517312-27.399885j  64.457122 -8.799502j 601.178097-82.081571j] corr=[-0.44 -0.45 -0.44] FOM Re=[158.89 158.96 158.82] Im=[21.64 21.81 21.91]
  n3f_1Hz col: PASS n=6 bad=0 z_mean=+0.44 z_std=0.55 max|z|=0.95; amp: 3 evaluated, max rel err 1.02e-04; phase: 3 evaluated, max 0.018 deg; mean=[201.505643-27.335196j  64.45327  -8.778793j 601.143289-81.889317j] ref=[201.517312-27.399885j  64.457122 -8.799502j 601.178097-82.081571j] corr=[-0.45 -0.45 -0.44] FOM Re=[158.75 158.14 157.82] Im=[21.74 21.88 21.94]
  n3f_1Hz generations/batch mean=1394 max=2319  T=7343s
  n3f_30Hz tl: PASS n=6 bad=0 z_mean=-0.40 z_std=0.94 max|z|=1.28; amp: 3 evaluated, max rel err 1.28e-03; phase: 3 evaluated, max 0.076 deg; mean=[ 50.269959 -87.24854j   16.076429 -27.908243j 149.956172-260.291692j] ref=[ 50.444949 -87.296825j  16.133542 -27.9216j   150.476169-260.417133j] corr=[-0.33 -0.36 -0.34] FOM Re=[5.77 5.75 5.67] Im=[36.55 35.94 36.07]
  n3f_30Hz col: PASS n=6 bad=0 z_mean=-0.39 z_std=0.94 max|z|=1.24; amp: 3 evaluated, max rel err 1.26e-03; phase: 3 evaluated, max 0.073 deg; mean=[ 50.271338 -87.250767j  16.077768 -27.90754j  149.955588-260.279946j] ref=[ 50.444949 -87.296825j  16.133542 -27.9216j   150.476169-260.417133j] corr=[-0.34 -0.34 -0.34] FOM Re=[5.74 5.7  5.69] Im=[36.12 35.77 36.16]
  n3f_30Hz generations/batch mean=3809 max=6543  T=22430s
  n3f_100Hz tl: PASS n=6 bad=0 z_mean=+0.82 z_std=0.47 max|z|=1.27; amp: 3 evaluated, max rel err 3.81e-03; phase: 3 evaluated, max 0.122 deg; mean=[ 5.944566 -33.820986j  1.898644 -10.817258j 17.738741-100.909171j] ref=[ 5.89605  -33.961086j  1.885367 -10.862228j 17.587535-101.308992j] corr=[0.09 0.12 0.07] FOM Re=[0.02 0.02 0.02] Im=[0.95 0.94 0.93]
  n3f_100Hz col: PASS n=6 bad=0 z_mean=+0.81 z_std=0.44 max|z|=1.23; amp: 3 evaluated, max rel err 3.69e-03; phase: 3 evaluated, max 0.120 deg; mean=[ 5.946387 -33.829479j  1.900203 -10.81833j  17.733396-100.913772j] ref=[ 5.89605  -33.961086j  1.885367 -10.862228j 17.587535-101.308992j] corr=[0.09 0.1  0.09] FOM Re=[0.03 0.02 0.03] Im=[0.94 0.93 0.93]
  n3f_100Hz generations/batch mean=12591 max=26980  T=99729s
Pooled z: n=72 mean=-0.24 std=1.04 |z|>=3: 0

Stage 7 prep (n3b,n3f): 6/6 checks passed
```
