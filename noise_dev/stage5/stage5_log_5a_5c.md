# Stage 5 log (TACC job 3592995, commit 1d41001dc; copied from slurm output)

```
Test 5a: TRIGA-like infinite medium with cancellation (4x4x4)
  5a_0.1Hz tl: PASS n=2 bad=0 z_mean=+0.17 z_std=3.04 max|z|=2.32; amp: 1 evaluated, max rel err 8.55e-03; phase: 1 evaluated, max 0.251 deg; mean=[1656.472286-606.384413j] ref=[1668.060305-618.9135j] corr=[-0.82] FOM Re=[11.93] Im=[1.86]
  5a_0.1Hz col: PASS n=2 bad=0 z_mean=+0.17 z_std=3.04 max|z|=2.32; amp: 1 evaluated, max rel err 8.55e-03; phase: 1 evaluated, max 0.251 deg; mean=[1656.473652-606.387102j] ref=[1668.060305-618.9135j] corr=[-0.82] FOM Re=[11.93] Im=[1.86]
  5a_0.1Hz generations/batch mean=1120 max=1933  T=6766s
    measured ratio 0.99465; physical K 0.99473; K~cb 0.99579
  5a_1Hz tl: PASS n=2 bad=0 z_mean=+0.13 z_std=3.36 max|z|=2.51; amp: 1 evaluated, max rel err 8.59e-03; phase: 1 evaluated, max 0.048 deg; mean=[1392.325743-165.538833j] ref=[1404.248942-168.146236j] corr=[-0.82] FOM Re=[13.52] Im=[4.99]
  5a_1Hz col: PASS n=2 bad=0 z_mean=+0.13 z_std=3.36 max|z|=2.51; amp: 1 evaluated, max rel err 8.58e-03; phase: 1 evaluated, max 0.048 deg; mean=[1392.333415-165.539557j] ref=[1404.248942-168.146236j] corr=[-0.82] FOM Re=[13.54] Im=[5.]
  5a_1Hz generations/batch mean=879 max=2073  T=5074s
    measured ratio 0.99293; physical K 0.99298; K~cb 0.99354
  5a_10Hz tl: PASS n=2 bad=0 z_mean=+0.11 z_std=0.13 max|z|=0.20; amp: 1 evaluated, max rel err 1.90e-04; phase: 1 evaluated, max 0.037 deg; mean=[1184.999896-502.544912j] ref=[1184.900447-503.404938j] corr=[-0.91] FOM Re=[10.87] Im=[2.37]
  5a_10Hz col: PASS n=2 bad=0 z_mean=+0.11 z_std=0.13 max|z|=0.20; amp: 1 evaluated, max rel err 1.92e-04; phase: 1 evaluated, max 0.037 deg; mean=[1184.996645-502.544998j] ref=[1184.900447-503.404938j] corr=[-0.91] FOM Re=[10.88] Im=[2.37]
  5a_10Hz generations/batch mean=860 max=1477  T=5768s
    measured ratio 0.99283; physical K 0.99285; K~cb 0.99413
  5a_100Hz tl: FAIL n=2 bad=1 z_mean=-1.42 z_std=2.42 max|z|=3.13; amp: 1 evaluated, max rel err 2.29e-03; phase: 1 evaluated, max 0.365 deg; mean=[74.001042-317.164144j] ref=[76.194233-317.412248j] corr=[0.07] FOM Re=[2.18] Im=[27.04]
  5a_100Hz col: FAIL n=2 bad=1 z_mean=-1.42 z_std=2.40 max|z|=3.12; amp: 1 evaluated, max rel err 2.25e-03; phase: 1 evaluated, max 0.364 deg; mean=[74.007889-317.174121j] ref=[76.194233-317.412248j] corr=[0.07] FOM Re=[2.17] Im=[26.86]
  5a_100Hz generations/batch mean=825 max=1614  T=5121s
    measured ratio 0.99240; physical K 0.99241; K~cb 0.99768
Test 5c: conditionally convergent cases with cancellation
  5c_3c_100Hz_eta1 tl: PASS n=2 bad=0 z_mean=+0.78 z_std=1.81 max|z|=2.05; amp: 1 evaluated, max rel err 2.53e-03; phase: 1 evaluated, max 0.264 deg; mean=[77.851501-317.858963j] ref=[76.194233-317.412248j] corr=[0.27] FOM Re=[1.99] Im=[27.29]
  5c_3c_100Hz_eta1 col: PASS n=2 bad=0 z_mean=+0.78 z_std=1.81 max|z|=2.06; amp: 1 evaluated, max rel err 2.53e-03; phase: 1 evaluated, max 0.265 deg; mean=[77.857485-317.857156j] ref=[76.194233-317.412248j] corr=[0.27] FOM Re=[1.99] Im=[27.06]
  5c_3c_100Hz_eta1 generations/batch mean=847 max=2095  T=4674s
  5c_3a_1e4Hz_eta1 tl: PASS n=2 bad=0 z_mean=-0.06 z_std=1.01 max|z|=0.77; amp: 1 evaluated, max rel err 1.18e-03; phase: 1 evaluated, max 0.036 deg; mean=[0.002536-0.074478j] ref=[0.002486-0.074392j] corr=[-0.84] FOM Re=[0.6] Im=[243.87]
  5c_3a_1e4Hz_eta1 col: PASS n=2 bad=0 z_mean=+0.58 z_std=0.65 max|z|=1.04; amp: 1 evaluated, max rel err 2.12e-03; phase: 1 evaluated, max 0.007 deg; mean=[0.00249-0.074234j] ref=[0.002486-0.074392j] corr=[0.4] FOM Re=[4.09] Im=[131.02]
  5c_3a_1e4Hz_eta1 generations/batch mean=560 max=1233  T=1819s

Stage 5 (5a,5c): 5/6 gate checks passed
  FAILED: 5a_100Hz
```
