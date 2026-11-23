# N3c log

```
N3c-b: infinite-medium adjoint with group windows (7a cases)
  n3 1 Hz (eta 1; window lower bounds [0.112 1.042 0.333] fastest first): adjoint with windows vs ref: PASS n=6 bad=0 z_mean=+0.44 z_std=0.04 max|z|=0.48; amp: 3 evaluated, max rel err 2.82e-03; phase: 3 evaluated, max 0.084 deg
    adjoint without windows vs ref (logged): n=6 bad=0 z_mean=+0.23 z_std=0.00 max|z|=0.23; amp: 0 evaluated, max rel err nan; phase: 3 evaluated, max 0.070 deg
    forward control vs ref (logged): n=6 bad=0 z_mean=-0.51 z_std=1.02 max|z|=1.46; amp: 3 evaluated, max rel err 2.57e-03; phase: 3 evaluated, max 0.006 deg
    FOM (source-group response): adj on 11 (sigma 0.0068, T 1994s, gens/batch 830); adj off 7.61 (sigma 0.0102, T 1259s, gens/batch 728); fwd 155 (sigma 0.0019, T 1866s, gens/batch 1186)
    FOM ratios [95 %]: adj on/off 1.44 [0.646, 3.43]; adj off/fwd 0.0492 [0.025, 0.104]; adj on/fwd 0.0711 [0.0344, 0.162]
  n3 30 Hz (eta 1; window lower bounds [0.112 1.042 0.333] fastest first): adjoint with windows vs ref: PASS n=6 bad=0 z_mean=+0.79 z_std=0.27 max|z|=1.05; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg
    adjoint without windows vs ref (logged): n=6 bad=0 z_mean=-1.15 z_std=0.26 max|z|=1.40; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg
    forward control vs ref (logged): n=6 bad=0 z_mean=-0.34 z_std=0.03 max|z|=0.39; amp: 3 evaluated, max rel err 5.92e-04; phase: 3 evaluated, max 0.083 deg
    FOM (source-group response): adj on 1.65 (sigma 0.0099, T 6210s, gens/batch 2100); adj off 0.786 (sigma 0.0190, T 3537s, gens/batch 1728); fwd 9.6 (sigma 0.0043, T 5600s, gens/batch 3122)
    FOM ratios [95 %]: adj on/off 2.1 [1.14, 3.7]; adj off/fwd 0.0818 [0.0477, 0.145]; adj on/fwd 0.171 [0.102, 0.299]
  n3 100 Hz (eta 5; window lower bounds [0.112 1.042 0.333] fastest first): adjoint with windows vs ref: PASS n=6 bad=0 z_mean=-0.04 z_std=0.21 max|z|=0.26; amp: 3 evaluated, max rel err 1.75e-03; phase: 0 evaluated, max nan deg
    adjoint without windows vs ref (logged): n=6 bad=0 z_mean=-1.15 z_std=0.64 max|z|=1.80; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg
    forward control vs ref (logged): n=6 bad=0 z_mean=-0.72 z_std=1.04 max|z|=1.71; amp: 3 evaluated, max rel err 5.49e-03; phase: 3 evaluated, max 0.004 deg
    FOM (source-group response): adj on 2.62 (sigma 0.0090, T 4699s, gens/batch 1508); adj off 3.07 (sigma 0.0104, T 2994s, gens/batch 1376); fwd 10.4 (sigma 0.0049, T 4016s, gens/batch 2005)
    FOM ratios [95 %]: adj on/off 0.852 [0.329, 1.73]; adj off/fwd 0.295 [0.138, 0.751]; adj on/fwd 0.252 [0.153, 0.408]
  cortex 1 Hz (eta 1; window lower bounds [0.002 0.333] fastest first): adjoint with windows vs ref: PASS n=4 bad=0 z_mean=-0.41 z_std=0.39 max|z|=0.75; amp: 0 evaluated, max rel err nan; phase: 2 evaluated, max 0.025 deg
    adjoint without windows vs ref (logged): n=4 bad=0 z_mean=+2.12 z_std=0.03 max|z|=2.15; amp: 0 evaluated, max rel err nan; phase: 2 evaluated, max 0.238 deg
    forward control vs ref (logged): n=4 bad=0 z_mean=-0.50 z_std=0.64 max|z|=1.07; amp: 2 evaluated, max rel err 4.05e-03; phase: 2 evaluated, max 0.003 deg
    FOM (source-group response): adj on 0.0663 (sigma 0.0195, T 39864s, gens/batch 586); adj off 1.61 (sigma 0.0143, T 3061s, gens/batch 592); fwd 3.25 (sigma 0.0037, T 22500s, gens/batch 1299)
    FOM ratios [95 %]: adj on/off 0.0413 [0.0165, 0.0903]; adj off/fwd 0.494 [0.159, 1.3]; adj on/fwd 0.0204 [0.00707, 0.0471]
  cortex 30 Hz (eta 1; window lower bounds [0.002 0.333] fastest first): adjoint with windows vs ref: PASS n=4 bad=0 z_mean=+0.49 z_std=2.08 max|z|=2.28; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg
    adjoint without windows vs ref (logged): n=4 bad=0 z_mean=-1.14 z_std=0.51 max|z|=1.59; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg
    forward control vs ref (logged): n=4 bad=0 z_mean=-0.43 z_std=1.15 max|z|=1.44; amp: 2 evaluated, max rel err 4.69e-03; phase: 2 evaluated, max 0.007 deg
    FOM (source-group response): adj on 0.0563 (sigma 0.0196, T 46355s, gens/batch 636); adj off 0.58 (sigma 0.0229, T 3283s, gens/batch 575); fwd 3.85 (sigma 0.0034, T 23078s, gens/batch 1331)
    FOM ratios [95 %]: adj on/off 0.097 [0.0427, 0.217]; adj off/fwd 0.151 [0.0673, 0.341]; adj on/fwd 0.0146 [0.00664, 0.0318]
  cortex 100 Hz (eta 5; window lower bounds [0.002 0.333] fastest first): adjoint with windows vs ref: PASS n=4 bad=0 z_mean=-0.23 z_std=0.29 max|z|=0.49; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg
    adjoint without windows vs ref (logged): n=4 bad=0 z_mean=-0.38 z_std=1.45 max|z|=1.63; amp: 0 evaluated, max rel err nan; phase: 0 evaluated, max nan deg
    forward control vs ref (logged): n=4 bad=0 z_mean=-0.30 z_std=0.71 max|z|=0.93; amp: 2 evaluated, max rel err 2.15e-03; phase: 2 evaluated, max 0.003 deg
    FOM (source-group response): adj on 0.0598 (sigma 0.0194, T 44359s, gens/batch 595); adj off 1.23 (sigma 0.0160, T 3196s, gens/batch 603); fwd 4.11 (sigma 0.0032, T 23601s, gens/batch 1351)
    FOM ratios [95 %]: adj on/off 0.0487 [0.0212, 0.133]; adj off/fwd 0.299 [0.108, 0.671]; adj on/fwd 0.0145 [0.00498, 0.0416]

N3c (b): 6/6 gate checks passed
```
