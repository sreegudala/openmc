# Stage 5d tuning log

```
Stage 5d tuning: slab 4a, 60 replicas x 20 batches, 20000 particles (x multiplier), seeds 300001 + 7919 r

1.0 Hz
  off: T=9508s gens/batch 947; var check replica/batch Re 0.71 Im 1.17
  h0.5: T=9826s gens/batch 926; var check replica/batch Re 0.97 Im 0.89; |sum w| after/before averaging 1.000, sites per occupied bin 71, removed by roulette 5.54e-08
  h1: T=9773s gens/batch 929; var check replica/batch Re 0.98 Im 1.07; |sum w| after/before averaging 1.000, sites per occupied bin 136, removed by roulette 1.91e-08
  h2: T=9663s gens/batch 928; var check replica/batch Re 0.86 Im 0.79; |sum w| after/before averaging 1.000, sites per occupied bin 264, removed by roulette 4.05e-09
  h4: T=9798s gens/batch 922; var check replica/batch Re 0.93 Im 0.94; |sum w| after/before averaging 1.000, sites per occupied bin 519, removed by roulette 1.60e-09
  h1_noroul: T=10102s gens/batch 921; var check replica/batch Re 1.02 Im 1.23; |sum w| after/before averaging 1.000, sites per occupied bin 137, removed by roulette 0.00e+00
  off_x4: T=32877s gens/batch 1167; var check replica/batch Re 0.86 Im 1.20
  h1_x4: T=35620s gens/batch 1130; var check replica/batch Re 0.76 Im 1.07; |sum w| after/before averaging 1.000, sites per occupied bin 437, removed by roulette 4.53e-09
  FOM h0.5/off: Re 0.90 [0.79, 1.02]  Im 2.66 [2.31, 3.00]  time off/on 0.968
  FOM h1/off: Re 0.95 [0.85, 1.08]  Im 2.91 [2.57, 3.25]  time off/on 0.973
  FOM h2/off: Re 0.88 [0.78, 0.99]  Im 2.76 [2.38, 3.15]  time off/on 0.984
  FOM h4/off: Re 0.92 [0.83, 1.03]  Im 2.69 [2.36, 3.07]  time off/on 0.970
  FOM h1_noroul/off: Re 0.92 [0.82, 1.02]  Im 2.84 [2.50, 3.21]  time off/on 0.941
  FOM h1_x4/off_x4: Re 0.98 [0.86, 1.12]  Im 2.92 [2.56, 3.36]  time off/on 0.923

0.1 Hz
  off: T=9910s gens/batch 1031; var check replica/batch Re 0.65 Im 0.75
  h0.5: T=9959s gens/batch 937; var check replica/batch Re 1.44 Im 1.55; |sum w| after/before averaging 0.999, sites per occupied bin 71, removed by roulette 1.88e-07
  h1: T=9768s gens/batch 942; var check replica/batch Re 1.12 Im 1.07; |sum w| after/before averaging 0.999, sites per occupied bin 135, removed by roulette 5.54e-08
  h2: T=10011s gens/batch 941; var check replica/batch Re 0.86 Im 0.99; |sum w| after/before averaging 0.999, sites per occupied bin 264, removed by roulette 1.19e-08
  h4: T=9845s gens/batch 933; var check replica/batch Re 0.91 Im 0.97; |sum w| after/before averaging 0.999, sites per occupied bin 523, removed by roulette 4.83e-09
  h1_noroul: T=10425s gens/batch 941; var check replica/batch Re 1.05 Im 0.96; |sum w| after/before averaging 0.999, sites per occupied bin 135, removed by roulette 0.00e+00
  FOM h0.5/off: Re 1.08 [0.98, 1.20]  Im 1.76 [1.54, 2.04]  time off/on 0.995
  FOM h1/off: Re 1.12 [0.98, 1.28]  Im 1.71 [1.49, 1.96]  time off/on 1.015
  FOM h2/off: Re 0.99 [0.88, 1.10]  Im 1.60 [1.39, 1.84]  time off/on 0.990
  FOM h4/off: Re 1.02 [0.91, 1.16]  Im 1.69 [1.49, 1.94]  time off/on 1.007
  FOM h1_noroul/off: Re 1.03 [0.91, 1.18]  Im 1.61 [1.41, 1.83]  time off/on 0.951
```
