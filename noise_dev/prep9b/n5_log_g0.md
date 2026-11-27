# N5 point-kinetic limit

```
TRIGA_DIR /home/sg44769/TRIGA_proposal_work/ce_mg_triga, TRIGA_LIB_DIR /home/sg44769/dev/triga_lib
N5 point-kinetic limit: total fission rate dF/F0 against G0(omega) drho, drho real, fitted over all frequencies
  IFP beta_eff 715.1 pcm, Lambda_eff 47.36 us; library physical beta 646.3 pcm; tallied lambda_k [0.01334 0.03272 0.12075 0.30291 0.85035 2.85514]

  IFP beta_eff: drho = -32.64 +- 0.07 pcm; chi2 72.4 / 13
       0.1 Hz: dF/F0 5.6700e-02 at  +158.94 deg; G0 drho 5.6534e-02 at  +159.21 deg; amp ratio - 1 +0.0029 (sigma 0.0065), dphase -0.27 deg (sigma 0.29); R 20, T 72803s
       0.3 Hz: dF/F0 4.9542e-02 at  +168.67 deg; G0 drho 4.8873e-02 at  +169.07 deg; amp ratio - 1 +0.0137 (sigma 0.0080), dphase -0.40 deg (sigma 0.18); R 20, T 60810s
         1 Hz: dF/F0 4.6963e-02 at  +173.37 deg; G0 drho 4.6057e-02 at  +173.61 deg; amp ratio - 1 +0.0197 (sigma 0.0071), dphase -0.25 deg (sigma 0.10); R 20, T 56886s
         3 Hz: dF/F0 4.5069e-02 at  +171.31 deg; G0 drho 4.5234e-02 at  +171.49 deg; amp ratio - 1 -0.0037 (sigma 0.0066), dphase -0.18 deg (sigma 0.12); R 20, T 55129s
        10 Hz: dF/F0 4.1578e-02 at  +156.46 deg; G0 drho 4.2042e-02 at  +157.04 deg; amp ratio - 1 -0.0111 (sigma 0.0062), dphase -0.58 deg (sigma 0.32); R 20, T 55395s
        30 Hz: dF/F0 2.8078e-02 at  +127.82 deg; G0 drho 2.8507e-02 at  +128.64 deg; amp ratio - 1 -0.0150 (sigma 0.0055), dphase -0.82 deg (sigma 0.32); R 20, T 56290s
       100 Hz: dF/F0 1.0521e-02 at  +103.30 deg; G0 drho 1.0665e-02 at  +103.51 deg; amp ratio - 1 -0.0135 (sigma 0.0054), dphase -0.21 deg (sigma 0.24); R 20, T 47360s

  library physical beta: drho = -30.22 +- 0.07 pcm; chi2 321.9 / 13
       0.1 Hz: dF/F0 5.6700e-02 at  +158.94 deg; G0 drho 5.7942e-02 at  +159.09 deg; amp ratio - 1 -0.0214 (sigma 0.0064), dphase -0.15 deg (sigma 0.29); R 20, T 72803s
       0.3 Hz: dF/F0 4.9542e-02 at  +168.67 deg; G0 drho 5.0084e-02 at  +168.93 deg; amp ratio - 1 -0.0108 (sigma 0.0078), dphase -0.27 deg (sigma 0.18); R 20, T 60810s
         1 Hz: dF/F0 4.6963e-02 at  +173.37 deg; G0 drho 4.7153e-02 at  +173.32 deg; amp ratio - 1 -0.0040 (sigma 0.0069), dphase +0.04 deg (sigma 0.10); R 20, T 56886s
         3 Hz: dF/F0 4.5069e-02 at  +171.31 deg; G0 drho 4.6230e-02 at  +170.73 deg; amp ratio - 1 -0.0251 (sigma 0.0065), dphase +0.58 deg (sigma 0.12); R 20, T 55129s
        10 Hz: dF/F0 4.1578e-02 at  +156.46 deg; G0 drho 4.2351e-02 at  +154.92 deg; amp ratio - 1 -0.0183 (sigma 0.0062), dphase +1.54 deg (sigma 0.32); R 20, T 55395s
        30 Hz: dF/F0 2.8078e-02 at  +127.82 deg; G0 drho 2.7384e-02 at  +125.86 deg; amp ratio - 1 +0.0253 (sigma 0.0057), dphase +1.97 deg (sigma 0.32); R 20, T 56290s
       100 Hz: dF/F0 1.0521e-02 at  +103.30 deg; G0 drho 9.9216e-03 at  +102.25 deg; amp ratio - 1 +0.0604 (sigma 0.0058), dphase +1.05 deg (sigma 0.24); R 20, T 47360s

  diagnostic fits (IFP beta_eff):
    drho, Lambda scale: p = -0.00032879 +- 8.4e-07, 1.0277 +- 0.0045; chi2 33.1 / 12
      amp ratio - 1 [%]: -0.4 +0.7 +1.3 -1.0 -1.4 -0.6 +0.5; dphase [deg]: -0.26 -0.38 -0.18 +0.01 -0.02 -0.06 +0.14
    drho, rho0: p = -0.00031987 +- 1.2e-06, 0.00017416 +- 2.4e-05; chi2 22.0 / 12
      amp ratio - 1 [%]: -0.5 +0.8 +1.5 -0.7 -1.2 -0.4 +0.5; dphase [deg]: +0.36 -0.11 -0.09 +0.03 -0.07 -0.13 +0.11

  collision-scheme cut-offs (ref_g0.cutoffs, IFP values):
    eta 1: {'kappa': np.float64(0.41421356237309515), 'omega_CC': np.float64(364.5315478257645), 'omega_D': np.float64(50978.893567867395), 'omega_VLF': np.float64(4.3354590870321894e-05)}
    eta 3: {'kappa': np.float64(0.16227766016837952), 'omega_CC': np.float64(930.4664047140975), 'omega_D': np.float64(130123.57393294353), 'omega_VLF': np.float64(1.6985155009618315e-05)}
```
