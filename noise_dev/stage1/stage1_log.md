# Stage 1 results

| check | result | detail |
|---|---|---|
| 1a tl survival=False Re bit-identical | PASS | 480 non-zero bins; max |diff| 0.0e+00 |
| 1a tl survival=False Im exactly zero | PASS | max |Im| 0.0e+00 |
| 1a col survival=False Re bit-identical | PASS | 480 non-zero bins; max |diff| 0.0e+00 |
| 1a col survival=False Im exactly zero | PASS | max |Im| 0.0e+00 |
| 1a tl survival=True Re bit-identical | PASS | 480 non-zero bins; max |diff| 0.0e+00 |
| 1a tl survival=True Im exactly zero | PASS | max |Im| 0.0e+00 |
| 1a col survival=True Re bit-identical | PASS | 480 non-zero bins; max |diff| 0.0e+00 |
| 1a col survival=True Im exactly zero | PASS | max |Im| 0.0e+00 |
| 1b tl survival=False Re = cos*base | PASS | max rel err 4.4e-16 |
| 1b tl survival=False Im = sin*base | PASS | max rel err 6.0e-15 |
| 1c tl survival=False corr = 1 | PASS | max |corr - 1| 1.8e-11 |
| 1c tl survival=False sigma_phi zero to rounding | PASS | max sigma_phi 4.2e-09 rad (rounding floor 10*sqrt(eps/(n-1)) = 5.0e-08) |
| 1b col survival=False Re = cos*base | PASS | max rel err 6.7e-16 |
| 1b col survival=False Im = sin*base | PASS | max rel err 7.1e-13 |
| 1c col survival=False corr = 1 | PASS | max |corr - 1| 5.3e-12 |
| 1c col survival=False sigma_phi zero to rounding | PASS | max sigma_phi 4.3e-09 rad (rounding floor 10*sqrt(eps/(n-1)) = 5.0e-08) |
| 1b tl survival=True Re = cos*base | PASS | max rel err 5.6e-16 |
| 1b tl survival=True Im = sin*base | PASS | max rel err 9.5e-15 |
| 1c tl survival=True corr = 1 | PASS | max |corr - 1| 8.9e-12 |
| 1c tl survival=True sigma_phi zero to rounding | PASS | max sigma_phi 4.3e-09 rad (rounding floor 10*sqrt(eps/(n-1)) = 5.0e-08) |
| 1b col survival=True Re = cos*base | PASS | max rel err 5.6e-16 |
| 1b col survival=True Im = sin*base | PASS | max rel err 3.6e-13 |
| 1c col survival=True corr = 1 | PASS | max |corr - 1| 6.3e-12 |
| 1c col survival=True sigma_phi zero to rounding | PASS | max sigma_phi 3.9e-09 rad (rounding floor 10*sqrt(eps/(n-1)) = 5.0e-08) |
| 1c two-source tl covariance vs per-batch | PASS | cov err/scale 9.8e-15, cross rel err 1.1e-16; corr range [-0.60, 0.69] |
| 1c two-source col covariance vs per-batch | PASS | cov err/scale 3.2e-14, cross rel err 3.3e-16; corr range [-0.49, 0.69] |
| 1d Python source file round trip with wgt_im | PASS | fields ('r', 'u', 'E', 'time', 'wgt', 'wgt_im', 'delayed_group', 'surf_id', 'particle') |
| 1d real-only file keeps legacy format | PASS | fields ('r', 'u', 'E', 'time', 'wgt', 'delayed_group', 'surf_id', 'particle') |
| 1d C++ reads wgt_im from file (== real file with phase) | PASS | max rel diff 0.0e+00 |
| 1d legacy file in non-noise mode | PASS | max rel diff vs rotated 2.9e-15 |
| 1d C++ written source file carries wgt_im | PASS | 258 sites, max |phase - theta| 3.3e-16, fields ('r', 'u', 'E', 'time', 'wgt', 'wgt_im', 'delayed_group', 'surf_id', 'particle') |
| 1d 1 vs 8 threads (fission on) | PASS | max rel diff 4.3e-15 |
| 1d event vs history | PASS | max rel diff 5.4e-15 |
| 1d 1 vs 8 threads (complex weights) | PASS | max rel diff 4.1e-15 |
| 1d event vs history (complex weights) | PASS | max rel diff 3.5e-15 |
