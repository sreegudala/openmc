# N3c gate record: weight windows with complex weights

Gate (design doc, "Stage 9b preparation", amendment 2026-11-22):
- **N3c-u:** unit tests.
- **N3c-a:**
  - STD of the windowed detector response against `ref_sn1d`;
  - two-sample test, windows on against off;
  - FOM ratio on/off, 95 % bootstrap lower bound above 1.
- **N3c-b:** STD of the windowed adjoint against `ref_inf`. FOM ratios are logged.

## Runs

| Run | Job | Commit | Settings | Result |
|---|---|---|---|---|
| N3c-u | local | 317578279 | `tests/unit_tests/test_noise.py` | 27/27 noise tests (complex-phase windows against the analytic answer; generators refused) |
| N3c-a | 3606525 | 317578279 | 20 × 20 × 20000 per set; on seeds 1 + 7919 r, off +500009 | **6/6** (`n3c_log_a.md`) |
| N3c-b | 3606526 | 317578279 | 20 × 20; N3 20000 particles, CORTEX adjoint 2e5, CORTEX forward 20000 | **6/6** (`n3c_log_b.md`) |

The non-noise path is bit-identical to the previous build (an MG window case, 1 thread, compared before and after).

## N3c-a: deep penetration

Water slab, 60 cm (15 fast / 73 thermal mfp), fast source at one face, thermal detector at the other; CADIS windows from S_N.

| f | η | detector vs S_N | on vs off (Welch t Re, Im) | σ on / off | FOM on/off [95 %] |
|---|---|---|---|---|---|
| 1 kHz (phase −32°) | 1 | amplitude 0.11 %, phase 0.00°, max\|z\| 0.42 | +1.49, +0.59 | 0.28 % / 6.6 % | **240** [133, 472] |
| 10 kHz (phase −101°, Re < 0) | 5 | amplitude 0.18 %, phase 0.07°, max\|z\| 0.50 | −1.01, −1.14 | 0.53 % / 151 % | **8600** [647, 1.4e5] |

- **Field with windows against S_N (logged, not gated).**
  - At 1 kHz it passes STD-C, with bins within 0.77 % and 0.32°.
  - At 10 kHz the integral and Hotelling tests pass, but 24 bins have |z| ≥ 3. Those bins are the thermal bins at 2–15 cm, next to the source.
- **Why those bins fail.** The windows are built for the far detector, so near the source their thermal lower bounds are 1e6–1e9. Almost every thermal particle there is rouletted, and those bins rest on a few high-weight survivors with unreliable σ. Windows tuned to one response starve the others; a field tally needs windows built for the field (FW-CADIS).
- **The source-face bin.** The fast bin at the source face misses at both frequencies: its MC σ is below 0.05 %, finer than the S_N discretisation error at a vacuum face with a source.

## N3c-b: adjoint efficiency in infinite media

Windows $C/|\psi_g|$, with ψ the importance for ⟨χ, R⟩. FOM is taken on the source-group response, summed over replicas; 95 % bootstrap intervals in brackets.

| Case | adj (no windows) / fwd | adj (windows) / fwd | windows on/off |
|---|---|---|---|
| N3 1 Hz | 0.049 [0.025, 0.10] | 0.071 [0.034, 0.16] | 1.44 [0.65, 3.4] |
| N3 30 Hz | 0.082 [0.048, 0.15] | 0.17 [0.10, 0.30] | 2.1 [1.14, 3.7] |
| N3 100 Hz | 0.30 [0.14, 0.75] | 0.25 [0.15, 0.41] | 0.85 [0.33, 1.7] |
| CORTEX 1 Hz | 0.49 [0.16, 1.3] | 0.020 [0.007, 0.047] | 0.041 [0.017, 0.090] |
| CORTEX 30 Hz | 0.15 [0.067, 0.34] | 0.015 [0.007, 0.032] | 0.097 [0.043, 0.22] |
| CORTEX 100 Hz | 0.30 [0.11, 0.67] | 0.015 [0.005, 0.042] | 0.049 [0.021, 0.13] |

- **The adjoint-versus-forward gap without windows.** This run puts it at 3–20× for N3 (Stage 7 prep: 9–25×) and 2–7× for CORTEX.
- **N3 with windows: the gap narrows at most modestly.** The gain is significant only at 30 Hz (2.1×); at 1 and 100 Hz it is within the bootstrap interval of 1. Windows control the weight spread from the thermal multiplicity of 2.4, but in a homogeneous medium the importance is flat within a group, so there is little for them to act on.
- **CORTEX with windows: 10–25 times worse.** The gap grows to 50–70×.
  - The cost is the rare thermal-to-fast transfer (probability 0.5 %), the only path into the fast group, the only group where adjoint fission happens.
  - The windows split particles after that transfer but cannot make it more likely. They then pay to carry the many low-weight fast particles they create.
  - Removing this cost needs biasing of the transfer itself (a reaction-channel bias), which windows on space and energy cannot provide.
- **All windowed results are unbiased.** STD passes in all 6 cases: max |z| 2.28, N3 amplitudes within 0.28 %.

## Verdict

N3c **passed**: unit tests, N3c-a 6/6, N3c-b 6/6.
- **The implementation is complex-safe.** The results are unbiased against S_N and `ref_inf`, including where Re w < 0 and where Re of the field changes sign.
- **Windows pay off for deep penetration.** The FOM gain is 240× at 1 kHz and about 8600× at 10 kHz, the 9b out-of-core detector problem.
- **Windows do not close the infinite-medium adjoint gap.** N3 gains at most about 2×, and CORTEX loses 10–25×. That gap comes from multiplicity and rare reaction branches, not from attenuation.
- **Lesson for 9b.** Build windows from the response being computed (detector importance for a detector, FW-CADIS for maps), and do not use group-only windows in regions without spatial attenuation.
