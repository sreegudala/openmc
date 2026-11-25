# Stage 8 gate record: vibrating component

Gate (design doc Stage 8, amendments 2026-11-22b and 2026-11-23):
- **8a:** STD against `ref_sn1d`.
- **8b:** n = 1 linear in ε within statistics.
- **8c:**
  - STD against Abeille on cells above 10 % of the maximum modulus;
  - the two cancellation meshes mutually consistent;
  - the Re of the first harmonic free of the seed-dependent bias of
    uncancelled runs.
- **Second estimator** for every gated detector response: the adjoint
  contraction ⟨χ, δS⟩ (P0 + P1 + fission, δS from tallies), paired test
  against the forward estimate.

Sampling rule and input: amendment 2026-11-23 (each side's own kernel and
fission data; c_n = (2/n) sin(n arccos u) e^{−inπ/2}, derived, and
Abeille's form). The runs were sized from a measured forward σ
(probe, jobs 3607700–3607703), as the amendment requires. All on TACC LS6,
20 concurrent replicas × 12 threads.

## Runs

| Run | Job | Commit | Settings | Result |
|---|---|---|---|---|
| probe: slab, 8c OpenMC, 8c adjoint | 3607700, 3607702, 3607703 | b5c51037d | slab 20 × 20 × 60000; 8c 20 × 60 × 20000, cancellation 170; adjoint 20 × 20 × 20000 | sizing (`stage8_log_probe.md`, `stage8c_log_*_probe.md`) |
| 8a, 8b | 3609592 | 99ee1d592 | 20 × 400 × 100000 static particles, cancellation 1 cm; adjoint 20 × 400 × 20000 | **10/10** (`stage8_log_8a_8b.md`) |
| Abeille n = 1, n = 2 | 3607704, 3607705 | b5c51037d | `noise_vibration.yaml` with approximate cancellation 170 × 170, 20 × 600 × 1e5 | n = 2: 19 replicas (one died on Abeille's lost-particle check, a corner reflection) |
| 8c diagnostics | 3609594, 3611332–3611334 | 99ee1d592, 610a09410 | OpenMC 240 batches at cancellation 0/170/340; half-pin adjoints | `diag/README.md` |
| 8c gate | 3612744 | b665bdfb9 | 20 × 960 × 20000, cancellation 170 and 340, n = 1 and 2 | **16/18**; the 2 misses are dispositioned below (`stage8c_log_analyze_gate.md`) |
| Abeille n = 1, no cancellation | 3612742 | b665bdfb9 | the shipped setting | running (see the end) |

## Results

**8a** (slab with a 4 cm water gap, interface at 10 cm, ε = 0.2 cm,
S_N k = 0.976094):

| n | detector dd/φ₀(D) against S_N | field (STD-C) | forward against contraction | contraction σ |
|---|---|---|---|---|
| 1 | within 0.05 %, σ 0.9 % | pass, phase within 0.16° | 0.16 % (paired t 0.04, −0.17) | 0.18 % (within 0.21 % of S_N) |
| 2 | consistent, σ 8.9 % | pass | consistent (t 1.67, −1.23) | 0.26 % (within 0.45 % of S_N) |

**8b:** (dd/φ₀)/ε at ε = 0.05 and 0.1 cm is 1.005 ± 0.014 and
0.990 ± 0.012 of the 0.2 cm value; each is consistent with S_N. S_N itself
is linear to 1e-4. n = 2 is O(ε²), since ∫c₂ du = 0.

**Deterministic checks.**
- The 8-bin P0 + P1 + fission contraction reproduces S_N's exact detector
  response to about 1e-3 (n = 1 and 2, ε = 0.05 and 0.2).
- The unit test reproduces the sampled-source integral of −iπ/2 (n = 1)
  and 0 (n = 2).

**8c, CORTEX pin vibrating along y** (Abeille's orientation; the paper's
x-vibrating pin is the same problem rotated, amendment 2026-11-22b):

| Check | n = 1 | n = 2 |
|---|---|---|
| fast field vs Abeille (two-sample STD-C) | pass at 170 and 340 | pass |
| thermal field vs Abeille | **fail** at 170 and 340 (Hotelling p 0.000) | pass |
| cancellation 170 vs 340 | pass | pass |
| forward vs contraction, upper/lower half of the pin cell (n = 1) or the pin cell (n = 2) | pass, 4/4 (\|t\| ≤ 1.2) | pass, 2/2 (\|t\| ≤ 2.0) |

Detector values (thermal, normalised by the static thermal flux at the
assembly centre):

| Detector | forward 170 | forward 340 | contraction | Abeille |
|---|---|---|---|---|
| n = 1, upper half | 0.585j ± 0.137 | 0.785j ± 0.210 | 0.737j ± 0.003 | 2.288j ± 0.203 |
| n = 1, lower half | −0.793j ± 0.140 | −0.567j ± 0.212 | −0.637j ± 0.002 | −1.780j ± 0.210 |
| n = 2, pin cell | −0.493 ± 0.101 | −0.394 ± 0.133 | −0.310 ± 0.002 | −0.684 ± 0.179 |

*Disposition of the n = 1 thermal misses against Abeille.*
- **The misses.** Near the pin the two codes have the same dipole shape,
  but Abeille's amplitude is about 3× OpenMC's.
- **The contraction is an independent arbiter.** It does not use OpenMC's
  vibration sampling, and it involves no cancellation. Its inputs are the
  static flux (OpenMC and Abeille agree to a median 0.06 %, 6c), the
  Stage 7 adjoint, and the contraction, which reproduces S_N to 1e-3.
- **It agrees with OpenMC.** OpenMC's forward estimate agrees with it at
  cancellation 0, 170 and 340. Abeille is 3.1× and 2.8× it (7.6σ and
  5.5σ).
- **Abeille is right in 1-D.** It matches S_N on the 8a slab (0.973) and
  on a 1-D vibrating plate with two opposite interfaces (dipole swing:
  S_N 0.126, OpenMC 0.129 ± 0.007, Abeille 0.120). Its excess is specific
  to the 2-D case. `diag/README.md` has the full chain.
- **n = 2 leans the same way.** Abeille is 2.1σ above the contraction.
- **Conclusion: the reference, not OpenMC, is off.**
  - A likely mechanism is Abeille's approximate cancellation on a sparse
    170 × 170 mesh. The published runs used cancellation; the shipped
    input has it off.
  - Abeille's run without cancellation (job 3612742) will show whether
    that is the cause, and is recorded below when it finishes.
  - The disposition does not depend on it.

**Seed-dependent bias of the first harmonic.** OpenMC's Re of the
first-harmonic detector responses is consistent with the contraction at
every cancellation setting, including none (|t| ≤ 1.2). No bias is seen.

**Cost.**
- **The forward estimate is the expensive one.** 8c needs 1.5–2.5 h per
  replica per run for 14–28 % at a half-pin detector. Abeille, on the same
  cells, is about 14× less efficient than OpenMC.
- **The contraction is far cheaper.** It gives 0.3 % per detector from a
  25-minute adjoint job plus the static moments. Even counting the adjoint
  run's time, its FOM is about 700× the forward one.
- **This is the route for vibration detector responses in the TRIGA
  application.**

## Verdict

Stage 8 **passed**:
- 8a 6/6;
- 8b 4/4;
- 8c 16/18, with the two thermal n = 1 misses against Abeille
  dispositioned as an error in the reference (7.6σ and 5.5σ against an
  independent adjoint estimator that agrees with OpenMC with and without
  cancellation).

The vibration source is verified against S_N in 1-D, against the adjoint
contraction in 1-D and 2-D, and against Abeille wherever Abeille is
reliable.

## Abeille without cancellation (job 3612742)

Pending.
