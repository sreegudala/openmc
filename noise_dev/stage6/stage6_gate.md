# Stage 6 gate record

Gate (design doc §5, Stage 6): 6a and 6b STD against the references; 6c
STD against Abeille on cells above 10 % of the maximum modulus (STD-C form
per the 2026-11-18 amendment, two-sample), and modulus within 2 % of
Abeille. 6c's secondary FEMFFUSION comparison is informational and was not
run (FEMFFUSION is not built).

All runs on TACC LS6, 20 concurrent replicas x 12 threads.

## Runs

| Run | Job | Commit | Settings | Result |
|---|---|---|---|---|
| 6a, 6b, nskip | 3596433 | 9dc11e56b | `run_stage6.py` defaults | 22/22 |
| 6a, 6b, nskip | 3596752 | 51e661299 | same (6a/6b code unchanged) | **22/22** (`stage6_log_6a_6b_nskip.md`) |
| 6c OpenMC | 3596962 | d7fa398fb | 20 x 600 batches x 20000, nskip 20, 200 static inactive, cancellation 68 x 68, seeds base 1 | `stage6c_openmc.npz` (not committed, 18 MB) |
| 6c Abeille static | | 3dc58ed29 | 1e5 particles, 1050 active generations | k = 0.99912 +- 0.00006, consistent with OpenMC 0.99921 +- 0.00020 (`stage6c_abeille_static.npz`) |
| 6c Abeille noise | 3596959 | d7fa398fb | 20 x 1e5 particles, 600 batches planned (100 ignored), 1260 min cap | 219-230 batches per replica (cap reached) (`stage6c_abeille.npz`, not committed, 12 MB) |
| 6c analysis | (local) | this commit | `run_stage6c.py analyze` | **6/6** (`stage6c_log.md`) |

## Results

**6a, 6b.** See `stage6_log_6a_6b_nskip.md`. 6a: amplitudes within 0.66 %
of $-G_0\varepsilon$ (100 Hz, eta 3; 0.3 % or better below), phases within
0.37 deg. 6b: STD-C at every frequency (Hotelling p 0.71-0.97), the sampled
source integral matches $-\delta\Sigma_t\phi_0$, and nskip 3 vs 20 gives the
same answer.

**6c, CORTEX assembly at 1 Hz vs Abeille.** Both codes use the shipped
$k_\text{eff}=0.99916$, so they solve the same equation; Abeille's static
run reproduces it (0.99912 +- 0.00006).

| Group | Two-sample STD-C (n = 38088 cells) | Assembly integral | Perturbed pin |
|---|---|---|---|
| fast | max\|z\| 1.16, integral z Re -0.94 / Im -0.58, Hotelling k=2 p=0.635 | modulus -1.45 %, phase +0.14 deg | -1.43 %, +0.13 deg |
| thermal | max\|z\| 1.37, integral z Re -0.94 / Im -0.59, Hotelling k=2 p=0.634 | -1.46 %, +0.14 deg | -1.11 %, +0.12 deg |

σ of the modulus ratio on the assembly integral: OpenMC 1.16 %, Abeille
2.22 % (including the 0.08 % static normalisation), combined 2.5 %.

- **Shape.** The per-cell modulus ratio OpenMC/Abeille has a 5-95 % range
  of [-1.6 %, -1.3 %] (fast) and [-1.8 %, -1.2 %] (thermal): the spatial
  pattern agrees to about 0.3 %, and the remaining difference is one
  global amplitude. The z histogram (std 0.2) and Hotelling k = 2 show the
  same thing: the replica scatter is one global (point-kinetic) mode.
- **Static shape** (informational): OpenMC and Abeille $\phi_0$ agree to a
  median 0.06 %, worst cell 0.9 %.

*Power of the 2 % check.* The 1.45 % global offset is 0.6 combined σ, so
it is consistent with zero. But with a combined σ of 2.5 % the check
cannot separate a true 2 % bias from none; it excludes offsets above
about 5 % (2σ). The limit is Abeille: its 20 replicas stopped at the
21 h cap after about 125 active batches each, against 500 planned. A
sharper test needs roughly 4x Abeille's batches (about 4 days of 20 x 12
threads) and is not planned. The shape agreement at 0.3 % and the phase
agreement at 0.14 deg do not depend on this.

## Verdict

Stage 6 **passed**: 6a and 6b 22/22, 6c 6/6. Gated amplitude agreement
with Abeille is resolved to about ±5 % on the global amplitude and 0.3 %
on the spatial shape.
