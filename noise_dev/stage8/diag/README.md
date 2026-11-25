# Stage 8c diagnostic: near-pin amplitude, OpenMC against Abeille

## Problem

In the CORTEX vibrating-pin benchmark at n = 1, OpenMC and Abeille give
the same dipole shape around the pin, but Abeille's amplitude is about 3×
larger.
- Thermal flux on the upper half of the pin cell: Abeille 2.29 ± 0.20 against
  OpenMC 0.95 ± 0.27, about 4σ.
- Both use cancellation on a 170 × 170 mesh.
- OpenMC at the 170 and 340 meshes agree with each other.

## Tests

| Test | Result |
|---|---|
| 8a slab, one interface, against S_N (OpenMC, job 3609592) | detector within 0.05 % (σ 0.9 %) |
| 8a slab, Abeille, local (`abeille_inputs/slab_*.yaml`, 600 batches) | 0.973 of S_N, per-bin 0.97–1.03 (batch σ ~13 %) |
| 1-D vibrating plate: a 0.73 cm fuel plate in the water gap, two opposite interfaces, 0.125 cm cancellation (`plate.py`, `plate_openmc.py`, `abeille_inputs/plate_*.yaml`) | dipole swing Im(12.44 cm) − Im(11.69 cm): S_N +0.1262, OpenMC +0.1292 ± 0.0072, Abeille +0.1197 |
| CORTEX half-pin detectors (`compare_halfpin.py`; adjoint jobs 3611332/3611333, OpenMC no-cancellation job 3611334) | see below |

CORTEX n = 1, thermal flux on half of the pin's lattice element (the 8 × 4
map bins above or below the pin centre), normalised by the static thermal
flux at the assembly centre:

| Estimate | upper half | lower half |
|---|---|---|
| adjoint contraction ⟨χ, δS⟩ (OpenMC static moments, Stage 7 adjoint, P0 + P1 + fission; no cancellation) | +0.740j ± 0.003 | −0.632j ± 0.003 |
| OpenMC forward, no cancellation | +0.28j ± 0.56 | −1.05j ± 0.56 |
| OpenMC forward, cancellation 170 | +0.94j ± 0.27 | −0.43j ± 0.27 |
| OpenMC forward, cancellation 340 | +0.45j ± 0.45 | −0.91j ± 0.45 |
| Abeille, cancellation 170 (job 3607704) | +2.29j ± 0.20 | −1.78j ± 0.21 |

- **The contraction does not depend on which forward run supplies its static
  moments:** 0.737–0.744 above, 0.630–0.637 below.
- **What the contraction shares with OpenMC's forward vibration
  sampling:** nothing.
  - Its inputs are the static flux (OpenMC and Abeille agree to a median
    0.06 % on this assembly, 6c), the adjoint (Stage 7), and the
    c_n/Δ contraction.
  - The contraction reproduces S_N to 1e-3, and 8a's forward to 0.16 %.

## Conclusion so far

**OpenMC's near-pin amplitude is right.**
- The forward estimate agrees with the contraction under every
  cancellation setting, including none.
- Abeille's is 3.1× (above) and 2.8× (below) the contraction, 7.6σ and
  5.5σ away.
- Abeille is right in 1-D, with one interface and with two opposite ones,
  so the problem is specific to the 2-D case.

Remaining Abeille-side suspect: its approximate cancellation, which I
switched on as in the published runs; the shipped input has it off. On a
170 × 170 mesh each cell holds far fewer sites than in the 1-D tests, and
mesh-averaging cancellation is known to bias when cells are sparse. Next:
Abeille n = 1 with the shipped setting (`run_stage8c.py abeille
--harmonics 1 --abeille-no-cancellation --tag nocanc`).

n = 2 (pin-cell sum, thermal): contraction −0.311 ± 0.003; OpenMC
−0.357 ± 0.216 (mesh 170) and −0.185 ± 0.392 (mesh 340); Abeille
−0.684 ± 0.179 (2.1σ above the contraction, same direction).
