# N5: the shared TRIGA library (ce_mg_triga CASMO-70, kinetics)

Brief `noise_stage7_brief.md` N5; design doc, Stage 9b preparation. The
library was regenerated in ce_mg_triga (commit 32ad772 checked out on LS6)
as `results/mgxs_ang/casmo70_hist32_kin.h5` (vector fission form, 86.7 MB)
and `casmo70_hist32_kin_matrix.h5` (matrix form, 117.0 MB). Local copies,
outside the repo, are in `~/dev/triga_lib/` (user decision 2026-11-26).

## Content audit (`n5_audit.py`; LS6 and local copies agree)

| Item | Result |
|---|---|
| Groups, delayed families | 70 (CASMO-70, 0 to 20 MeV), 6 |
| Materials | 185 (114 fissile), one temperature (294 K) |
| Scattering | 32-bin μ histogram, isotropic representation |
| Inverse velocity | present in all 185 |
| Delayed data, vector form | `prompt-nu-fission`, `delayed-nu-fission` (6 × 70), `chi-prompt`, `chi-delayed` (6 × 70) |
| Delayed data, matrix form | `prompt-nu-fission` (70 × 70), `delayed-nu-fission` (6 × 70 × 70); χ inside the matrices |
| Decay constants | tallied, 14 slightly different sets of about (0.0133, 0.0327, 0.121, 0.303, 0.850, 2.85) s⁻¹ |
| β | stored as `beta-tallied`, which OpenMC does not read. Not needed: the delayed ν-fission carries β |
| Other datasets | `chi-reference`, `nu-fission-reference`, `kappa-fission`, `fission` (not used by OpenMC's MG transport) |

**Flags.**
- **Empty scatter groups in Air and Aluminum.** 8 Air and Aluminum
  materials have groups with zero scattering but nonzero total (Σt of
  order 1e-4 cm⁻¹). OpenMC's MG collision still chooses scattering with
  probability (Σt − Σa)/Σt. The empty row's outgoing-group loop ends at its
  last entry, so such a "scatter" acts as in-group scattering instead of
  being lost. This is negligible for k and noise at these densities. The
  ce_mg_triga donor repair (`repair_mgxs.py`) covered water, not these
  materials.
- **Multiplicity varying along scatter rows.** Row 0 in 122 materials,
  rows 0–1 in 9: (n,2n) above about 10 MeV. This is OpenMC's MG
  double-counting issue (STATUS open items; ce_mg_triga
  `openmc_issues/ISSUE_mg_multiplicity_row.md`). ce_mg_triga quantified it
  as 0.46 % of secondaries in Zircaloy, 0.14 % in steel and 0.12 % in
  fuel, far below its MC noise. Forward runs only: the N1 transposition
  makes the multiplicity constant along rows.

## Transposed library (N1 + N2) and infinite-medium checks

`adjoint_noise_lib.make_adjoint_noise_library` builds the transposed
library from either form in about 1 s (matrix fission form; 185 materials;
114 fissile). Checks are on the operator OpenMC transports
(`effective_operator`):

| Check | Result |
|---|---|
| Forward balance (scatter + absorption = total) | max error 5.0e-4 (fuel_2951_4772), absorbed by N1's rebalance; 4152 groups with absorption′ set to 0 |
| k∞ transposed against forward, 114 fuel materials | max \|Δk/k\| 2.3e-5 (the same material), both forms |
| k∞ forward, 114 fuel materials | 1.46274 to 1.57710 |
| k∞ matrix against vector form | −1.3e-5 to +1.7e-5 per material, though the matrix χ depends on the incident group (χ spread up to 0.29). Production per g_in and scattering are identical in the two forms |
| Adjoint second-moment radius (N3c) | above 1 in every material (water 1.40, graphite 1.53, fuel up to 1.90; max multiplicity 3.08): adjoint TRIGA runs need the N3c weight windows |

## Still to do

- **Forward MG k against CE.** Take it from ce_mg_triga's own check of
  this library (user decision 2026-11-26), with the status of its
  vector-against-matrix k gap (commit 32ad772): `results/mgxs/kinetics_lib.md`
  and any notes on the gap. The infinite-medium comparison above shows the
  two forms agree per material, so a full-core gap would come from χ(E_in)
  against leakage. The form used in 9b follows ce_mg_triga's resolution.
- **Transposed k in the full core:** an MG eigenvalue run of the TRIGA
  model on the transposed library, against the forward one.
- **Point-kinetic limit:** the TRIGA noise response at low ω against
  G₀(ω) with the IFP β_eff = 715 ± 7 pcm and Λ_eff = 47.36 ± 0.10 µs
  (`results/kinetics/kinetics.json`), on LS6.
