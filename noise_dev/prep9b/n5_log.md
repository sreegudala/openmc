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

## Forward MG k (ce_mg_triga, user decision 2026-11-26)

From ce_mg_triga `results/mgxs/kinetics_lib.md` (commit e26de6a, "Task I:
vector-matrix k gap settled"), sections 8–9: multigroup MC, mean of 5
seeds, 100k × (100 + 200) per seed.

| | k | against CE |
|---|---|---|
| CE | 1.0008 ± 0.0005 | |
| vector form | 1.00137 ± 0.00009 | +57 ± 51 pcm |
| matrix form | 1.00152 ± 0.00009 | +72 ± 51 pcm |

- **The gap is settled.** Vector − matrix is −16 ± 13 pcm (−1.2σ). The
  earlier −46 pcm came from one seed.
- **Infinite medium.** ce_mg_triga's own k∞ comparison
  (`results/mgxs_ang/gapcheck.json`) matches the one above: vector −
  matrix within 2.5 pcm per material, k∞ 1.46274–1.57710.
- **Hand-off.** Transpose the readable keys only; this is what
  `adjoint_noise_lib` does. Use the tallied λ_k, and these k for
  `noise.keff`.
- **β.** The library's physical β is 646 pcm (core), against the IFP β_eff
  of 715 ± 7 pcm (β_eff/β = 1.107). χ_d per family is in the library, so
  the multigroup response should carry the delayed-spectrum importance.

## Transposed k on the TRIGA model (`run_n5.py kadj`)

- **Without variance control the run is unusable.** A local smoke run
  (5000 × 20 active) gave 0.586 ± 0.082: the transposed library's second
  moment is unbounded (multiplicity up to 3.08).
- **A flat weight window fixes it.** One mesh cell, bounds 0.5–2.0, gives
  1.014 ± 0.012, consistent with the forward k.
- **LS6 runs** (`n5_log_kadj.md`): 20 replicas × 2·10⁵ × (100 + 400) per
  form with the flat window, against the forward k of the same form.

  | form | transposed | forward | transposed − forward |
  |---|---|---|---|
  | vector | 1.00104 ± 0.00015 | 1.00137 ± 0.00009 | −33 ± 17 pcm (−1.9σ) |
  | matrix | 1.00152 ± 0.00017 | 1.00152 ± 0.00009 | 0 ± 19 pcm |

  - Both are consistent within 2σ. The replica spread (57 and 75 pcm)
    matches the within-run σ (67 and 66 pcm), so the window leaves no hidden
    variance.
  - Transposed vector − matrix is −48 ± 23 pcm (−2.1σ), against the forward
    −16 ± 13 pcm. The two transposed libraries differ only in how the
    vector form's χ is folded into the matrix (k∞ per material within
    1.7e-5), so this is taken as a fluctuation. It does not affect 9b, which
    uses one form throughout: the vector form, whose forward k is
    `noise.keff`.

## Point-kinetic limit (`run_n5.py g0`, `g0-analyze`; `n5_log_g0.md`)

- **Setup.** One LS6 job per F7 frequency: 20 replicas × 25 batches × 10⁴,
  nskip 5, cancellation mesh 25 × 25 × 8 over the fuel box, η = 1 up to
  30 Hz and 3 at 100 Hz, `noise.keff` 1.00137 (vector form).
  - Perturbation: absorption in the central thimble water raised by 10 %
    (cell 3789, r = 1.69 cm, full height). Response: δF/F₀ of the total
    fission rate.
  - Fit: one real δρ over all 7 frequencies (14 real numbers), so the shape
    in frequency and the phase are predictions.
  - Cost: 112 node-hours, 39–61 min per replica (longest at 0.1 Hz).
  - The local 1 Hz probe (20 × 5000) gave δρ = −34.1 pcm with an imaginary
    remainder of 0.6 pcm; each batch ran 400–1300 noise fission generations.
- **Result: the library's kinetics are the IFP β_eff, not the physical β.**

  | G₀ with | δρ [pcm] | χ² / 13 | amplitude residual | phase residual |
  |---|---|---|---|---|
  | IFP β_eff,k (715 pcm), Λ_eff 47.36 µs | −32.64 ± 0.07 | 72 | −1.5 to +2.0 % | −0.8 to −0.2° |
  | library physical β_k (646 pcm) | −30.22 ± 0.07 | 322 | −2.5 to +6.0 % | −0.3 to +2.0° |

  Per frequency against IFP G₀ (σ about 0.6 % and 0.1–0.3°):

  | f [Hz] | 0.1 | 0.3 | 1 | 3 | 10 | 30 | 100 |
  |---|---|---|---|---|---|---|---|
  | amplitude − 1 [%] | +0.3 | +1.4 | +2.0 | −0.4 | −1.1 | −1.5 | −1.4 |
  | phase [°] | −0.27 | −0.40 | −0.25 | −0.18 | −0.58 | −0.82 | −0.21 |

  - The physical-β curve misses by 6 % and 2° at 30–100 Hz, where G₀ is
    set by β/Λ. So the multigroup transport, through the library's χ_d per
    family, reproduces the delayed neutrons' extra importance
    (β_eff/β = 1.107) to within a few per cent.
- **The IFP residual is a 3–4 % tilt, not noise (χ² 72 / 13).**
  Diagnostic fits, one parameter added:
  - **β/Λ scale:** Λ × 1.028 ± 0.005, or equivalently β × 0.973 ± 0.004;
    χ² 33 / 12. The multigroup β_eff/Λ would be 2.7 % below the CE IFP
    value. The IFP β_eff itself is known only to ±1 %.
  - **Reactivity offset ρ₀**, δF/F₀ = δρ / (1/G₀ − ρ₀): ρ₀ = +17 ± 2 pcm;
    χ² 22 / 12, the best. That is the signature of `noise.keff` (1.00137 ±
    0.00009) sitting about 17 pcm below this model's own k (1.9σ of the
    forward k). At 1 Hz a 9 pcm k error moves the amplitude by about 1.5 %.
  - Adding both does not help (Λ × 0.996 ± 0.010, χ² 22 / 11). A constant
    complex space-dependent term does not either (χ² 59 / 11).
  - The remaining scatter (+1.5 % at 1 Hz, −1.2 % at 10 Hz; χ² 22 / 12) is
    at the 2σ level. The space-dependent share the design doc expects
    above about 24 Hz (its paragraph on the point-kinetic plateau) is not resolved here.
- **Verdict: N5's point-kinetic check passes.** The brief sets no numeric
  tolerance. The multigroup response follows G₀ with the IFP β_eff and
  Λ_eff to within 2 % and 0.8° from 0.1 to 100 Hz, and rejects the physical
  β.
- **For 9b.** `noise.keff` uncertainty (±9 pcm) is a 1–2 % effect on
  absolute low-frequency amplitudes. 9b's gates (forward–adjoint
  reciprocity, two mesh sizes, Σ|w|) compare runs that use the same
  `noise.keff`, so it largely cancels there. If 9b's absolute responses
  are to be quoted below 3 Hz, a dedicated vector-form k of this model to
  about ±2 pcm should come first.
