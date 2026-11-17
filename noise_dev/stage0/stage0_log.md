# Stage 0 — Preflight log (2026-11-16)

## Task 1: by-eye checks of arXiv:2207.09822 (Belanger et al.)

Read from the ar5iv HTML LaTeX alt-text (the arXiv PDF could not be text-extracted
by the fetch tool). LaTeX alt-text is the authors' source, so the transcription is
reliable; nothing was checked against the typeset PDF by eye.

| Item | Finding |
|---|---|
| Eq. 10 | `g(ω) = ∫ e^{−iωt} g(t) dt` — confirms the convention used here. |
| Eq. 15 | `Σα(r,E,t) = Σα(r,E)[1 + εα sin(ω0 t) 1_P(r)]` |
| Eq. 17 | `δΣα = −iπ εα Σα 1_P [δ(ω−ω0) + δ(ω+ω0)]` — the paper prints **plus**. The transform of sin(ω0t) under e^{−iωt} is −iπ[δ(ω−ω0) − δ(ω+ω0)], so the printed "+" is a typo in the paper, not an extraction error. Only the +ω0 component is used in a single-frequency run, so this does not affect Stage 6. |
| Eq. 26 | `c_n = (2/n) sin(n arcsin((x−x0)/ε)) e^{−inπ/2}` — the paper prints **arcsin**; Abeille codes arccos. Resolve numerically before Stage 8 (see note below). |
| Eq. 26 resolved numerically | Direct Fourier coefficients of H(x−x0) − H(x−x0−ε sin ω0t) (2×10⁶-point quadrature, u = ±0.5, 0.2, n = 1–3) match **(2/n) sin(n arccos u) e^{−inπ/2}** to 4 digits for every case; the arcsin form is wrong for odd n (e.g. u=0.5, n=1: direct −1.7321i, arccos −1.7321i, arcsin −1.0000i). Use Abeille's arccos form in Stage 8. The paper's parity sentence is also reversed: odd harmonics are imaginary and **even** in u; even harmonics are real and **odd** in u. |
| Eq. 27 | `c_0 = π − 2 arcsin(u)` for x ≥ x0, `−π − 2 arcsin(u)` for x < x0. |
| Parity sentence | "All even harmonics are purely real and have even spatial parity about x0, whereas all odd harmonics are purely imaginary and have odd parity." |
| Vibrating pin | (Summary only, not verbatim) pin "two cells down and two cells left from the centre", displacement along **x**, ε = 0.2 cm, ω0 = 2π rad/s (4π also tested). Abeille's shipped input vibrates along y for a pin at x=+2.52: still a discrepancy to settle before Stage 8. |
| Cancellation meshes | (Summary only) 170×170 coarse and 340×340 fine; 138×138 scoring mesh. No statement on culling after mesh averaging was found in the portion read. |

## Task 2: search for prior OpenMC work

GitHub issue/PR search (public API, `repo:openmc-dev/openmc`) for "noise",
"frequency domain", "complex weight", "negative weight", "alpha eigenvalue":
no issue or PR about frequency-domain noise or complex weights. Related only:

- openmc-dev/openmc#3971 Delta tracking (open PR) — relevant to exact cancellation (Section 6.3).
- openmc-dev/openmc#1596 Introduce the alpha-k iteration (open PR).
- "noise" hits are unrelated (statistical noise, numerical robustness).

GitHub Discussions and openmc.discourse.group were **not searched**.

## Task 3: external references

- **Abeille** cloned at `c6aaa5f` into `~/dev/external/abeille` (outside the repo,
  GPL-3.0, read and run only) and built (Release, no GUI, no LTO).
  A k-eigenvalue variant of `input_files/noise_oscillation.yaml` (100 000
  particles, 300 generations, 50 inactive) gives k_col = 0.99890 ± 0.00011,
  k_trk = 0.99882 ± 0.00013 against the shipped `keff: 0.99916`:
  a difference of 2.3σ (k_col) with batch σ that is known to be optimistic.
  **Gate: marginal, not a clean pass.** Matters only for the Stage 6c
  comparison, which uses Abeille's own noise run in any case.
- **FEMFFUSION** not built (needs deal.II + PETSc + SLEPc). Its 1-D
  group-2 column re-read was not done. Its role (informational SP3/diffusion
  cross-checks in 4b, 6c, 9a) is not gating; recorded as not done.

## Task 4: library check (CASMO-70, TRIGA)

`~/TRIGA_proposal_work/ce_mg_triga/results/mgxs/mgxs_CASMO-70_*.h5`: 185 data
sets, 114 fissionable, `delayed_groups = 0`. **Every fissionable material lacks
`inverse-velocity`, `decay-rate`, and any delayed data (`beta` /
`delayed-nu-fission`, `chi-delayed`).** Representation is isotropic in the
`_iso` file. The library must be regenerated with these quantities before
Stage 9b. Slowest group speed could not be recorded (no inverse velocity).
Kinetics: IFP β_eff = 715 ± 7 pcm, six β_k available, Λ = 47.36 µs;
λ_k were **not** tallied — Keepin U-235 thermal values are used in Stages 3+.

## Task 5: reference scripts (`noise_dev/reference/`)

| Script | Gate | Result |
|---|---|---|
| `ref_g0.py --test` | Reproduce the design doc's G0 table and f_CC table | PASS (all entries within rounding) |
| `ref_inf.py --test` | Ferrero case → 2.486e-3 − 7.439e-2 i | PASS: 2.4864e-3 − 7.4392e-2 i |
| `ref_sn1d.py --test inf` | Reflective S_N = infinite medium to 1e-6 | PASS: ≤ 4e-13 (fwd and adjoint, 1 and 2 groups, 0.1 Hz–10 kHz) |
| `ref_sn1d.py --test diff` | Thick weakly absorbing slab vs diffusion within 2 % away from source/boundary | PASS: 0.011 % (1 Hz), 0.35 % (100 Hz) |
| `ref_sn1d.py --test conv` | S32 vs S64 and mesh halving to 0.1 % | PASS: ≤ 0.07 % (mesh), ≤ 0.002 % (angle) |
| `ref_diff1d.py --test` | (own check) closed form vs finite differences | PASS: ≤ 1e-6 |

CORTEX data consistency: νΣf,1 = 0.01494071 and β = 0.00535 reproduce
FEMFFUSION's values. CORTEX fuel alone has k∞ = 0.6179 (moderation is in the water).
