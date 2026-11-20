# Brief: noise branch, preparing Stage 7 (adjoint noise) and Stage 9b

**For:** the Claude Code agent on the `mg-complex-weight-noise-solver` branch (`noise_dev/`).
**Date:** 2026-11-19.
**Based on:** STATUS.md of 2026-11-19 (Stages 0–5 passed, 6c waiting for Abeille).

## Purpose

Stage 7 builds the adjoint noise solution from a transposed multigroup library, the same way the TRIGA repository (`ce_mg_triga`) builds its fission-free adjoints in `step3a_adjoint.py`. The goal is that both projects use **one** transposition, not two that drift apart.

For now that function stays a script in the TRIGA repo. It will move into `openmc.mgxs` later, not in this brief. This branch uses a **copied, frozen** version of it, plus an extension for the parts the noise operator needs and step3a does not have: fission and delayed neutrons.

This brief covers:
- N1, copying the function in;
- N2, extending it;
- N3, three checks on the noise branch that the forward stages never exercised;
- N4, changes to the Stage 7 plan;
- N5, a dependency for Stage 9b.

Stage 6c finishes first. N1–N3 can start now, since they don't touch 6c.

Verification rules are unchanged:
- 20 independent replicas;
- gates as in the design doc;
- misses are rerun with fresh seeds before being called bias;
- criteria change only by a dated amendment.

---

## N1. Copy step3a's transposition into this branch

1. The user copies the transposition function(s) from `ce_mg_triga/step3a_adjoint.py` into `noise_dev/common/adjoint_lib.py`. If you are asked to do the copy, copy only the library-transposition code, with no TRIGA-specific paths or I/O.
2. **Header:** start the file with a comment block giving:
   - the source repo and file;
   - the commit hash it was copied from;
   - the date;
   - the line "Copied verbatim. Do not edit here; extend in adjoint_noise_lib.py."
3. **Do not edit the copied code.** Any change belongs in the TRIGA repo first and is then copied over again.
4. **Sync test** (`tests` or `noise_dev/common/test_adjoint_lib.py`):
   - The TRIGA repo provides a small fixture: a random 4-group library with histogram scattering, and its transposed output from step3a.
   - The test transposes the same input with the copy here and asserts the output is identical.
   - When the copy is refreshed, the fixture is refreshed with it.

**What the copied function does** (check against the code, and record anything that differs in the header):
- It transposes the scattering matrix. Each transfer g→g′ becomes g′→g and keeps its 32-bin μ histogram; the scattering cosine is unchanged under time reversal.
- It rebalances absorption: Σ′_a,g = max(0, Σ_t,g − S_in,g), where S_in,g = Σ_g′ Σ_s,g′→g.
- It sets a per-row multiplicity c_g = S_in,g / (Σ_t,g − Σ′_a,g), which is constant across the row (up to 3.08 in the TRIGA library).
- It keeps Σ_t unchanged.
- It writes **no fission data**, because the TRIGA adjoints are fission-free.

## N2. Extend it to the full noise operator (`adjoint_noise_lib.py`)

Write a function that calls the copied scattering transposition and adds the following.

1. **Prompt fission.** The forward operator is F_p[g←g′] = χ_p,g νΣ_f,p,g′, and its transpose is νΣ_f,p,g χ_p,g′.
   - **Matrix route.** If the branch's fission path accepts a nu-fission matrix, write the transposed matrix. This is the general option.
   - **Rank-1 route.** If it accepts only the χ-vector form, the transpose is still rank 1:
     - χ′_p = νΣ_f,p / Σ_g νΣ_f,p,g;
     - νΣ′_f,p = χ_p · Σ_g νΣ_f,p,g;
     - this is per material, and zero where there is no fission.
   - Use whichever the branch supports, and say which in the log.
2. **Delayed families.** Each family k is transposed the same way, from χ_d,k,g β_k νΣ_f,g′ to β_k νΣ_f,g χ_d,k,g′. The factor λ_k/(λ_k + iω) is a scalar per family and carries through unchanged.
3. **Carried unchanged:** inverse velocity (diagonal) and decay rates.
4. **Fission and absorption keys.** Define and document what the transposed library's fission-related keys mean:
   - the noise fission path must read production from the transposed nu-fission (row sums, in the matrix form);
   - absorption follows N1's rebalance.

   Check that nothing in the noise path reads a non-nu `fission` key in a way that would be wrong after transposition.
5. **Checks** (unit level, no replicas):
   - transposing twice returns the original fission and delayed data;
   - the transposed library's infinite-medium k equals the forward k exactly, from the eigenvalue of the transposed matrices.

## N3. Checks on the noise branch

None of the forward stages exercised these, because forward libraries have multiplicity 1, χ-vector fission, and no weight windows. Each check is a code audit plus a small test.

**N3a. Complex weight under multiplicity.**
- **Audit:** every place that scales `wgt` must scale `wgt_im` by the same real factor. That includes:
  - scatter multiplicity;
  - implicit capture with the rebalanced Σ′_a;
  - survival biasing;
  - the copy reaction;
  - split and roulette;
  - cancellation.

  List each place you checked, with file and line.
- **Test:** a 2- or 3-group infinite medium whose transposed library has c_g of about 2–3.
  - Run the noise solver at three frequencies: one on the plateau, one near the break and one above it.
  - Compare against `reference/ref_inf.py`'s **adjoint** solution, using the Stage 2–3 gate criteria.

**N3b. Fission data in the form N2 produces.**
- Confirm the noise fission path, both prompt and delayed, handles the form N2 writes: a matrix, or rank 1 with χ′ and νΣ′_f.
- If it needs a matrix and supports only vectors, either add matrix support or use the rank-1 route. Record the choice.
- **Test:** the transposed-library infinite medium from N3a with fission and delayed data on, against `ref_inf.py`'s adjoint.

**N3c. Weight windows in noise mode.**
- **Why:** the TRIGA adjoints needed hand-built C_i/φ windows. Without them, G was biased 79% low from weight tails in the pool. Stage 9b cannot reach the out-of-core detectors without windows.
- **Audit:** do weight windows apply in the noise generation loop?
  - Splitting and roulette should act on |w| and keep the phase.
  - Split particles carry both parts.
  - Note how they interact with the cancellation mesh, which acts at generation boundaries.
- **Test:** a deep-penetration slab, absorbing or water-like, 30–50 mean free paths. Compare a run with windows to one without:
  - the means must agree in amplitude and phase (two-sample test on replicas);
  - the FOM must rise;
  - report both.
- If windows don't work in noise mode, say what is missing before changing anything. This is a feature addition and needs a design-doc amendment.

## N4. Changes to the Stage 7 plan

**Convention: a dated design-doc amendment.** The detector signal is written as the bilinear form

  δd(ω) = Rᵀ δφ(ω) = Rᵀ L(ω)⁻¹ δS(ω) = χ(ω)ᵀ δS(ω), with L(ω)ᵀ χ = R.

Consequences:
- The adjoint run uses the transposed library at the **same +ω**, with R as a real external source. No conjugation happens in the run.
- If results are reported as ψ† with the sesquilinear product ⟨a, b⟩ = Σ conj(a) b, then ψ† = conj(χ). Conjugation happens only in post-processing.
- State this in the amendment, and correct any earlier text that says the run uses −ω.

**Stage 7 tests.** These replace the current definitions where they differ; keep the existing pass criteria.
- **7a. Infinite medium, full operator.** Adjoint from the N2 library against `ref_inf.py`'s adjoint, with delayed neutrons on, at about five frequencies from 0.1 to 100 Hz. Check amplitude and phase. This is N3b promoted to a gate.
- **7b. Slab reciprocity, external source.**
  - **Setup:** take a Stage 4 slab and frequency set. δS is the known external source and R is a detector region.
  - **Compare:**
    - forward: δd = Rᵀ δφ;
    - adjoint: χᵀ δS, with χ from an adjoint run with source R.

    Use the two-sample test on the replicas.
  - **Frequencies:** include a frequency below the break (~24 Hz for the TRIGA, but use the slab's own), one at it and one above it.
- **7c. Reciprocity with a sampled source.**
  - **Setup:** use Stage 6's δΣ perturbation, with δS = −δL φ₀.
  - **Absorber first** (capture/total only), so the scalar static flux in the perturbed cells is enough for the contraction.
  - **Then a scattering perturbation.** It needs angular moments of φ₀ and χ (the P0 and P1 terms). Copy the contraction from `ce_mg_triga/step5e_pool_pt.py` the same way as N1:
    - verbatim, with a header and the source commit;
    - with χ now complex;
    - with the same sign convention and the factor 3 in the P1 term.
- **7d (optional). S_N cross-check.** If `reference/ref_sn1d.py` can solve the transposed problem, use it as an independent adjoint reference for the 7b slab.

**Normalization.** Use the TRIGA repo's adjoint source convention, S_i = V_i Σ_R,i, so that G(ω) columns are comparable with the static G at ω → 0.

**Low-frequency check.** As ω → 0 without delayed neutrons, and on a subcritical fixed-source problem, the adjoint noise solution should approach a static adjoint run with the same transposed library. Add this as a sanity check in 7b.

## N5. Dependency for Stage 9b: one shared TRIGA library

STATUS already notes that the CASMO-70 TRIGA library lacks inverse velocity, decay rates and delayed data.
- **Regenerate it once, in the TRIGA repo** (`step1d_mgxs.py`), with inverse velocity, delayed nu-fission per family, χ_d per family, β and decay rates, all tallied from CE. The λ_k come from the tally, replacing the Keepin values.
- Both projects then use that library, and its transposed version comes from N1 + N2.
- **Checks when it arrives:**
  - **Static k:** forward and transposed k in multigroup eigenvalue mode against the CE value.
  - **Point-kinetic limit:** the TRIGA noise response at low ω against G₀(ω) with the IFP values β_eff = 715 ± 7 pcm and Λ_eff = 47.36 µs. The library's β is physical β, not β_eff, so expect a small, explainable difference and report it.
- **Weight windows for 9b:** apply the TRIGA repo's existing C_i/φ windows to |w|, once N3c passes.

## Deliverables

- `common/adjoint_lib.py` (the copy) with its sync test, and `common/adjoint_noise_lib.py`.
- `stage7/stage7_prep_log.md`:
  - the N3 audits with file and line;
  - the N3a–c test results;
  - the fission-form choice (matrix or rank 1).
- The design-doc amendment for the N4 convention, and the updated Stage 7 definitions.
- STATUS.md updated: a Stage 7 prep row, and N3c and N5 under "Open items".

## Not in scope

- Moving the transposition into `openmc.mgxs` (later, its own branch and pull request).
- MPI support for the sampled source.
- Stage 8 (vibration).
