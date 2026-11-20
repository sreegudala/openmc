# Frequency-domain neutron noise in OpenMC

This branch (`mg-complex-weight-noise-solver`) adds a **frequency-domain
neutron noise solver** to OpenMC's multigroup fixed-source mode. It solves
the first-order noise equation for the complex flux response
$\delta\phi(\mathbf r, E, \omega)$ to a small periodic perturbation, using
Monte Carlo particles that carry **complex weights**.

The plan, theory, conventions and verification gates are in
[`../OpenMC_neutron_noise_transport_design.md`](../OpenMC_neutron_noise_transport_design.md)
(the "design doc"). This directory, `noise_dev/`, holds the verification
work: reference solvers, stage drivers, logs and gate records. Progress is
tracked in [`STATUS.md`](STATUS.md).

## Method in brief

- **Complex weights.** Every particle carries `wgt` (Re) and `wgt_im` (Im).
  Random walks depend only on $|w|$. Roulette and splitting use the modulus.
- **The $i\omega/v$ term** uses the collision-based scheme of Rouchon et al.
  The total cross section becomes $\tilde\Sigma_t=\Sigma_t+\eta\omega/v$. At
  each collision a "copy" reaction (probability $\Sigma_\omega/\tilde\Sigma_t$,
  weight $w(1-i/\eta)$) is banked. Capture is implicit.
- **Fission** is divided by the static $k_\text{eff}$ (`noise.keff`).
  - Delayed sites are multiplied by $\lambda_k/(\lambda_k+i\omega)$; no
    emission times are sampled.
  - Each batch runs fission generations until the bank is empty.
  - A growth guard aborts divergent games.
- **Weight cancellation** (optional) averages fission-bank weights on a
  mesh at each generation. This makes the game convergent above the
  frequency cut-off where it otherwise diverges.
- **Noise sources.** There are two kinds:
  - An external source with a complex phase.
  - Sampled from cross-section perturbations $\delta\Sigma_x/\Sigma_x$ in
    chosen cells. Each batch runs `nskip` static generations, and the last
    of these banks complex source sites at its collisions in perturbed
    cells.
- **Adjoint.** The adjoint (bilinear convention: transposed operator at the
  same +ω) is the same solver run on a transposed library built by
  `common/adjoint_noise_lib.py`. The transport code does not change.
- **Tallies.** `ComplexPartFilter` splits each score into Re and Im bins. A
  cross-product accumulator gives the Re–Im covariance, so amplitude and
  phase uncertainties are available through `openmc.noise.complex_tally_stats`.
  `NoisePhaseFilter` separates the static flux from the noise response in
  the same run.

## Minimal input

```python
settings.run_mode = "fixed source"
settings.noise = {
    "frequency": 1.0,            # Hz
    "eta": 1.0,                  # Sigma_omega = eta * omega / v
    "keff": 0.99916,             # static k from an eigenvalue run
    # optional: weight cancellation on a mesh
    "cancellation_mesh": mesh,
    # optional: sampled source instead of the external source
    "perturbations": [{"cells": [cell], "total": -0.041j * np.pi,
                       "scatter": -0.034j * np.pi,
                       "fission": -0.021j * np.pi}],
    "nskip": 20, "static_inactive": 50,
}
tally.filters = [openmc.MeshFilter(m), openmc.NoisePhaseFilter(),
                 openmc.ComplexPartFilter()]
```

Then, after the run:
```python
stats = openmc.noise.complex_tally_stats(sp.get_tally(name="flux"))
# stats["value"], ["amplitude"], ["amplitude_std"], ["phase"], ["phase_std"] ...
```

The full key list is in the `openmc.Settings` docstring (`noise`). The
library must provide inverse velocities, and delayed data for delayed
neutrons. Noise mode is multigroup only.

## Code map (OpenMC changes)

| Area | Files |
|---|---|
| Noise driver: generation loop, guard, cancellation, sampled source, statepoint output | `src/noise.cpp`, `include/openmc/noise.h` |
| Complex weight in particle and bank | `include/openmc/particle_data.h`, `particle.h`, `src/particle.cpp`, `src/source.cpp`, `include/openmc/source.h` |
| Copy reaction, implicit capture, noise fission | `src/physics_mg.cpp`, `src/physics.cpp`, `src/physics_common.cpp`, `src/mgxs.cpp`, `include/openmc/mgxs.h` |
| Settings (`<noise>` element) | `src/settings.cpp`, `include/openmc/settings.h`, `openmc/settings.py` |
| Complex tallies and covariance | `src/tallies/filter_complex_part.cpp`, `src/tallies/filter_noise_phase.cpp`, `src/tallies/tally.cpp`, `tally_scoring.cpp`, `src/state_point.cpp`, `openmc/filter.py`, `openmc/tallies.py` |
| Post-processing | `openmc/noise.py` |
| Python bindings | `openmc/lib/filter.py`, `openmc/lib/core.py` |
| Unit tests | `tests/unit_tests/test_noise.py` |

## This directory

| Path | Contents |
|---|---|
| `reference/` | Independent NumPy references. `ref_g0.py`: point-kinetic $G_0(\omega)$ and cut-offs. `ref_inf.py`: multigroup infinite medium, forward and adjoint. `ref_sn1d.py`: 1-D slab $S_N$ static and noise solver. `ref_diff1d.py`: slab diffusion Green's function. `noise_common.py`: shared data. |
| `common/harness.py` | Runs 20 independent replicas (concurrently on TACC). Implements the pass criteria: STD, STD-C (Hotelling $T^2$ on principal components for correlated mesh tallies), the paired on/off test, the two-sample test against another Monte Carlo code, and FOM. |
| `common/build.py` | Builders for multigroup libraries, slab models, sources and meshes. |
| `common/adjoint_lib.py` | Verbatim copy of the ce_mg_triga library transposition (scattering only, fission removed). Do not edit it here; `test_adjoint_lib.py` and `fixtures/adjoint_lib/` keep it in sync. |
| `common/adjoint_noise_lib.py` | Transposed library for the adjoint noise operator: adds transposed prompt and delayed production matrices to the copy. Also gives the effective operator OpenMC transports and the second-moment diagnostic. |
| `briefs/` | Task briefs from the user (e.g. `noise_stage7_brief.md`). |
| `stageN/run_stageN.py` | Stage driver: builds the cases, runs the replicas, scores against the references, writes `stageN_log_*.md`. |
| `stageN/stageN_gate.md` | Gate record: which runs count, disposition of any misses, and the verdict. |
| `stage6/cortex.py`, `run_stage6c.py` | The CORTEX assembly benchmark and the comparison with Abeille. |
| `tacc/` | SLURM job script, the Abeille build script, and TACC instructions (`tacc/README.md`). |

## Verification approach

Each stage adds one piece of physics and must pass its gate before the next
stage starts. Pass/fail always uses the spread over **20 independent
replicas**, never the batch error bars, which are optimistic for these
weight games. A gate miss is either explained by a rerun with fresh seeds
(recorded in the gate file) or fixed in the code. A criterion changes only
by a dated amendment in the design doc.

Running a stage locally:
```bash
source ~/venvs/openmc-dev/bin/activate
cd noise_dev
python stage4/run_stage4.py            # or --only 4a, see --help
```
Large runs go to TACC LS6; see `tacc/README.md`.
