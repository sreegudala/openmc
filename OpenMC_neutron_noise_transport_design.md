# Building frequency-domain neutron noise transport in OpenMC

## Executive summary

**Recommended design.** Implement the collision-based complex-weight scheme of Rouchon et al. in OpenMC's multigroup fixed-source mode. Carry the weight as two doubles (`wgt`, `wgt_im`), run each noise batch as a loop over fission generations, cancel weights on a mesh at each generation boundary, and use a norm-based roulette. This is the best-documented option and the only one with an open reference implementation: Abeille, GPL-3.0, read at commit `c6aaa5f` ([repository](https://github.com/HunterBelanger/abeille)).

**Convention.** Every source that could be opened uses the forward transform $\hat g(\omega)=\int g(t)e^{-i\omega t}dt$, so the operator gains $+i\omega/v$ and delayed production is multiplied by $\lambda_k/(\lambda_k+i\omega)$ ([Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)). This report uses that convention throughout.

**Findings that change the suggested stage order.**

- **Cut-off frequency.** The collision-based scheme needs cancellation above $\omega_{CC}\approx\beta/(\kappa\Lambda)$, with $\kappa=\sqrt{1+\eta^2}-\eta$ ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)). With your $\beta_\text{eff}=715$ pcm and $\Lambda_\text{eff}=47.4\ \mu$s this gives **58 Hz at $\eta=1$**, 102 Hz at $\eta=2$ and 148 Hz at $\eta=3$ (derived here; a point-kinetics estimate). The top of your band therefore needs $\eta\ge3$ or cancellation.
- **Vibration after cancellation.** A vibrating interface produces nearly equal positive and negative source weight, and the CEA group calls cancellation "essentially mandatory" for it ([Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)). Vibration sources move to Stage 8, after cancellation (Stage 5).
- **Roulette.** The 2025 CEA analysis attributes earlier convergence trouble to roulette applied separately to real and imaginary parts and recommends a norm-based roulette that preserves phase ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)). It goes in from Stage 1.
- **Re–Im covariance.** Belanger et al. could not put error bars on modulus and phase because they did not score it ([Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)). It is scored from Stage 1.

**What does not exist.** No published adjoint noise Monte Carlo was found, so Stage 7 is new work; the search was limited and this is not proof of absence. No merged OpenMC commit mentions noise, frequency domain or complex weights ([openmc-dev/openmc](https://github.com/openmc-dev/openmc)), but open issues, unmerged pull requests, Discussions and Discourse were not searched.

**Adjoint convention.** With an unconjugated inner product the adjoint is the transposed library at the **same** $+i\omega$. The $i\omega\to-i\omega$ form in your prompt belongs to the conjugated inner product. The two adjoint functions are complex conjugates and give the same $\delta d(\omega)$ (derivation, Section 1.4).

**OpenMC already holds the data.** Inverse velocity, decay rates, prompt and delayed $\nu\Sigma_f$ and $\chi$ are read and exposed in multigroup mode at commit `013cfbb999ed42da13ec10f1948198c5b9db0765` ([xsdata.h](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/include/openmc/xsdata.h#L66-L93)). No library-format change is needed for forward noise.

**Check before coding.** Several items need a look by eye or a search first; they are listed in Section 6.5 and most of them form Stage 0.

**How to read the labels.** "Sourced" statements carry a link to something the research notes record as read. "Derivation" marks algebra done in the notes or here and not found printed in a source. "Unverified" marks anything not checked against a source. Sources that could not be opened are listed in Section 3.5 and in the reference list.

---

## 1. Theory, with explicit conventions

### 1.1 Conventions used in this report

| Item | Choice here | Where sources differ |
|---|---|---|
| Fourier transform | $\hat g(\omega)=\int_{-\infty}^{\infty}g(t)\,e^{-i\omega t}\,dt$, so $\partial_t\to+i\omega$ | Belanger et al. Eq. 10 and the 2025 Chalmers thesis Eq. 2.50 print this form. Salih (1980) states $e^{+i\omega t}$ but writes equations consistent with $e^{-i\omega t}$; the scan is ambiguous. |
| Opposite convention | Every $i\omega\to-i\omega$; every solution is the complex conjugate (derivation) | No source opened uses $-i\omega/v$. |
| Fission scaling | All fission production divided by the static $k_\text{eff}$ | Same in CEA and Chalmers sources. |
| Delayed factor | $\Upsilon(\omega)=1-\beta+\sum_k\dfrac{\lambda_k\beta_k}{\lambda_k+i\omega}$ | Chalmers writes $1-\sum_k\dfrac{i\omega\beta_k}{i\omega+\lambda_k}$; the two are identical (derivation). |
| Inner product | Bilinear, $\langle f,g\rangle=\sum_g\!\int f_g\,g_g\,d\mathbf r\,d\hat{\boldsymbol\Omega}$, no conjugate | See Section 1.4. |
| Spectra | $\text{CPSD}_{mn}\propto\langle\delta d_m^*\,\delta d_n\rangle$, conjugate on the first index | Phase sign flips with the transform convention and with which signal is conjugated. |

The transform conventions of Yamamoto (2013), Rouchon et al. (2017), Rouchon's thesis, the Pázsit and Demazière handbook chapter and the CORE SIM model report were not read. That they also use $+i\omega/v$ is unverified. The flight factor $e^{-i\omega s/v}$ attributed to Yamamoto second-hand is consistent with it ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)).

### 1.2 First-order noise equation

**Sourced continuous-energy form.** Belanger et al. eliminate the precursors as a convolution, split the operator and flux into static part plus fluctuation, drop the second-order term $\delta\mathcal B\,\delta\varphi$ (the "orthodox linearization") and transform (their Eqs. 1, 8–12; [Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)):

$$\mathcal B(\omega)\,\delta\varphi(\mathbf r,\hat{\boldsymbol\Omega},E,\omega)=-\delta\mathcal B(\omega)\,\varphi_c(\mathbf r,\hat{\boldsymbol\Omega},E)$$

$$\mathcal B(\omega)=i\frac{\omega}{v}+\hat{\boldsymbol\Omega}\cdot\nabla+\Sigma_t-\iint\nu_s\Sigma_s f_s-\frac1{k_\text{eff}}\iint\nu_{f,p}\Sigma_f f_{f,p}-\frac1{k_\text{eff}}\sum_j\frac{\lambda_j}{\lambda_j+i\omega}\iint\nu_{f,d}^{j}\Sigma_f f_{f,d}^{j}$$

**Sourced multigroup forms.** The 2025 CEA paper prints the 1-D multigroup operator with a single $\chi_g$ and $\Upsilon$ as defined above (its Eqs. 70–71; [Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)). The Chalmers diffusion form separates $\chi^p_g[1-\tilde\beta]$ from $\sum_i\chi^d_{g,i}\lambda_i\tilde\beta_i/(i\omega+\lambda_i)$ (Eq. 2.57, reconstructed from a garbled extraction; [Hansen & Hansson 2025](https://odr.chalmers.se/server/api/core/bitstreams/c88f0408-122b-40a7-b838-be90e2d0a4ec/content)). TRIPOLI-4 writes the delayed multiplicity as $\nu^j_{d,\omega}=\dfrac{\lambda_j^2-i\lambda_j\omega}{\lambda_j^2+\omega^2}\nu^j_d$ ([Rouchon et al., PHYSOR 2020](https://www.epj-conferences.org/articles/epjconf/pdf/2021/01/epjconf_physor2020_21002.pdf)), and Abeille codes the same expression ([noise_fission_operator.hpp](https://github.com/HunterBelanger/abeille/blob/master/include/simulation/collision_operators/noise_fission_operator.hpp)).

**Target equation for OpenMC (derivation).** No single source prints exactly this form; it combines the three above, with isotropic fission:

$$\Big[\hat{\boldsymbol\Omega}\cdot\nabla+\Sigma_{t,g}+\frac{i\omega}{v_g}\Big]\delta\psi_g=\sum_{g'}\!\int\!\Sigma_{s,g'\to g}(\hat{\boldsymbol\Omega}'\!\cdot\!\hat{\boldsymbol\Omega})\,\delta\psi_{g'}\,d\hat{\boldsymbol\Omega}'+\frac{1}{4\pi k_\text{eff}}\Big[\chi_{p,g}(1-\beta)+\sum_k\chi_{d,k,g}\frac{\lambda_k\beta_k}{\lambda_k+i\omega}\Big]\sum_{g'}\nu\Sigma_{f,g'}\,\delta\phi_{g'}+\delta S_g$$

$$\delta S_g=-\delta\Sigma_{t,g}\,\psi_{0,g}+\sum_{g'}\!\int\!\delta\Sigma_{s,g'\to g}\,\psi_{0,g'}\,d\hat{\boldsymbol\Omega}'+\frac{1}{4\pi k_\text{eff}}\Big[\chi_{p,g}(1-\beta)+\sum_k\chi_{d,k,g}\frac{\lambda_k\beta_k}{\lambda_k+i\omega}\Big]\sum_{g'}\delta(\nu\Sigma_f)_{g'}\,\phi_{0,g'}$$

Where each requested term enters:

- $1/v_g$ appears only as $+i\omega/v_g$ on the diagonal.
- $\chi_p$ carries the prompt yield $(1-\beta)$.
- Each $\chi_{d,k}$ carries $\lambda_k\beta_k/(\lambda_k+i\omega)$.
- The delayed part of the **source** carries the same factor as the delayed part of the operator.
- If $\beta_k$ depends on the incident group, keep $\nu_{d,k,g'}\Sigma_{f,g'}$ inside the $g'$ sum, as OpenMC's `delayed_nu_fission[dg][gin]` already does.

The operator is exactly critical at $\omega=0$. The static flux and the noise operator must therefore use the same $k_\text{eff}$ (derivation).

### 1.3 Noise sources

All three cases are instances of Belanger Eq. 13, $\delta S=-\delta\mathcal B\,\varphi_c$. They differ only in the space and frequency structure of $\delta\Sigma$.

**Absorber of varying strength.** The CEA model is $\Sigma_\alpha(\mathbf r,E,t)=\Sigma_\alpha[1+\varepsilon_\alpha\sin(\omega_0t)\,\mathbb 1_{\mathcal P}(\mathbf r)]$ (Eqs. 15–16; [Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)). For a pure absorber the multigroup source reduces to

$$\delta S_g(\mathbf r,\hat{\boldsymbol\Omega},\omega)=-\delta\Sigma_{a,g}(\mathbf r,\omega)\,\psi_{0,g}(\mathbf r,\hat{\boldsymbol\Omega}).$$

Their Eq. 17 was transcribed as $\delta\Sigma_\alpha=-i\pi\varepsilon_\alpha\Sigma_\alpha\mathbb 1_{\mathcal P}[\delta(\omega-\omega_0)+\delta(\omega+\omega_0)]$. Under the paper's own convention the transform of $\sin\omega_0t$ has a **minus** between the two deltas (derivation). The "+" is an extraction error or a typo; check the PDF. Abeille uses $-i\pi\varepsilon_xe^{\pm i\phi}$ at $n=\pm1$ ([square_oscillation_noise_source.cpp](https://github.com/HunterBelanger/abeille/blob/master/src/square_oscillation_noise_source.cpp)). For OpenMC, take the complex amplitude $\delta\Sigma_{x,g}/\Sigma_{x,g}$ per region as the user input. Results then differ from Abeille's by the constant factor $-i\pi\varepsilon$, which must be applied when comparing.

**Vibrating component.** For a planar interface at $x_0+\varepsilon\sin\omega_0t$ between materials L and R, the exact perturbation is confined to the swept layer $|x-x_0|\le\varepsilon$ and contains every harmonic (Eqs. 21–27; [Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)):

$$\delta[\nu_\alpha\Sigma_\alpha f_\alpha](x,\omega)=\Delta[\nu_\alpha\Sigma_\alpha f_\alpha]\Big\{c_0\,\delta(\omega)+\sum_{n\ge1}c_n(x,x_0)\big[\delta(\omega-n\omega_0)+(-1)^n\delta(\omega+n\omega_0)\big]\Big\},\qquad\Delta[\cdot]=[\cdot]^L-[\cdot]^R$$

The notes disagree on $c_n$. With $u=(x-x_0)/\varepsilon$:

| Origin | $c_n$ for $n\ge1$ | $n=1$ value |
|---|---|---|
| Paper Eq. 26 as transcribed by a fetch tool | $\frac2n\sin(n\arcsin u)\,e^{-in\pi/2}$ | $-2iu$ |
| Abeille code, read directly ([flat_vibration_noise_source.cpp](https://github.com/HunterBelanger/abeille/blob/master/src/flat_vibration_noise_source.cpp)) | $\frac2n\sin(n\arccos u)\,e^{-in\pi/2}$ | $-2i\sqrt{1-u^2}$ |

Direct evaluation of the $n=1$ Fourier coefficient of the Heaviside difference gives $-2i\sqrt{1-u^2}$, which matches the code (derivation here). The arcsin form would integrate to zero across the layer, whereas the notes expect a net layer source as $\varepsilon\to0$. The parity statement transcribed from the paper (even harmonics even about $x_0$, odd harmonics odd) matches neither form for all $n$, so it also needs checking. Treat the code as authoritative and confirm Eq. 26 in the PDF before coding Stage 8. Other sourced properties:

- Each harmonic is an independent fixed-source problem $\mathcal B_n\delta\varphi_n=-\delta\mathcal B_n\varphi_c$, with $\delta\varphi_{-n}=\delta\varphi_n^*$ (Eqs. 35–37).
- The $n=0$ term is the one the orthodox linearization drops.
- In a symmetric system with a central vibrating rod, the $2\omega_0$ component can exceed the $\omega_0$ one (abstract only; [Zoia et al. 2021](https://research.chalmers.se/en/publication/522444)).
- The linear model over-predicts the second harmonic by more than an order of magnitude for $\varepsilon=0.5$ cm against a non-perturbative solution ([Ferrero thesis](https://webthesis.biblio.polito.it/41032/1/tesi.pdf)).

In multigroup, $\Delta[\nu\Sigma f]$ is the difference of the two materials' group matrices, so the fictitious-union-material sampling that continuous energy needs (Belanger Eqs. 30–34) is unnecessary (derivation). The deterministic "$\varepsilon/d$" model keeps only the fundamental, transcribed as $\delta\Sigma_\alpha=-\frac{i\pi}{2}\Delta\Sigma_\alpha\,\delta(\omega-\omega_0)$ on the swept layer ([Vinai et al. 2023](https://research.chalmers.se/publication/533337/file/533337_Fulltext.pdf)). The point-absorber form $\delta\Sigma_a\approx-\gamma\,\epsilon(\omega)\,\delta'(x-x_p)$ is a derivation from Salih's Eq. 3.54 ([Salih 1980](https://behost.lib.iastate.edu/DR/Salih_ISU-1980-Sa33.pdf)). Its usual name and Pázsit's printed formula are unverified.

**Moderator density fluctuations.** The general source perturbs $\Sigma_t$, scattering and fission kernels together, and a density change enters through $\delta N_k(\mathbf r,\omega)$ multiplying microscopic data (Eqs. 13, 29–31; [Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)). For a pure density change of one constituent with unchanged microscopic group constants, this gives $\delta\Sigma_{x,g}=\Sigma^{(k)}_{x,g}\,\delta N_k/N_k$ for every reaction $x$ of that constituent (derivation). Three gaps remain:

- The notes found no verified mapping from moderator density or temperature to multigroup $\delta\Sigma$ that includes spectral shift. That needs cross-section derivative tables from your CASMO-70 generation step.
- The propagating form $\delta\Sigma(z,\omega)=\delta\Sigma(0,\omega)e^{-i\omega z/u}$ for a coolant transit velocity $u$ is unverified.
- Coupled CORE SIM demonstrates inlet-temperature noise at 0.5 Hz converted to cross-section sources ([manual](https://publications.lib.chalmers.se/records/fulltext/204034/204034.pdf)); its conversion formulas were not extracted.

### 1.4 Adjoint equation and reciprocity

Only one source with explicit adjoint noise equations could be read, a 1980 one-group diffusion thesis whose scan loses daggers and asterisks ([Salih 1980](https://behost.lib.iastate.edu/DR/Salih_ISU-1980-Sa33.pdf)). Everything below is therefore a **derivation** by standard linear algebra, unverified against Pázsit and Demazière, Williams, or Bell and Glasstone.

Write the forward problem as $\mathcal L(\omega)\,\delta\psi=\delta S$ and the detector response as $\delta d(\omega)=\langle\Sigma_d,\delta\psi\rangle$.

**Bilinear product (recommended).** The adjoint operator is the plain transpose: streaming reverses, group matrices transpose, and $i\omega$ keeps its sign.

$$\Big[-\hat{\boldsymbol\Omega}\cdot\nabla+\Sigma_{t,g}+\frac{i\omega}{v_g}\Big]\delta\psi^\dagger_g=\sum_{g'}\!\int\!\Sigma_{s,g\to g'}\,\delta\psi^\dagger_{g'}\,d\hat{\boldsymbol\Omega}'+\frac{\nu\Sigma_{f,g}}{4\pi k_\text{eff}}\sum_{g'}\Big[\chi_{p,g'}(1-\beta)+\sum_k\chi_{d,k,g'}\frac{\lambda_k\beta_k}{\lambda_k+i\omega}\Big]\delta\phi^\dagger_{g'}+\Sigma_{d,g}$$

$$\delta d(\omega)=\langle\Sigma_d,\delta\psi\rangle=\langle\delta\psi^\dagger,\delta S\rangle\quad\text{(no conjugate)}$$

**Sesquilinear product.** With $\langle f,g\rangle_H=\sum_g\!\int f_g^*\,g_g$ the adjoint operator is the conjugate transpose, so $i\omega\to-i\omega$ everywhere, including the delayed factor, and $\delta d=\int(\delta\psi^\dagger_H)^*\,\delta S$. Because $\Sigma_d$ is real, $\delta\psi^\dagger_H=(\delta\psi^\dagger)^*$ and both conventions give the same $\delta d$.

Consequences for the design:

- Your prompt's "transposed library plus $i\omega\to-i\omega$" is the sesquilinear adjoint. It is correct only if the reciprocity integral conjugates the adjoint function. The OpenMC code-map note also proposes running at $-\omega$, while the theory and methods notes use the same $+i\omega$. The two are consistent once the inner product is stated.
- The sesquilinear adjoint at $\omega$ equals the bilinear adjoint at $-\omega$. One transposed-data solver serves both.
- The noise transfer map is $G(\mathbf r_d;\mathbf x',\omega)=\delta\psi^\dagger(\mathbf x',\omega)$: the forward Green's function viewed as a function of its source argument. The abstract of Pázsit and Dykin (2015) attributes this link to the dynamic adjoint ([record](https://research.chalmers.se/en/publication/229494)); the paper was not read.
- One-group diffusion is complex-symmetric, which is why Salih's adjoint equation has the same $s=j\omega$ as his forward one.

In the transposed game the delayed factor attaches to the group in which the adjoint particle collides, through that group's mix of $\chi_p$ and $\chi_{d,k}$ (derivation; Stage 7 gives the sampling rule). The Chalmers thesis uses the static adjoint as a weight with no conjugation shown ([Hansen & Hansson 2025](https://odr.chalmers.se/server/api/core/bitstreams/c88f0408-122b-40a7-b838-be90e2d0a4ec/content)).

### 1.5 Zero-power transfer function and the point-kinetic component

Under the $e^{-i\omega t}$ convention (Eqs. 2.68–2.69; [Hansen & Hansson 2025](https://odr.chalmers.se/server/api/core/bitstreams/c88f0408-122b-40a7-b838-be90e2d0a4ec/content)):

$$\frac{\delta P(\omega)}{P_0}=G_0(\omega)\,\delta\rho(\omega),\qquad G_0(\omega)=\frac{1}{i\omega\Big(\Lambda+\sum_k\dfrac{\beta_k}{i\omega+\lambda_k}\Big)}$$

The same source gives a plateau $|G_0|\approx1/\beta$ with phase near zero for $\lambda<\omega<\beta/\Lambda$. Under the opposite convention $G_0\to G_0^*$ and phases become positive (derivation).

For your kinetics, with one family at $\lambda=0.08\ \text{s}^{-1}$ (an illustrative value, unverified; use your six-family data), the notes and a recomputation here give:

| $f$ (Hz) | $\lvert G_0\rvert$ | $\arg G_0$ |
|---|---|---|
| 0.1 | 140.9 | −7.5° |
| 1 | 139.7 | −3.1° |
| 10 | 129.1 | −22.7° |
| 24 | 98.9 | −45.0° |
| 100 | 32.6 | −76.5° |

The break is at $\beta/\Lambda=150.8\ \text{s}^{-1}$, or 24.0 Hz. Your band covers both the plateau and the $1/(\omega\Lambda)$ roll-off. With six families the plateau is not flat near 0.1 Hz, so compare against the multi-family formula.

The noise splits into a point-kinetic part $\phi_0(\mathbf r)\,G_0(\omega)\,\delta\rho(\omega)$ and a space-dependent part from shape fluctuation ([Hansen & Hansson 2025](https://odr.chalmers.se/server/api/core/bitstreams/c88f0408-122b-40a7-b838-be90e2d0a4ec/content)). A numerical estimator of $G_0$ for a space-dependent solver is their Eq. 2.74:

$$G_0^{dyn}(\omega)\approx\frac1{\delta\rho}\frac{\int\sum_g\frac1{v_g}\phi^+_{g,0}\,\delta\phi_g\,d\mathbf r}{\int\sum_g\frac1{v_g}\phi^+_{g,0}\,\phi_{g,0}\,d\mathbf r}$$

In a small, tightly coupled core the point-kinetic term should dominate on the plateau unless $\delta\rho\approx0$ by symmetry, as for a centred vibration. Above 24 Hz the space-dependent share grows as $|G_0|$ falls. This is a qualitative expectation from the notes, not a sourced result; no source read quantifies it for a TRIGA or for out-of-core detectors.

### 1.6 Spectra

**Sourced.** In one-group diffusion (Eqs. 68–69, reconstructed; [Pázsit group SSM report](https://www.stralsakerhetsmyndigheten.se/contentassets/b1a81ffee6da4b9d970b54b302109cfb/0028-research-and-development-program-in-reactor-diagnostics-and-monitoring-with-neutron-noise-methods)):

$$\text{CPSD}_{\delta\phi}(\mathbf r_1,\mathbf r_2,\omega)=\frac1{D_0^2}\iint G^*(\mathbf r_1,\mathbf r',\omega)\,G(\mathbf r_2,\mathbf r'',\omega)\,\text{CPSD}_{\delta\Sigma_a}(\mathbf r',\mathbf r'',\omega)\,\phi_0(\mathbf r')\phi_0(\mathbf r'')\,d\mathbf r'd\mathbf r''$$

The APSD is the case $\mathbf r_1=\mathbf r_2$. The report's point-kinetic limit prints "$G_0^2$", which must mean $|G_0|^2$ for a real result (a reading of the extraction, not confirmed).

**Derivation for transport and multigroup.** With $\delta d_m=\langle\delta\psi^\dagger_m,\delta S\rangle$:

$$\text{CPSD}_{mn}(\omega)=\iint\delta\psi^{\dagger*}_m(\mathbf x')\,\delta\psi^{\dagger}_n(\mathbf x'')\,\text{CPSD}_{\delta S}(\mathbf x',\mathbf x'',\omega)\,d\mathbf x'd\mathbf x''$$

For one coherent source $\delta S=s(\mathbf x)\,\xi(\omega)$, with $\delta d^{(1)}$ the response to unit $\xi$:

$$\text{APSD}_m=|\delta d_m^{(1)}|^2\,\text{APSD}_\xi,\qquad\text{CPSD}_{mn}=\delta d_m^{(1)*}\,\delta d_n^{(1)}\,\text{APSD}_\xi$$

Coherence is 1 and the phase is $\arg\delta d_n^{(1)}-\arg\delta d_m^{(1)}$. One solve per source per frequency therefore gives all spectra. Tally the static response $d_0$ in the same run and report $\delta d/d_0$, which removes the arbitrary flux normalisation.

**Counting noise and zero-power noise (unverified).** No primary source could be opened; Pázsit and Pál, Williams, and Thie were not accessible. The notes record only general knowledge: measured spectra contain three additive, mutually uncorrelated parts.

- Parametric noise, as above, scaling as the square of the mean signal.
- Detection shot noise, white, scaling as the mean count rate. It appears in each APSD but not in the CPSD of two different detectors.
- Branching (zero-power) noise, correlated between detectors, with a Lorentzian-like shape breaking at $\alpha=\beta/\Lambda$ and relative weight falling with power.

If this holds, zero-power noise shares the 24 Hz break with the plateau edge and the two can be confused at low power. Confirm against one of those references before using any formula.

---

## 2. Monte Carlo algorithms for complex weights

### 2.1 Two ways to handle $i\omega/v$

**Flight-based (Yamamoto; second-hand).** Flights are sampled with the real $\Sigma_t$ and the weight is multiplied by $e^{-i\omega s/v}$ after a flight of length $s$ ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062); [CEA slides](https://www.cea.fr/energies/tripoli-4/Documents/8_rouchon_manu.pdf)). The track-length estimator must then integrate the rotating weight:

$$TL_j(\omega)=\frac{w_{j-1}}{i\omega/v_g}\Big[1-e^{-i\omega s_j/v_g}\Big]$$

MGMC had an error in this estimator. Before the fix it returned $\Re=1.70\times10^{-2}$ against a reference of $2.49\times10^{-3}$ at $10^4$ Hz ([Ferrero thesis](https://webthesis.biblio.polito.it/41032/1/tesi.pdf)). Whether Yamamoto rotates continuously along the track or at collision points in his MCNP implementation is unverified.

**Collision-based (Rouchon et al.).** Add $\eta\frac{\omega}{v}\delta\varphi$ to both sides. Flights use the real $\tilde\Sigma_t=\Sigma_t+\Sigma_\omega$ with $\Sigma_\omega=\eta\,\omega/v$. A "copy" reaction, chosen with probability $\Sigma_\omega/\tilde\Sigma_t$, leaves position, direction and energy unchanged and multiplies the weight by $\nu_\omega=(\eta-i)/\eta$ ([Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822); [Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)). The expected copy weight per collision is $(\eta-i)\frac{\omega}{v}w/\tilde\Sigma_t$, which is exactly the added right-hand-side term, so the scheme is unbiased for any $\eta>0$ (derivation). Standard track-length and collision estimators are unchanged.

Abeille's per-collision sequence, read from code ([noise_branching_collision.hpp](https://github.com/HunterBelanger/abeille/blob/master/include/simulation/collision_operators/noise_branching_collision.hpp)):

1. Fission branching (Section 2.2).
2. Copy with probability $\Sigma_\omega/\tilde\Sigma_t$ and weight $w(1-i/\eta)$, pushed to the secondary bank.
3. Implicit capture: $w\leftarrow w\,\Sigma_s/\tilde\Sigma_t$.
4. Roulette.
5. Scatter.

The default is $\eta=1$, and $k_\text{eff}$ is a mandatory input ([noise.cpp](https://github.com/HunterBelanger/abeille/blob/master/src/noise.cpp)). A forced copy at every collision with weight $w\frac{\Sigma_\omega}{\tilde\Sigma_t}\frac{\eta-i}{\eta}$ would be equivalent, but it is unverified as a published variant.

### 2.2 Fission, delayed neutrons and population control

- **Branching.** Each collision produces $n=\lfloor\nu_f\sigma_f/(\sigma_tk_\text{eff})+\xi\rfloor$ fission noise particles, independent of the parent weight; each inherits the parent's complex weight (Eq. 20; [Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)). In Abeille the denominator already includes $\Sigma_\omega$.
- **Delayed.** A fission neutron is delayed with probability $\nu_d/\nu$. Its weight is multiplied by $\lambda/(\lambda+i\omega)$ and it is emitted immediately from the delayed spectrum ([noise_fission_operator.hpp](https://github.com/HunterBelanger/abeille/blob/master/include/simulation/collision_operators/noise_fission_operator.hpp)).
- **Generations.** Fission neutrons go to the next generation. Copies and other secondaries are processed in the current one ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)). Abeille iterates until the fission bank is empty, with no generation cap ([noise.cpp](https://github.com/HunterBelanger/abeille/blob/master/src/noise.cpp)). OpenMC should add a cap.
- **Roulette.** TRIPOLI-4 and Abeille roulette $|\Re w|$ and $|\Im w|$ separately and kill only when both die ([Rouchon et al., PHYSOR 2020](https://www.epj-conferences.org/articles/epjconf/pdf/2021/01/epjconf_physor2020_21002.pdf); [russian_roulette.hpp](https://github.com/HunterBelanger/abeille/blob/master/include/utils/russian_roulette.hpp)). The norm-based rule is: if $|w|\le w_0$, survive with probability $|w|$ and take new weight $w/|w|$. It preserves phase, does not raise the spectral radius, and is recommended at no significant cost ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)).

The OpenMC code-map note proposes a different fission rule: expected number from $|w|\,|\nu_\text{eff}(\omega)|\Sigma_f/\Sigma_t$ with site weight rescaled. That keeps OpenMC's existing weight-dependent site count but is not a published scheme. Start with the published Abeille rule.

### 2.3 Convergence regimes and the cut-off for your TRIGA

Fauvel et al. define spectral radii for the physical complex operator ($\mathfrak c$ within a generation, $\mathfrak K$ between fission generations) and for the "extended" game on absolute weights ($\tilde{\mathfrak c}$, $\tilde{\mathfrak K}$), with $\mathfrak K\le\tilde{\mathfrak K}$ ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)):

| Regime | Condition | Behaviour |
|---|---|---|
| Convergent | $\tilde{\mathfrak K}<1$, $\tilde{\mathfrak c}<1$ | Works without cancellation |
| Conditionally convergent | $\mathfrak K<1\le\tilde{\mathfrak K}$, $\tilde{\mathfrak c}<1$ | Tally mean can be undefined without enough cancellation; fixed-source runs "did not end" |
| Divergent | $\mathfrak K>1$ or $\tilde{\mathfrak c}\ge1$ | End-of-generation cancellation cannot cure it |

For a one-group infinite medium with one precursor family, $\Lambda=1/(v\Sigma_a)$, $\tau=1/(v\Sigma_t)$, $c_\infty=\Sigma_s/\Sigma_t$ and $\tilde\Upsilon=1-\beta+\beta\lambda/|\lambda+i\omega|$ (their Eqs. 58, 61, 64, 68):

$$\tilde{\mathfrak c}^{cb}=\frac{c_\infty+\sqrt{1+\eta^2}\,\omega\tau}{1+\eta\omega\tau},\qquad\tilde{\mathfrak K}^{cb}=\frac{\tilde\Upsilon}{1-\kappa\omega\Lambda},\qquad\tilde{\mathfrak K}^{fb}=\tilde\Upsilon,\qquad\tilde{\mathfrak c}^{fb}=c_\infty$$

$$\omega^{cb}_{CC}\approx\frac{\beta}{\kappa\Lambda},\qquad\omega^{cb}_D=\frac1{\kappa\Lambda},\qquad\omega^{cb}_{VLF}\approx\frac{2\kappa\lambda^2\Lambda}{\beta},\qquad\kappa=\sqrt{1+\eta^2}-\eta$$

**Evaluation for $\beta_\text{eff}=715$ pcm, $\Lambda_\text{eff}=47.4\ \mu$s (derivation).** The paper says the cut-offs can be evaluated from point-kinetics parameters during power iteration; substituting effective parameters for a heterogeneous 70-group core is an estimate, not a guarantee.

| $\eta$ | $\kappa$ | $f_{CC}$ (cancellation needed above) | $f_D$ (divergent above) |
|---|---|---|---|
| 0.5 | 0.618 | 38.8 Hz | 5.4 kHz |
| 1 | 0.414 | **58.0 Hz** | 8.1 kHz |
| 2 | 0.236 | 101.7 Hz | 14.2 kHz |
| 3 | 0.162 | **147.9 Hz** | 20.7 kHz |
| 5 | 0.099 | 242 Hz | 33.9 kHz |

- Covering 100 Hz without relying on cancellation needs $\eta\ge1.96$. Use $\eta=3$ above about 50 Hz for margin.
- At 100 Hz with $\eta=1$, $\tilde{\mathfrak K}^{cb}\approx1.0053$. Taking the paper's flight-based expression $|\Upsilon/(1+i\omega\Lambda)|\approx0.9924$ as the physical radius, the efficiency model $(1-\alpha)\tilde{\mathfrak K}<\mathfrak K$ requires only $\alpha\gtrsim1.3\%$. That is far milder than the paper's $\alpha>0.94$ case at $10^4$ Hz, where $\tilde{\mathfrak K}=7.50$. This number assumes $\lambda=0.08\ \text{s}^{-1}$, which is unverified.
- $\omega_{VLF}\approx3.5\times10^{-5}$ rad/s at $\eta=1$ with the same $\lambda$, far below 0.1 Hz.
- Across 0.1–24 Hz, $\tilde{\mathfrak K}\approx1-0.9\beta\approx0.994$. Each source particle therefore drives of order $1/\beta\approx140$ fission generations. This matches the "hundreds or thousands" of inner generations per batch reported by Belanger et al.
- Raising $\eta$ costs little in your band. At 100 Hz, $\omega/v\approx3\times10^{-3}\ \text{cm}^{-1}$ for a thermal group at $2.2\times10^5$ cm/s, so copies stay below about 1 % of collisions for $\eta\le3$ (derived estimate; check against your slowest CASMO-70 group).

**The notes disagree on scheme emphasis.** The methods note recommends collision-based with $\eta=1$. The cancellation note recommends "flight-based or frequency-checked collision-based", because the flight-based scheme with norm-based roulette is formally convergent at all frequencies. This report picks collision-based with $\eta$ chosen per frequency, because your band sits far below $f_D$ and the standard estimators stay valid. Flight-based remains a later option for kHz work.

**Earlier claims, second-hand.**

- CEA slides report that with Yamamoto's scheme the particle count "explodes" below 0.01 Hz and above 1 kHz ([CEA slides](https://www.cea.fr/energies/tripoli-4/Documents/8_rouchon_manu.pdf)).
- A third-party paper attributes the same band to Rouchon's method ([Yang et al. 2023](https://www.kns.org/files/int_paper/paper/RPHA-2023_2023_2/RPHA-2023-XuranYang.pdf)).
- The two notes attribute the band to different schemes. Neither is verified against the 2013 or 2017 primary papers.
- The 2025 paper does not support the earlier claim that pseudo-analog games extend collision-based convergence. Such games can end with probability 1 yet have an undefined expectation, which only independent replicas reveal ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)).

### 2.4 Noise-source sampling

**General scheme (CEA).** The source $-\delta\mathcal B\varphi_c$ is sampled on the fly at collisions of one extra power-iteration generation inside the perturbed region, because its sign varies in space, energy and angle. A few plain generations then decorrelate before the next noise batch ([Rouchon et al., PHYSOR 2020](https://www.epj-conferences.org/articles/epjconf/pdf/2021/01/epjconf_physor2020_21002.pdf)). Abeille runs `nskip-1` plain generations, one sampling generation, then the noise simulation ([noise.cpp](https://github.com/HunterBelanger/abeille/blob/master/src/noise.cpp)).

**Rules at a static collision of weight $w$ inside the perturbed region** ([noise_maker.cpp](https://github.com/HunterBelanger/abeille/blob/master/src/noise_maker.cpp)):

| Source particle | Number | Weight |
|---|---|---|
| Copy (same $\mathbf r,\hat{\boldsymbol\Omega},g$) | 1 | $-\dfrac{\delta\Sigma_t}{\Sigma_t}\,w$ |
| Fission | $\lfloor\nu\sigma_f/(\sigma_tk_\text{eff})+\xi\rfloor$ | $w\,\dfrac{\delta\Sigma_f}{\Sigma_f}$; delayed ones also $\times\dfrac{\lambda_j}{\lambda_j+i\omega}$ |
| Scatter (sampled outgoing state) | 1 | $w\cdot\text{yield}\cdot\Big(1-\dfrac{\sigma_a}{\sigma_t}\Big)\dfrac{\delta\Sigma_s}{\Sigma_s}$ |

For a pure absorber only the copy is non-zero. TRIPOLI-4 neglects perturbation of the scattering distribution and keeps only the scattering cross-section fluctuation ([Vinai et al. 2023](https://research.chalmers.se/publication/533337/file/533337_Fulltext.pdf)).

**Normalisation.** Abeille divides source weights by the mean modulus before transport and multiplies the batch tallies back ([noise.cpp](https://github.com/HunterBelanger/abeille/blob/master/src/noise.cpp)). Because $|\delta\Sigma/\Sigma|$ is small, this keeps roulette thresholds meaningful.

**Vibration.** The source lives on the swept layer; each harmonic is a separate run with spatial weight $c_n$. Belanger et al. give exact sampling through a fictitious union material; TRIPOLI-4's earlier scheme was approximate (Eq. 28). For a vibrating pin the source is positive on one side and negative on the other, giving what the authors call "catastrophic convergence issues" ([Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)).

**Alternative.** NECP-MCX tallies the static source and converts it to a complex source with a script ([Yang et al. 2023](https://www.kns.org/files/int_paper/paper/RPHA-2023_2023_2/RPHA-2023-XuranYang.pdf)). That route needs the angular flux for $-\delta\Sigma_t\psi_0$ and a spatial mesh; on-the-fly sampling needs neither.

How Yamamoto samples his source, and Rouchon 2017's own source equations, were not read.

### 2.5 Weight cancellation

**Why it is needed.** Signed-weight games carry a second, non-physical spectrum whose dominant eigenvalue is at least the physical one. Without cancellation the sum of $|w|$ follows the non-physical eigenvalue while the net signal follows the physical one ([Belanger, Mancusi, Zoia 2021](https://ar5iv.labs.arxiv.org/html/2103.13891)). In the plateau both schemes are formally convergent, so cancellation there is variance reduction over chains of about $1/\beta$ generations. Above $\omega_{CC}$, and for mixed-sign sources, it is required.

**Where it is applied.** At a rendezvous at the end of each fission generation, to all banked fission particles, never inside a generation. Fission is chosen because emission is isotropic and nearly independent of incident energy, so binning can key on position alone. Copy particles must stay within the generation ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)). A fixed-source noise batch must therefore be restructured into inner generations ([Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)).

| Method | Rule | Bias | Evidence on cost and gain |
|---|---|---|---|
| Mesh averaging, all particles kept (Zhang-style, used in MGMC) | Sort fission particles into a rectilinear mesh; assign each the complex mean weight of its cell | Biased by intra-cell averaging; cannot be bounded a priori without repeating at several mesh sizes ([Belanger et al. 2022](https://ar5iv.labs.arxiv.org/html/2205.12062)) | FOM gain 234–304 on total time and 1252–1624 on noise time at 170×170; 122–129 and 461–487 at 340×340 ([Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)) |
| Mesh resampling (Yamamoto-style; second-hand) | Sum weights per cell, resample particles of equal weight | Same mechanism | Exact rule, bin sizes and Re/Im treatment unverified (2013 paper not read) |
| Exact regional (Booth and Gubernatis; Belanger et al.) | Split each weight into a point part $w(\zeta-\beta)/\zeta$ and a uniform part $w\beta/\zeta$; sum uniform parts per region and re-emit uniformly | Unbiased for any $\beta$ independent of the landing point, provided $\zeta>0$ over the region ([Belanger et al. 2022](https://ar5iv.labs.arxiv.org/html/2205.12062)) | FOM gain about 8 (total) and 14 (noise) on the same mesh; "cancels much less total weight"; "limited production-code potential" ([Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)) |
| Pairwise point estimator (Booth 2003; second-hand) | Cancel between pairs | — | Quadratic in walkers ([Belanger, Mancusi, Zoia 2021](https://ar5iv.labs.arxiv.org/html/2103.13891)) |

**Convergence conditions.**

- Imperfect cancellation with efficiency $\alpha$ gives $\tilde{\mathfrak K}^{(\alpha)}=\max(\mathfrak K,(1-\alpha)\tilde{\mathfrak K})$. Convergence needs $\tilde{\mathfrak c}<1$ and $(1-\alpha)\tilde{\mathfrak K}<\mathfrak K$. The efficiency of real algorithms is "generally unknown" ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)).
- A minimum efficiency and a critical net weight exist below which the non-physical mode wins ([Belanger, Mancusi, Zoia 2021](https://ar5iv.labs.arxiv.org/html/2103.13891)).
- Fission-only cancellation is sufficient whenever $\tilde{\mathfrak c}<1$. That holds always for flight-based and for $\omega<\omega_D$ for collision-based.

**Known failure modes.**

- Too few particles per cell: efficiency drops, as the fine-mesh gain of 487 against 1624 shows.
- Too coarse a mesh: an undetected bias. The coarse mesh in the benchmark carried a "slight but undetected bias" ([Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)).
- Intra-generation divergence above $\omega_D$.
- "Symmetric divergence": the mean is compensated but the variance grows and the run does not end.
- Under-estimated error bars from correlated batches. The no-cancellation baseline differed between seeds.

**Cost.**

- Published implementations gather the whole fission bank on one node at every inner generation; that gather is the main parallel-scaling liability. A 2024 paper addresses distributed algorithms (abstract only; [Grablevsky & Belanger](https://epjwoc.epj.org/articles/epjconf/abs/2024/12/epjconf_snamc2024_09007/epjconf_snamc2024_09007.html)).
- For mesh averaging only the per-cell complex sum and count are global, so an `MPI_Allreduce` over a mesh array avoids the gather. This is an inference, not a published design.
- Exact cancellation stores three extra items per fission particle. With surface tracking it needs an optical-path ray trace per $\zeta$ evaluation.
- Part of the published gain comes from the inner-generation restructuring alone and was not separated.

**Mesh guidance.** The only noise benchmark used cells of 0.06–0.13 cm (10–20 per 1.26 cm pitch) with $10^6$ particles per generation. No guidance exists for hexagonal or TRIGA lattices, for 0.1 Hz, or above the plateau in a heterogeneous system.

### 2.6 Adjoint noise Monte Carlo

No published Monte Carlo solution of the adjoint frequency-domain noise equation with complex weights was found. Coverage was limited, so this is a negative search result and not proof of absence.

- The 2025 CEA paper, the MGMC thesis and the Chalmers MCNP method do not treat an adjoint. Abeille's noise mode is forward only ([noise.cpp](https://github.com/HunterBelanger/abeille/blob/master/src/noise.cpp)).
- Yamamoto and Sakamoto's 2022 paper on sensitivity analysis of neutron noise (doi:10.1016/j.pnucene.2022.104130) could not be opened. Whether it contains a Monte Carlo adjoint is unknown ([KAKENHI report](https://kaken.nii.ac.jp/file/KAKENHI-PROJECT-17K07011/17K07011seika.pdf)).
- The nearest published idea is the Chalmers "modified Green's function" method: forward, region-discretised and real-valued, on unmodified MCNP6.2. It was tested at 0.0135 Hz and 0.159 Hz against CORE SIM at about 10 % in amplitude and 2° in phase ([Demazière et al., PHYSOR 2020](https://research.chalmers.se/publication/524729/file/524729_Fulltext.pdf)).

That a transposed-library adjoint walk can reuse the collision-based $\eta$ trick unchanged is an inference; no source does it.

### 2.7 Tallies and statistics

- All codes score real and imaginary parts as two real tallies. Abeille exposes `real-flux`, `imag-flux`, `real-source` and `imag-source`, with one sample per noise batch ([noise_vibration.yaml](https://github.com/HunterBelanger/abeille/blob/master/input_files/noise_vibration.yaml)). TRIPOLI-4 uses independent replicas.
- Detector scoring is a linear functional $R=\int h\,\psi$ of the complex field ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)).
- For $A=\sqrt{x^2+y^2}$ and $\varphi=\operatorname{atan2}(y,x)$, standard first-order propagation (not from a noise paper) gives

$$\sigma_A^2=\frac{x^2\sigma_x^2+y^2\sigma_y^2+2xy\,\text{cov}}{A^2},\qquad\sigma_\varphi^2=\frac{y^2\sigma_x^2+x^2\sigma_y^2-2xy\,\text{cov}}{A^4}$$

- Batch error bars can be under-estimated because noise batches share a power-iteration lineage. The MGMC thesis used 20 skipped generations against 3 in Belanger et al., and the 2025 paper used 20 independent replicas.
- No source read reports FOM against frequency. Rouchon 2017 is said to compare the two schemes; that is unverified.

### 2.8 Comparison and recommendation

| Method | $i\omega/v$ | Fission and delayed | Source sampling | Cancellation | Adjoint | Codes | Demonstrated range |
|---|---|---|---|---|---|---|---|
| Yamamoto flight-based (2013, 2018); second-hand | Real $\Sigma_t$ flights; $w\times e^{-i\omega s/v}$ | Complex delayed multiplicity; two-part roulette; details unverified | Unverified | Binning; needed outside the plateau | None found | In-house; MCNP-based | About 0.01–1000 Hz without cancellation (slides) |
| Rouchon collision-based (2017) | Real $\tilde\Sigma_t=\Sigma_t+\eta\omega/v$; copy $\times(\eta-i)/\eta$ | $n=\lfloor\nu\Sigma_f/(\tilde\Sigma_tk)+\xi\rfloor$; delayed $\times\lambda/(\lambda+i\omega)$ | On the fly at collisions of an extra power-iteration generation | Not in the plateau for oscillations; yes above $\omega_{CC}$ and for vibrations | None found | TRIPOLI-4, MGMC, Abeille | 0.1 Hz–$10^4$ Hz multigroup tests; 1–3 Hz assembly |
| Belanger et al. vibration sampling (2023) | As Rouchon, $\eta=1$ | Same | Exact union-material sampling per harmonic | "Essentially mandatory" | None | MGMC, Abeille | 1 and 2 Hz, harmonics 1–2 |
| NECP-MCX k–s iteration (2023) | Rouchon pseudo cross section | Eigenvalue-style cycles | Tallied, converted by script | Claimed unnecessary | None | NECP-MCX | 1 Hz only; no errors or timing |
| Chalmers modified Green's functions (2020) | No complex weights; Re/Im coupled outside the code | Standard subcritical fixed source | Dirac-like per region and group | No | Forward only | MCNP6.2 unmodified | 0.0135 and 0.159 Hz |

**Recommendation for a first multigroup implementation.**

1. Collision-based scheme, $\eta$ user-settable; default 1, with a warning when $\omega>\beta/(\kappa\Lambda)$ if kinetics parameters are supplied.
2. Abeille's fission branching with $1/k_\text{eff}$ and the delayed factor.
3. Fission-generation loop with a generation cap and a $\sum|w|$ monitor.
4. Norm-based roulette on $|w|$.
5. Mesh-averaged cancellation at the rendezvous, all particles kept, followed by a norm-based roulette of the bank. This last step is a design choice; the notes could not establish how MGMC culls after averaging.
6. On-the-fly source sampling with mean-modulus normalisation.
7. Re and Im tallies with the cross term accumulated.
8. Exact regional cancellation later, as a bias-free cross-check only.

The recommendation rests on CEA-authored sources and the Abeille code. Yamamoto's side was seen only second-hand, so any statement that his scheme is less stable is the CEA group's characterisation.

---

## 3. Existing implementations and references for verification

### 3.1 Codes

| Code | Type | Availability and licence | Use as reference |
|---|---|---|---|
| **Abeille** ([repo](https://github.com/HunterBelanger/abeille), commit `c6aaa5f`, 2026-02-16) | Monte Carlo, CE and multigroup, complex weights, oscillation and flat-vibration sources, approximate and exact cancellation | Public, GPL-3.0; C++20, OpenMP, optional MPI | Tally-for-tally Monte Carlo reference. Ships the 2-group assembly noise inputs. Read and run it, but do not copy code into MIT-licensed OpenMC. |
| MGMC ([repo](https://github.com/HunterBelanger/mgmc), commit `1231c27`, 2022-11-03) | Monte Carlo, multigroup predecessor | Public, CeCILL-2.1 | Same noise inputs; used in the CEA papers |
| TRIPOLI-4 | Monte Carlo, CE, original production noise solver | NEA Data Bank, single-user R&D licence, no source ([CEA page](https://www.cea.fr/energies/tripoli-4/distribution)) | Published numbers only |
| Yamamoto's MCNP-based code | Monte Carlo, flight-based | No public release found (unverified) | — |
| **FEMFFUSION** ([repo](https://github.com/Zonni/FEMFFUSION), commit `1611f90`, 2025-07-28) | Deterministic FEM diffusion and SPN, frequency-domain noise, any group count, forward or adjoint | Public, GPL-3.0; deal.II, PETSc, SLEPc | Runnable diffusion and SP3 reference for the slab, CORTEX assembly and C5G7 |
| CORE SIM / CORE SIM+ | Two-group diffusion, MATLAB, noise and adjoint | GPL, but distributed on request from Chalmers under an agreement ([Demazière 2022](https://research.chalmers.se/publication/530864/file/530864_Fulltext.pdf)) | Diffusion reference if obtained |
| NOISE-SN | Frequency-domain $S_N$, any groups and families | Availability and licence not stated ([Yi et al.](https://research.chalmers.se/publication/525456/file/525456_Fulltext.pdf)) | Published numbers only |
| APOLLO3 (IDT) | $S_N$ short characteristics | Availability not found (unverified; presumed CEA-restricted) | Published numbers only |

- SCONE has no noise mode: a grep of commit `a4a02bf` returned no hits ([repo](https://github.com/CambridgeNuclear/SCONE)).
- "No noise mode" for Serpent, McCARD, RMC and MC/DC is unverified; it rests on absence of evidence.
- Chenille was not found as a public repository.
- SIMULATE-3K, PARCS, DYN3D, PROTEUS and Griffin were not researched.
- The cancellation note could not access GitHub and gives Abeille's canceller file paths as an unverified recollection. The other two notes cloned the repository and confirm an `exact_mg_cancelator` class and a `cancelator: {type: approximate, shape: [170,170,170]}` input block; individual canceller source files were not read.

### 3.2 OpenMC prior work

Scope of the search: the full commit log to `013cfbb` (2026-11-11), which covers merged pull requests only. Open issues, unmerged or closed pull requests, GitHub Discussions and openmc.discourse.group were **not searched**. The source tree itself was not grepped for "noise".

- No commit message matches "noise", "frequency", "negative weight" or "weight cancel". "complex" matches only geometry and `std::complex` for multipole data.
- IFP kinetics: [#3133](https://github.com/openmc-dev/openmc/pull/3133), multigroup IFP [#3425](https://github.com/openmc-dev/openmc/pull/3425), fixes [#3580](https://github.com/openmc-dev/openmc/pull/3580), and IFP restricted to eigenvalue mode [#3681](https://github.com/openmc-dev/openmc/pull/3681).
- Delayed-neutron emission time and combing: [#2992](https://github.com/openmc-dev/openmc/pull/2992). Multigroup delayed $\chi$ fix: [#4067](https://github.com/openmc-dev/openmc/pull/4067).
- Random-ray adjoint: [#3191](https://github.com/openmc-dev/openmc/pull/3191); local adjoint sources: [#3717](https://github.com/openmc-dev/openmc/pull/3717).
- Shared secondary bank: PR #3863, per the 0.16.0 release notes.
- Outside the repository: an event-tracking feature for zero-power (Rossi-alpha) noise, a different problem (abstract only; [Saliba et al. 2025](https://www.epj-conferences.org/10.1051/epjconf/202533804001)); and a static alpha-eigenvalue capability "in a modified version" of OpenMC ([arXiv:2606.00907](https://arxiv.org/abs/2606.00907)), which is point-kinetics post-processing and has no complex-weight machinery.

### 3.3 Analytic solutions

**Zero-power transfer function.** Section 1.5, sourced.

**One-group infinite medium (derivation, checked numerically).** With $k=\nu\Sigma_f/\Sigma_a$ and $\Lambda=1/(v\Sigma_a)$, the operator reduces to $\Sigma_a/G_0(\omega)$, so for a uniform source

$$\delta\phi(\omega)=\frac{G_0(\omega)}{\Sigma_a}\,\delta S(\omega),\qquad\frac{\delta\phi}{\phi_0}=-G_0\frac{\delta\Sigma_a}{\Sigma_a}\ \text{(absorber)},\qquad\frac{\delta\phi}{\phi_0}=G_0\,\Upsilon\,\frac{\delta\Sigma_f}{\Sigma_f}\ \text{(fission)}$$

Check: for the MGMC case ($\Sigma_t=1$, $\Sigma_s=0$, $\Sigma_f=0.4\ \text{cm}^{-1}$, $\beta=7\times10^{-3}$, $\lambda=0.08\ \text{s}^{-1}$, $v=3\times10^5$ cm/s, $f=10^4$ Hz, $\varepsilon_f=5\times10^{-3}$) the reference is $\Re=2.49\times10^{-3}$, $\Im=-7.44\times10^{-2}$ ([Ferrero thesis](https://webthesis.biblio.polito.it/41032/1/tesi.pdf)). The fission formula with amplitude $\pi\varepsilon_f$ gives $2.486\times10^{-3}-7.439\times10^{-2}i$. The factor $\pi$, with no $-i$, implies that case uses a cosine-type perturbation normalised as $\pi\varepsilon$; that is an inference from the numbers. Scattering cancels out of the infinite-medium result.

**Multigroup infinite medium (derivation).** Solve the $G\times G$ complex system $\big[\mathrm{diag}(\Sigma_{t,g}+i\omega/v_g)-S^T-\tfrac1k\boldsymbol\chi(\omega)(\nu\Sigma_f)^T\big]\delta\boldsymbol\phi=\delta\mathbf S$, with $\boldsymbol\chi(\omega)=\chi_p(1-\beta)+\sum_k\chi_{d,k}\beta_k\lambda_k/(\lambda_k+i\omega)$ and $k$ the dominant eigenvalue of the same matrices at $\omega=0$.

**One-group slab, diffusion (sourced equation, derived solution).** The SSM report writes $\nabla^2\delta\phi+B^2(\omega)\delta\phi=\delta\Sigma_a\phi_0/D_0$ with $B^2(\omega)=B_0^2[1-1/(\rho_\infty G_0(\omega))]$ and $\nabla^2G+B^2G=\delta(\mathbf r-\mathbf r')$ (Eqs. 51–55, reconstructed; [SSM report](https://www.stralsakerhetsmyndigheten.se/contentassets/b1a81ffee6da4b9d970b54b302109cfb/0028-research-and-development-program-in-reactor-diagnostics-and-monitoring-with-neutron-noise-methods)). For a bare slab with $G(\pm a)=0$, the notes derive

$$G(x,x',\omega)=-\frac{\sin B(a+x_<)\,\sin B(a-x_>)}{B\,\sin2Ba},\qquad\delta\phi(x,\omega)=\frac1{D_0}\int G(x,x',\omega)\,\delta\Sigma_a(x',\omega)\,\phi_0(x')\,dx'$$

and for an infinite medium $G(r,\omega)=-e^{-iB r}/(4\pi r)$ with $\operatorname{Im}B<0$. Sources that put $-\delta$ on the right or absorb $1/D_0$ into $G$ differ by a sign and a factor $D_0$. Pázsit and Analytis (1980) was not opened, so sign and normalisation relative to it are unverified.

**Two-group slab.** Explicit analytic Green's functions could not be retrieved. Jonsson et al. (2012, doi:10.3139/124.110258) is the pointer; only its abstract was read. Do not quote two-group formulas without the paper. Use the $S_N$ script of Stage 0 as the two-group reference.

### 3.4 Published benchmarks

| Benchmark | Geometry and data | Availability | Published agreement |
|---|---|---|---|
| **CORTEX 2-D UOX assembly, cross-section oscillation, 1 Hz** | Reflective 17×17, pitch 1.26 cm, 264 square pins 0.7314 cm, 25 water holes, 0.08 cm water blade, 2 groups without up-scatter, 1 precursor family; $\varepsilon$ = 4.1 % total, 3.4 % scatter, 2.1 % fission in one pin | Complete in Abeille [`noise_oscillation.yaml`](https://github.com/HunterBelanger/abeille/blob/master/input_files/noise_oscillation.yaml); cross-check FEMFFUSION [`2D_UOX_FA`](https://github.com/Zonni/FEMFFUSION/tree/master/examples/2D_UOX_FA) | TRIPOLI-4 against APOLLO3 modulus within 2 %, $\sigma<4$ % ([PHYSOR 2020](https://www.epj-conferences.org/articles/epjconf/pdf/2021/01/epjconf_physor2020_21002.pdf)) |
| Same assembly, isotropic source $S=-1+i$, 3 Hz, thermal group, one pin | As above | Pin location only in the paper's figure (not readable in the notes) | Modulus within 1 %, $\sigma<0.45$ % |
| Vibrating pin, harmonics 1 and 2 | Same assembly; $\omega_0=2\pi$ rad/s, $\varepsilon=0.2$ cm | Abeille and MGMC [`noise_vibration.yaml`](https://github.com/HunterBelanger/abeille/blob/master/input_files/noise_vibration.yaml) | Agrees with APOLLO3 within error bars with cancellation ([Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)) |
| 1-D 2-group slab | 300 cm, 60 cells, vacuum, 1 Hz, $v$ = 1.25e7 and 2.5e5 cm/s, $\lambda=0.0784$, $\beta=0.0065$; perturbation in cell 24 | FEMFFUSION [`test/1D_noise_SPN`](https://github.com/Zonni/FEMFFUSION/tree/master/test/1D_noise_SPN) | Diffusion and SPN only |
| C5G7 7-group, 5 % capture in one MOX pin, 1 Hz | 2×2 assemblies with reflector, 8 families | Cross sections from NEA/NSC/DOC(2003)16; kinetics in FEMFFUSION [`c5g7.prec`](https://github.com/Zonni/FEMFFUSION/blob/master/examples/2D_C5G7-TD/c5g7.prec) and Abeille `c5g7_td.py`; pin location only in the paper figure | NOISE-SN results ([Yi et al.](https://research.chalmers.se/publication/525456/file/525456_Fulltext.pdf)) |
| C4V 2-group, 5 % capture, 1 Hz | $v$ = 1.9919e7 and 3.8184e5 cm/s, $\beta=0.0049$, $\lambda=0.0797$ | C4V cross sections not obtained | Diffusion differs from $S_N$ by up to 16 % |
| C5G7-TD Phase I-B (noise) | Official specification | NEA working area; conditions-for-release form required ([NEA page](https://oecd-nea.org/jcms/pl_32145/deterministic-time-dependent-neutron-transport-benchmark-without-spatial-homogenisation-c5g7-td)) | — |

**CORTEX 2-group data, read from Abeille.**

| Quantity | Fuel | Water |
|---|---|---|
| $\Sigma_t$ (cm⁻¹) | 0.3779, 0.55064 | 0.25411, 1.2182 |
| $\Sigma_a$ | 0.025755, 0.15788 | 0.00079457, 0.029316 |
| $\Sigma_{s,1\to1}$, $\Sigma_{s,1\to2}$, $\Sigma_{s,2\to2}$ | 0.35128029, 0.00086471, 0.39276 | 0.22519143, 0.028124, 1.188884 |
| $\Sigma_f$ | 0.0057671, 0.10622 | — |
| $\nu_p$, $\nu_d$ (both groups) | 2.576819862, 0.013860138 | — |
| $\chi$, $\chi_d$ | [1, 0], [1, 0] | — |

One family with $\lambda=0.0851\ \text{s}^{-1}$; speeds 1.82304e7 and 4.13067e5 cm/s; $k_\text{eff}=0.99916$. The perturbed pin is the box $x\in[2.1543,2.8857]$, $y\in[-2.8857,-2.1543]$. FEMFFUSION's $\nu\Sigma_{f,1}=0.01494071$ and $\beta=0.00535$ reproduce from these numbers, so the two repositories encode the same benchmark (inference from the arithmetic).

**Cautions.**

- FEMFFUSION's 1-D file prints 1.09620e-01 in the first column for group 2 beside $\Sigma_{a,2}=1.07186\text{e-}01$. Re-read the file header before using it.
- The paper says the vibrating pin moves along $x$, two cells down and two left of centre. The shipped Abeille input vibrates along $y$ for a pin at $x=+2.52$. Check the paper figure before comparing fields. *Resolved by symmetry (amendment 2026-11-22b, Stage 8): the two set-ups are related by a 90° rotation of a fully symmetric assembly; 8c uses Abeille's.*
- The M&C 2021 multi-solver CORTEX comparison (Vinai et al.) was not read. Frequencies other than 1 and 3 Hz are unverified.
- Abeille's `rod.yaml` is not a noise input as shipped. Its link to any published rod-model case is speculative.
- Fauvel et al. Problem 2 (two-group rod, half-length 44.43 cm) is only partly recorded in the notes and cannot be set up from them.
- Experiments: the first COLIBRI campaign in CROCUS oscillated 18 rods at 0.1–2 Hz; the public record is a paper, not a data set ([Zenodo 3989875](https://zenodo.org/records/3989875)). No usable TRIGA transfer-function data set was identified.

### 3.5 Sources requested but not directly read

Everything attributed to these is second-hand or marked unverified:

- Yamamoto 2013 (Ann. Nucl. Energy 58) and 2018 (Ann. Nucl. Energy 115).
- Yamamoto and Sakamoto 2019, 2021 and 2022.
- Rouchon, Zoia and Sanchez 2017 (Ann. Nucl. Energy 102) and its HAL copy; Rouchon's 2016 thesis; Rouchon, Jarrah and Zoia, M&C 2019.
- The Pázsit and Demazière handbook chapter; Pázsit and Analytis 1980; Pázsit and Dykin 2015 (abstract only); Jonsson et al. 2012 (abstract only); Zoia et al. 2021 (abstract only).
- The CORE SIM model report; Williams; Bell and Glasstone; Thie; Pázsit and Pál.
- Booth 2003; Booth and Gubernatis 2010 (abstract only); Zhang, Lee and Lee 2016; Fauvel et al. 2024 (conference); H. Belanger's thesis.
- Vinai et al., M&C 2021; CORTEX deliverables; the FEMFFUSION manual.

Equations from Belanger et al. passed through an HTML extraction or a summarising tool. No equation was checked against a PDF by eye. Equations taken from Abeille code are the highest-confidence items. The journal reference for arXiv:2207.09822 (Nucl. Sci. Eng. 197(4), 534–557, 2023, doi:10.1080/00295639.2022.2126719) appears in two notes and was not checked by the other two.

---

## 4. OpenMC code map

Source read: `openmc-dev/openmc`, branch `develop`, commit **`013cfbb999ed42da13ec10f1948198c5b9db0765`** (2026-11-11). Read-only; nothing was built or run. All links below are permalinks at that commit. Sketches are proposals, not existing code.

### 4.1 Particle weight

**What exists.**

- Weight is one `double` everywhere: `ParticleData::wgt_`, `wgt_last_`, `wgt_born_`; `SourceSite::wgt`; `NuBank::wgt`; `TrackState::wgt` ([particle_data.h L33-L103](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/include/openmc/particle_data.h#L33-L103), [L491-L573](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/include/openmc/particle_data.h#L491-L573)).
- "Alive" is literally `wgt_ != 0.0`, and killing is "set weight to zero" ([particle_data.h L619-L635](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/include/openmc/particle_data.h#L619-L635)).
- `collision_mg` → `sample_reaction` → `create_fission_sites`, `absorption`, `scatter`, roulette ([physics_mg.cpp L26-L97](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/physics_mg.cpp#L26-L97)).
- `create_fission_sites` puts the parent weight into the expected number of sites and gives each site weight `1/weight` ([physics_mg.cpp L99-L242](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/physics_mg.cpp#L99-L242)). This is the hardest real-positive assumption to remove:

```cpp
double nu_t = p.wgt() / simulation::keff * weight * p.macro_xs().nu_fission / p.macro_xs().total;
...
site.wgt = 1. / weight;
```

- `simulation::keff` is 1.0 in fixed-source mode unless set ([simulation.cpp L334](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/simulation.cpp#L334)).
- Delayed sites get an emission time sampled from `DECAY_RATE`, with temperature and angle indices hard-coded to `0, 0` ([physics_mg.cpp L167-L179](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/physics_mg.cpp#L167-L179)).
- `absorption`: survival-biasing branch `wgt -= wgt * Σa/Σt`; analog branch kills ([physics_mg.cpp L244-L264](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/physics_mg.cpp#L244-L264)).
- `russian_roulette` compares `p.wgt()` with `<` ([physics_common.cpp L12-L36](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/physics_common.cpp#L12-L36)); weight windows do the same and split with `ceil(weight/…)` ([weight_windows.cpp L971-L1053](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/weight_windows.cpp#L971-L1053)).
- Kill-by-zero sites sit in `event_advance` (time cutoff), `event_check_limit_and_revive`, `cross_vacuum_bc` and `mark_as_lost` ([particle.cpp L277-L335](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/particle.cpp#L277-L335), [L546-L563](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/particle.cpp#L546-L563), [L737-L757](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/particle.cpp#L737-L757)).
- `SourceSite` has three hand-maintained mirrors: the MPI struct type ([initialize.cpp L175-L219](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/initialize.cpp#L175-L219)), the HDF5 bank type ([state_point.cpp L564-L596](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/state_point.cpp#L564-L596)) and the ctypes struct ([lib/core.py L23-L38](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/openmc/lib/core.py#L23-L38)).
- Event-based mode routes entirely on `alive()` ([event.cpp L20-L201](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/event.cpp#L20-L201)), so it needs no structural change.

**Minimal change.**

- Add `wgt_im_` and `wgt_im_last_` to `ParticleData`, and `wgt_im` to `SourceSite`, `NuBank` and `TrackState`, with accessors `wgt_im()` and `wgt_abs()`. Keep `wgt_` as the real part so real-mode code compiles unchanged.
- Replace the roughly eight scattered `wgt() = 0.0` kills by a `kill()` helper and an explicit `alive_` flag.
- Give `create_secondary` and `split` a `wgt_im` argument defaulting to 0.
- Roulette and weight windows act on `hypot(wgt, wgt_im)` and rescale both components by the same real factor.
- In noise mode, `create_fission_sites` takes a separate branch: Abeille's count, complex site weights, the delayed factor, no emission-time sampling.
- The copy reaction goes in `sample_reaction` next to `absorption`. The modified total goes in `Mgxs::calculate_xs`.
- `time_cutoff` must be off in noise mode.

Not checked in detail: `particle_restart.cpp`, the track, collision-track and MCPL writers, and whether DAGMC or surface-source paths compare weights (unverified).

### 4.2 Multigroup cross sections

**What exists.**

- `XsData` stores, per temperature, real tensors `total`, `absorption`, `nu_fission`, `prompt_nu_fission`, `inverse_velocity` `[angle][g]`; `decay_rate` `[angle][dg]`; `delayed_nu_fission` `[angle][dg][g]`; `chi_prompt` `[angle][gin][gout]`; `chi_delayed` `[angle][dg][gin][gout]` ([xsdata.h L66-L93](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/include/openmc/xsdata.h#L66-L93)).
- Readers accept `beta` plus `nu-fission`, or separate `prompt-nu-fission` and `delayed-nu-fission`, in vector or matrix form ([xsdata.cpp L126-L476](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/xsdata.cpp#L126-L476)).
- Access is through `Mgxs::get_xs(MgxsType, …)` with types `INVERSE_VELOCITY`, `DECAY_RATE`, `PROMPT_NU_FISSION`, `DELAYED_NU_FISSION`, `CHI_PROMPT`, `CHI_DELAYED` ([constants.h L295-L311](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/include/openmc/constants.h#L295-L311), [mgxs.cpp L424-L521](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/mgxs.cpp#L424-L521)).
- `sample_fission_energy` already selects prompt or a delayed group on real probabilities and samples the outgoing group from the matching $\chi$ ([mgxs.cpp L525-L574](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/mgxs.cpp#L525-L574)).
- `calculate_xs` caches only `total`, `absorption` and `nu_fission` ([mgxs.cpp L587-L614](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/mgxs.cpp#L587-L614)).
- `Particle::speed()` already uses the group inverse velocity in multigroup mode ([particle.cpp L55-L69](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/particle.cpp#L55-L69)).
- Histogram scattering is `ScattDataHistogram`; every `sample` ends with `wgt *= mult[gin][i_gout]` ([scattdata.cpp L531-L556](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/scattdata.cpp#L531-L556)). `multiplicity_matrix` is read when present ([xsdata.cpp L478-L540](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/xsdata.cpp#L478-L540)).
- Python setters `set_inverse_velocity`, `set_decay_rate`, `set_beta`, `set_chi_delayed` and others exist ([mgxs_library.py L707-L1036](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/openmc/mgxs_library.py#L707-L1036)).

**Three traps.**

- If `inverse-velocity` is absent, `get_xs` silently substitutes a default built from the group bounds ([mgxs_interface.cpp L240-L251](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/mgxs_interface.cpp#L240-L251)). Noise mode must raise a fatal error instead.
- `XsData::combine` **adds** density-weighted decay rates over constituents and copies inverse velocity from the first constituent only ([xsdata.cpp L589-L640](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/xsdata.cpp#L589-L640)). Use macroscopic `XSdata`, one per material, for noise work, or fix `combine`.
- The hard-coded `0, 0` at `physics_mg.cpp` L171.

**Minimal change.** Cache `inverse_velocity` in the particle's macro cross-section cache. Precompute $\Re$ and $\Im$ of $\lambda_k/(\lambda_k+i\omega)$ per material and delayed group at initialisation. No library-format change is required for forward noise.

### 4.3 Fixed-source mode and sources

**What exists.**

- A history is `sample_external_source` → `from_source` → event loop. Fission neutrons go to a per-particle LIFO secondary bank or, optionally, to a shared secondary bank processed generation by generation ([simulation.cpp L948-L965](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/simulation.cpp#L948-L965), [L980-L1094](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/simulation.cpp#L980-L1094)).
- The shared-bank loop redistributes with `synchronize_global_secondary_bank` ([bank.cpp L195-L296](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/bank.cpp#L195-L296)) and has no generation cap.
- The `max_secondaries` guard exists only in continuous-energy `physics.cpp` ([physics.cpp L120-L135](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/physics.cpp#L120-L135)). Multigroup fixed-source runs near critical have no brake.
- `Source` rejects negative strength ([source.cpp L82-L90](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/source.cpp#L82-L90)). `FileSource` and compiled sources fill a `SourceSite` directly ([source.cpp L496-L580](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/source.cpp#L496-L580)).
- Fixed-source-only settings are parsed at [settings.cpp L1233-L1291](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/settings.cpp#L1233-L1291).
- The Python source-file dtype is at [source.py L2196-L2215](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/openmc/source.py#L2196-L2215).

**Minimal change.**

- New settings `noise_mode`, `noise_frequency`, `noise_eta`, `noise_keff`, `noise_adjoint`, `noise_max_generations`, mirroring the `random_ray` dict pattern in Python ([settings.py L1447-L1475](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/openmc/settings.py#L1447-L1475)).
- Validation: multigroup only, fixed source, no time cutoff, no analog estimators, IFP off.
- Add `wgt_im` to the MPI type, the HDF5 bank type, ctypes, `SourceParticle` and the docs.
- `IndependentSource` gets an optional `phase`. Source selection stays on $|\text{strength}|$.
- `read_source_bank` matches members by name, but whether a missing `wgt_im` is tolerated is unverified ([state_point.cpp L660-L766](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/state_point.cpp#L660-L766)).
- Reuse the shared-secondary-bank generation loop as the noise fission-generation loop: fission sites go to the shared bank, copies stay in the local bank, and the rendezvous is where cancellation and the bank roulette go.

### 4.4 Tallies

**What exists.**

- `results_` is a real tensor `[filter bin][score][VALUE, SUM, SUM_SQ (, SUM_THIRD, SUM_FOURTH)]`, accumulated once per batch with `norm = total_source / (n_particles · gen_per_batch)` ([tally.cpp L829-L905](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/tallies/tally.cpp#L829-L905)).
- Multigroup scores are `score * filter_weight`, built from `p.wgt()` or `p.wgt_last()` ([tally_scoring.cpp L1621-L1640](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/tallies/tally_scoring.cpp#L1621-L1640), [L2527-L2545](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/tallies/tally_scoring.cpp#L2527-L2545)).
- Filters return bins and a real weight per bin (`FilterMatch::weights_`); `EnergyFunctionFilter` is the precedent for a filter that only multiplies ([filter_match.h L10-L16](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/include/openmc/tallies/filter_match.h#L10-L16), [filter_energyfunc.cpp L93-L103](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/tallies/filter_energyfunc.cpp#L93-L103)).
- The `higher_moments_` flag extends the third dimension from 3 to 5; that is the pattern to copy for a cross-term slot.
- Python `Tally` has `mean`, `std_dev` and `vov` but no cross-score covariance ([tallies.py L577-L652](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/openmc/tallies.py#L577-L652)).
- Trigger statistics use `std_dev/|mean|` ([trigger.cpp L30-L48](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/tallies/trigger.cpp#L30-L48)), which is meaningless for a component whose mean crosses zero.

**Minimal change.**

- A two-bin `ComplexPartFilter` (`['real','imag']`) plus component-aware weight accessors `p.wgt_tally()` and `p.wgt_last_tally()` in the roughly fifteen multigroup scoring lines. The driver calls `score_general_mg` once per component. Every existing linear score (flux, absorption, fission, inverse-velocity and so on) then becomes available as Re and Im, with no change to the statepoint format or `results_` shape.
- Add a `SUM_CROSS` slot holding $\sum_\text{batches}\Re\cdot\Im$ when the filter is present, and Python helpers `amplitude`, `phase` and their standard deviations using Section 2.7.
- Until that lands, per-batch statepoints or `openmc.lib` give the covariance by post-processing.
- Disallow analog estimators and the `WeightFilter` in noise mode.
- Scores with complex cross-section multipliers are post-processed in Python from Re/Im tallies with real multipliers.

`score_general_mg` was sampled by grep, not read in full, and the mesh-filter track-length path was not read.

### 4.5 Random-ray adjoint as a template

- The only adjoint in OpenMC is in the random-ray solver. It is implemented as `transpose_scattering_matrix()` plus `nu_sigma_f_.swap(chi_)`, followed by a normal solve ([random_ray_simulation.cpp L368-L396](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/random_ray/random_ray_simulation.cpp#L368-L396), [flat_source_domain.cpp L2349-L2368](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/random_ray/flat_source_domain.cpp#L2349-L2368)).
- It flattens only `TOTAL`, `NU_FISSION`, `FISSION`, `CHI_PROMPT`, `KAPPA_FISSION` and `NU_SCATTER`; a grep found no inverse velocity, decay rate or delayed data in `src/random_ray/` ([flat_source_domain.cpp L2093-L2169](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/random_ray/flat_source_domain.cpp#L2093-L2169)). It cannot serve as a noise reference solver as it stands.
- Monte Carlo transport has no adjoint flag in `physics_mg.cpp`, `mgxs.cpp` or `particle.cpp`. FW-CADIS is restricted to random ray ([weight_windows.cpp L788-L806](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/src/weight_windows.cpp#L788-L806)).

**Sketch for adjoint noise.** Do the transposition at library level in Python, as your existing tool does:

- Transposed scatter matrix per histogram, with the row-sum ratio carried by `multiplicity_matrix`.
- Fission in matrix form with $\chi$ and $\nu\Sigma_f$ exchanged; the matrix readers accept separate prompt and delayed matrices.
- The forward $\Sigma_t$ kept from the file.

Whether `ScattData::base_init` tolerates a transposed matrix whose normalisation no longer sums to $\Sigma_s$ was not traced (unverified). The useful conventions to borrow are the two-stage run on one simulation object and the output tagging.

### 4.6 Complex type, two doubles, or two particles

| Aspect | `std::complex<double>` weight | Two doubles `wgt`, `wgt_im` (recommended) | Two real particles |
|---|---|---|---|
| Existing real-mode code | Every `p.wgt()` comparison and `double flux = p.wgt()*d` breaks | Unchanged; noise code opts in | No struct change |
| MPI type, HDF5 bank, ctypes | Member type changes in all three; ctypes has no portable complex | One more `double` in each | None, but source files need twice the sites and a tag |
| Tallies | Real `results_` cannot take complex | Re and Im routed by a filter | Re and Im from independent histories; higher amplitude and phase variance |
| Roulette, weight windows | `std::abs` everywhere | Explicit `hypot` in noise branches | Each carrier can be negative, so sign handling remains |
| Physics | Natural | Natural | Every complex multiplier rotates weight between parts, so a Re-carrier spawns Im-carriers at each such event: bank growth and extra variance |
| GPU / SoA | Operator support inside target regions depends on compiler; unverified for the fork | Two flat `double` arrays, the form the `ParticleData` comment describes | Flat, but doubles the in-flight count |

**Use two doubles.** The `ParticleData` comment states the class exists so GPU branches can switch to a structure-of-arrays layout through accessors ([particle_data.h L459-L490](https://github.com/openmc-dev/openmc/blob/013cfbb999ed42da13ec10f1948198c5b9db0765/include/openmc/particle_data.h#L459-L490)). `std::complex<double>` remains fine for local helper arithmetic; OpenMC already uses it in `wmp.cpp`, `math_functions.cpp` and `source.cpp`.

The GPU port is a separate fork using OpenMP target offload and event-based queues ([Tramm et al.](https://arxiv.org/html/2403.12345v1)). Its concrete constraints (virtual functions, `std::vector`, multigroup or `std::complex` support on device) were **not inspected and are unverified**.

There is no delta or Woodcock tracking in OpenMC at this commit: a grep returned nothing. That matters for exact cancellation (Section 6.3).

---

## 5. Staged implementation plan with verification gates

**Rules for every stage.**

- **Build on passes only.** A stage starts only when all earlier gates pass.
- **Statistics.** Report batch statistics, and also run **20 independent replicas** (different seeds) and use the replica standard deviation for pass/fail. Batch error bars are known to be optimistic in this problem.
- **Standard pass criterion ("STD").** For each tallied quantity, $z=(\hat q-q_\text{ref})/\sigma$ is computed for Re and Im separately. Pass if $|z|<3$ for at least 99 % of entries, and the $z$ histogram is consistent with a standard normal (mean within ±0.3, standard deviation within 0.8–1.2 when there are at least 30 entries). In addition amplitude must agree within 2 % and phase within 1°, evaluated only where the replica $\sigma$ is below one third of that tolerance.
  - *Amendment (2026-11-18, Stage 4): correlated mesh tallies ("STD-C").* Mesh bins of one tally share fission chains, and in the Stage 4 slabs their replica correlation is 0.9–1.0. The z histogram of the raw bins then measures essentially one random number (the common amplitude) and fails for unbiased estimators (4a at 0.1 Hz: all $|z|<1.13$, histogram std 0.19). Normalising out the common amplitude is not enough: the shape still fluctuates in a few smooth modes, and with the measured covariance the null 99 % range of the z-histogram std is about 0.1–2.3, so no histogram test has useful power. For spatial or other strongly correlated tallies, STD-C replaces the histogram test: (i) $|z|<3$ for at least 99 % of raw bins and the amplitude/phase limits on raw bins, as in STD; (ii) $|z|<3$ for Re and Im of the integral $\sum_j v_j$ over the gated bins; (iii) a Hotelling $T^2$ test of mean = reference on the $k$ leading principal components of the replica covariance of the stacked, standardised (Re, Im) bins that explain 99 % of the variance ($k\le R-3$), with $(R-k)T^2/(k(R-1))\sim F(k,R-k)$ and pass if $p>0.0027$ (the two-sided 3σ level). The raw and shape ($v_i/\sum_j v_j$) z histograms are logged, not gated. Implemented as `harness.std_correlated`. Single-quantity tallies (Stage 3) keep STD.
- **Always log.** Per frequency: $\eta$, $k_\text{eff}$ used, number of fission generations to extinction, $\sum|w|$ per generation, bank size per generation, wall time, FOM $=1/(\sigma^2T)$ for Re and Im, and the Re–Im correlation coefficient.
- **Always plot.** Bode plots (amplitude and phase against frequency) with error bars against the reference; the $z$ histogram; $\sum|w|$ against generation on a log axis.
- **Frequency grid "F7".** 0.1, 0.3, 1, 3, 10, 30, 100 Hz.

### Stage 0 — Preflight (no OpenMC changes)

**Inputs.** This document; Abeille at `c6aaa5f`; FEMFFUSION at `1611f90`; OpenMC at `013cfbb`; your CASMO-70 library and transposition tool.

**Tasks.**

1. Open the PDF of arXiv:2207.09822 and record: Eq. 26 (arcsin or arccos); Eq. 17 (sign between the deltas); the orientation and position of the vibrating pin in the figure.
2. Search OpenMC GitHub issues, open and closed pull requests, Discussions and openmc.discourse.group for "noise", "frequency domain", "complex weight", "negative weight", "alpha eigenvalue". Record hits.
3. Build Abeille and run `input_files/noise_oscillation.yaml` unchanged. Build FEMFFUSION and run `test/1D_noise_SPN` and `examples/2D_UOX_FA`. Re-read the group-2 columns of the FEMFFUSION 1-D file.
4. Check the library: `inverse-velocity`, `decay-rate`, and delayed data (`beta` or `delayed-nu-fission`, `chi-delayed`) present for every fissile material; macroscopic `XSdata` only; record the slowest group speed.
5. Write four reference scripts (NumPy, `complex128`):
   - `ref_g0.py`: $G_0(\omega)$ for any number of families, and the cut-offs $\omega_{CC}$, $\omega_D$ for given $\eta$, $\beta$, $\Lambda$.
   - `ref_inf.py`: the multigroup infinite-medium solve of Section 3.3, forward and transposed.
   - `ref_sn1d.py`: 1-D slab multigroup $S_N$ (Gauss–Legendre, diamond difference, isotropic scattering, vacuum or reflective boundaries). It solves the static eigenproblem, then the complex fixed-source problem with the fission source scaled by **its own** $k$, by direct sparse solve or GMRES. It accepts a prescribed source, or builds $-\delta\Sigma\,\psi_0$ from its own static angular flux. It has a transposed (adjoint) mode.
   - `ref_diff1d.py`: the slab Green's function of Section 3.3.

**Gate.**

- `ref_inf.py` returns $2.486\times10^{-3}-7.439\times10^{-2}i$ for the case in Section 3.3.
- `ref_sn1d.py` with reflective boundaries and a uniform source equals `ref_inf.py` to $10^{-6}$ relative.
- `ref_sn1d.py` on a thick, weakly absorbing slab approaches `ref_diff1d.py` within 2 % away from boundaries and source.
- `ref_sn1d.py` is converged in mesh and quadrature to 0.1 % (S32 against S64, mesh halving).
- The Abeille run reproduces its shipped $k_\text{eff}=0.99916$ within statistics.

**Log.** A short file of the three by-eye answers and the search hits.

### Stage 1 — Complex weights, complex tallies and covariance; no physics change

**Change.**

- `wgt_im` plumbing (Section 4.1): structs, MPI type, HDF5 bank, ctypes, `SourceParticle`; `kill()` and `alive_`; modulus-based roulette.
- `settings.noise` block with `frequency`, defaulting so that no physics changes in this stage.
- Source `phase`; file sources with `wgt_im`.
- `ComplexPartFilter`; `SUM_CROSS` accumulator; Python `amplitude`, `phase` and their standard deviations.

Fission handling is **not** changed in this stage.

**Tests.**

- **1a, identity.** Problem P1: your existing multigroup fixed-source regression problem, or a 2-group slab 40 cm thick with CORTEX fuel data, vacuum boundaries and an isotropic group-1 source in the central 2 cm. Fission stays on if the problem is clearly subcritical; otherwise set `create_fission_neutrons = false`. Run with noise mode off and with noise mode on at $\omega=0$, same seed. The Re bins must equal the baseline tallies **bit for bit** in `sum` and `sum_sq`. The Im bins must be exactly zero.
- **1b, phase rotation.** P1 with fission off and source phase $\theta=\pi/3$. Re $=\cos\theta\times$ baseline and Im $=\sin\theta\times$ baseline to $10^{-12}$ relative. All operations are linear in $w$ and random numbers depend only on $|w|$.
- **1c, covariance.** In test 1b the Re–Im correlation must equal 1 to $10^{-10}$, and $\sigma_\varphi$ from the helper must be zero to rounding. Then mix two sources with phases 0 and $\pi/2$ in different regions, and compare the accumulated covariance with the value computed from per-batch statepoints: agreement to $10^{-12}$.
- **1d, plumbing.** Source-file round trip with `wgt_im`. Identical results for 1 and 4 MPI ranks and for 1 and 8 threads, to the level OpenMC already guarantees. Event-based against history-based. The existing regression suite unchanged with noise mode off.

**Pass.** Exact criteria as stated; no statistical tolerance in this stage.

**Log and plot.** Diff report for 1a; table of ratios for 1b; correlation values for 1c.

### Stage 2 — The $i\omega/v$ term (copy reaction), no fission

**Change.** In noise mode `calculate_xs` adds $\Sigma_\omega=\eta\,\omega/v_g$ to the total. `sample_reaction` creates a copy with probability $\Sigma_\omega/\tilde\Sigma_t$ and weight $w(1-i/\eta)$ into the local secondary bank, then applies implicit capture $w\leftarrow w\,\Sigma_s/\tilde\Sigma_t$ (forced in noise mode), norm-based roulette, and scatter. Fatal error if inverse velocity is missing. The $\Sigma_\omega$ addition must be excluded from physical reaction-rate scores.

**Tests.** Infinite medium: a 20 cm cube with reflective faces and a uniform isotropic real unit source.

- **2a, one group.** $\Sigma_t=1.0$, $\Sigma_s=0.5$, $\Sigma_a=0.5\ \text{cm}^{-1}$, $v=2.2\times10^5$ cm/s. Frequencies chosen so that $x=\omega/(v\Sigma_a)\in\{0.01,0.1,0.5,1,5\}$, that is $f=x\,v\Sigma_a/2\pi\approx175$ Hz to 87.5 kHz. These are unit-test frequencies, deliberately outside your band so the term is not negligible. Also run F7 to confirm stability where the term is tiny. $\eta\in\{0.5,1,5\}$. Reference (derivation): $\delta\phi=S/(\Sigma_a+i\omega/v)$.
- **2b, two groups.** CORTEX water data, source in group 1, same $x$ values defined with the group-2 quantities. Reference: `ref_inf.py`.
- **2c, estimators.** Collision and track-length flux agree within statistics.

**Pass.** STD, tightened to 1 % in amplitude and 0.5° in phase. Results at the three $\eta$ values mutually consistent ($|z|<3$).

**Log and plot.** Bode plot against $x$; copies per collision against $\eta$; measured within-generation growth of $\sum|w|$ per collision against $\tilde{\mathfrak c}^{cb}=(c_\infty+\sqrt{1+\eta^2}\,\omega\tau)/(1+\eta\omega\tau)$ as an informational check.

### Stage 3 — Fission with $1/k_\text{eff}$, delayed neutrons and the fission-generation loop

**Change.**

- Noise branch in `create_fission_sites`: $n=\lfloor\nu\Sigma_f/(\tilde\Sigma_tk)+\xi\rfloor$ with $k$ = `noise_keff`. Each site inherits the complex weight; delayed sites are multiplied by the precomputed $\lambda_k/(\lambda_k+i\omega)$. No emission-time sampling.
- Fission sites go to the shared secondary bank, forming the next generation. Copies stay local.
- Generation loop until the bank is empty, with `noise_max_generations` and an abort when $\sum|w|$ grows for 50 consecutive generations.
- Record $\sum|w|$ and the bank size per generation.

**Tests.** Infinite medium as in Stage 2, uniform real unit source.

- **3a, published point data.** $\Sigma_t=1$, $\Sigma_s=0$, $\Sigma_f=0.4\ \text{cm}^{-1}$, $\beta=7\times10^{-3}$, $\lambda=0.08\ \text{s}^{-1}$, $v=3\times10^5$ cm/s. $\nu$ is not recorded in the notes; choose $\nu=2.5$ and set `noise_keff` $=\nu\Sigma_f/\Sigma_a$ exactly. Frequencies F7 at $\eta=1$; for this data $f_{CC}(\eta=1)\approx807$ Hz, so all are convergent (derived). Reference: $\delta\phi=G_0(\omega)\,S/\Sigma_a$. Add the $10^4$ Hz point with source $\pi\varepsilon_f\Upsilon\Sigma_a$ and $\eta=20$, since $\eta\gtrsim15$ is needed to be convergent there without cancellation (derived). Reference: $2.49\times10^{-3}-7.44\times10^{-2}i$.
- **3b, TRIGA-like kinetics.** One group, $\Sigma_t=1.0$, $\Sigma_s=0.9$, $\Sigma_a=0.1\ \text{cm}^{-1}$, $\nu\Sigma_f=0.1$, $v=2.1097\times10^5$ cm/s so that $\Lambda=47.4\ \mu$s. Use your six-family $\beta_k$, $\lambda_k$ summing to 715 pcm. Frequencies F7 with $\eta=1$ up to 30 Hz and $\eta=3$ at 100 Hz. Reference: multi-family $G_0$ from `ref_g0.py`.
- **3c, predicted failure.** Case 3b at 100 Hz with $\eta=1$. The measured generation ratio of $\sum|w|$ must match $\tilde{\mathfrak K}^{cb}=\tilde\Upsilon/(1-\kappa\omega\Lambda)\approx1.005$ and the growth guard must fire. This confirms the cut-off estimate and tests the guard.
- **3d, two groups.** CORTEX fuel data with $k_\infty$ from `ref_inf.py`, F7, source in group 2. Reference: `ref_inf.py`.

**Pass.**

- 3a, 3b, 3d: STD.
- Measured generation ratio in the convergent cases within 1 % of $\tilde{\mathfrak K}^{cb}$; generations to extinction of order $1/\beta$.
- 3c: ratio above 1 and within 1 % of the formula; run aborted by the guard.

**Log and plot.** Bode plots of $\delta\phi\,\Sigma_a/S$ against $G_0$; $\sum|w|$ against generation; generations to extinction against frequency.

### Stage 4 — Space dependence in a slab

**Change.** None to physics. Add a helper that runs the eigenvalue calculation and passes its $k_\text{eff}$ to the noise run, and a mesh-tally post-processor.

**Tests.** Slabs with reflective transverse faces, vacuum at $x=\pm a$, isotropic scattering. Each code uses its own $k_\text{eff}$.

- **4a, one group.** $2a=60$ cm, $\Sigma_t=1.0$, $\Sigma_s=0.9$, $\nu\Sigma_f$ adjusted with `ref_sn1d.py` so that $k\approx1$, $v=2.2\times10^5$ cm/s, $\beta=0.0065$, $\lambda=0.0784\ \text{s}^{-1}$. Source: isotropic, real, unit strength in $x\in[14,15]$ cm. Frequencies 0.1, 1, 10, 100 Hz. Tally: Re and Im flux on a 1 cm mesh. Reference: `ref_sn1d.py`. Secondary: `ref_diff1d.py`, with the difference recorded and not gated; published transport–diffusion differences near a source reach 9–16 % in another problem ([Yi et al.](https://research.chalmers.se/publication/525456/file/525456_Fulltext.pdf)).
- **4b, two groups.** The FEMFFUSION `1D_noise_SPN` data after the Stage 0 re-read; otherwise CORTEX fuel data in a slab sized by `ref_sn1d.py`. 1 Hz and 10 Hz. Source: a prescribed isotropic source in one 5 cm cell. Reference: `ref_sn1d.py`; FEMFFUSION SP3 as an informational cross-check.
- **4c, point-kinetic limit.** For 4a at 1 Hz, evaluate $G_0^{dyn}$ of Section 1.5 with the static adjoint from `ref_sn1d.py` and compare with $G_0$.

**Pass.** STD on mesh cells where $|\delta\phi|>10\ \%$ of its maximum. For 4c, amplitude within 2 % and phase within 1°.

**Log and plot.** Amplitude and phase profiles against $x$ with both references; $k_\text{eff}$ of OpenMC and of the script.

**Note on $k_\text{eff}$.** The noise operator is singular at $\omega=0$, so an error $\delta k$ acts like a spurious reactivity amplified by $|G_0|$. Converge the static $k_\text{eff}$ so that $\delta k/k\ll\beta|\delta\rho|$ for your perturbation, and log the sensitivity by rerunning 4a at 0.1 Hz with $k\pm2\sigma_k$. This is a derived caution, not a published criterion.

### Stage 5 — Weight cancellation

**Change.** At each fission-generation rendezvous:

1. Bin banked fission sites on a user mesh restricted to fissionable cells. Position only, provided sites are banked before the outgoing group is sampled; otherwise bin on outgoing group as well.
2. Compute per-cell complex sum and count with thread-local accumulation and `MPI_Allreduce`.
3. Assign each site the cell mean.
4. Apply a norm-based roulette to the bank.

Options: mesh, on or off. Bank order must be canonical before binning, for reproducibility.

**Tests.**

- **5a, plumbing.** Case 3b. In an infinite medium mesh averaging is exactly unbiased. Results unchanged within statistics; log where the generation ratio falls between the physical $\mathfrak K$ and $\tilde{\mathfrak K}$.
- **5b, bias.** Case 4a at 1 Hz with cancellation off, and on with mesh sizes $h=2,1,0.5$ cm. Bin-wise $t_i=(\varphi_{i,\text{on}}-\varphi_{i,\text{off}})/\sqrt{\sigma^2_\text{on}+\sigma^2_\text{off}}$ for Re and Im, excluding bins above 20 % relative error. Repeat with the mesh origin shifted by $h/2$.
- **5c, conditionally convergent.** Case 3c (100 Hz, $\eta=1$) now terminates and meets STD. The $10^4$ Hz point of 3a at $\eta=1$ meets STD with at least $5\times10^5$ particles per generation.
- **5d, variance.** FOM on against off for 4a at 0.1 and 1 Hz.
- **5e, histogram scattering.** Repeat 5b with a 32-bin histogram scattering library for the same material (linearly anisotropic, $\bar\mu=0.3$, converted to histogram form), against the same case with cancellation off.

**Pass.**

- 5a, 5c: STD.
- 5b, 5e: $t$ histograms standard-normal at the two finest meshes; no trend with origin shift. Any failure at $h=2$ cm is recorded as the bias scale, not as a defect.
- 5d: FOM ratio above 1 at both frequencies; value recorded.
  - *Amendment (2026-11-19, Stage 5): 5d gates the quadrature component.* Mesh averaging can only remove phase dispersion among sites in a bin. In slab 4a at 0.1 and 1 Hz the field is nearly in phase with the source (|Im/Re| 0.06–0.14) and the in-phase weights all carry the same sign: $\sum|w|$ after averaging is 1.000 of the value before, and the bank roulette removes fewer than $2\times10^{-7}$ of sites. The Re variance comes from population (branching) fluctuations, which averaging does not touch; a tuning study (h = 0.5–4 cm, bank roulette on/off, 4× particles; 60 replicas × 20 batches, FOM from pooled batch variances) found Re FOM ratios 0.88–0.98 at 1 Hz and 0.99–1.12 at 0.1 Hz for every setting, against Im 2.66–2.92 and 1.60–1.76. 5d therefore passes when the 95 % bootstrap lower bound of the median Im FOM ratio exceeds 1 at both frequencies; the Re ratio is recorded and expected near 1.
- $\sum|w|$ plateaus in all runs.

**Log and plot.** $\sum|w|$ and net $|\sum w|$ against generation, on and off; FOM table; $t$ histograms; particles per occupied cell.

### Stage 6 — Noise source for an absorber of varying strength

**Change.**

- Two-stage driver: eigenvalue run to convergence; then for each noise batch `nskip−1` plain generations and one sampling generation.
- In the sampling generation, at every collision inside a perturbed region, bank complex source sites by the table of Section 2.4, with user input $\delta\Sigma_{x,g}/\Sigma_{x,g}$ (complex) per region.
- Normalise by the mean modulus and carry the factor into the tallies.
- Tally the static flux and static detector response in the same run for $\delta d/d_0$. Tally the sampled source (Re and Im).
- Default `nskip` 20.

**Tests.**

- **6a.** Infinite-medium case 3b with a uniform $\delta\Sigma_a/\Sigma_a=10^{-3}$. Reference: $\delta\phi/\phi_0=-G_0\,\delta\Sigma_a/\Sigma_a$. F7.
- **6b.** Slab 4a with $\delta\Sigma_a/\Sigma_a=0.01$ in $x\in[14,15]$ cm. Reference: `ref_sn1d.py` with source $-\delta\Sigma_a\psi_0$ from its own static angular flux. Compare $\delta\phi/\phi_0(x_\text{ref})$. Frequencies 0.1, 1, 10, 100 Hz.
- **6c, published benchmark.** CORTEX assembly oscillation at 1 Hz with the data of Section 3.4, and a 138×138 Re/Im flux tally per group. Reference: Abeille run from `noise_oscillation.yaml`, with the factor $-i\pi$ applied to map conventions and both normalised by the static group-2 flux at the assembly centre. Secondary: FEMFFUSION `2D_UOX_FA` for the pattern.

**Pass.**

- 6a, 6b: STD.
- 6c: STD against Abeille on cells above 10 % of the maximum modulus, and modulus within 2 % of Abeille, the level TRIPOLI-4 and APOLLO3 reached with each other.
- The sampled-source tally integrates to $-\delta\Sigma_a\times$ (static collision rate in the region)$/\Sigma_t$ within statistics.

**Log and plot.** 2-D amplitude and phase maps with difference maps; source tally; sensitivity to `nskip` (3 against 20).

### Stage 7 — Adjoint noise and reciprocity

*Amendment (2026-11-20, brief `noise_dev/briefs/noise_stage7_brief.md`, N4). This replaces the original Change and Tests of this stage; the original text is in git history (commit 255bb94f5).*

**Convention.** The detector signal is the bilinear form

$$\delta d(\omega)=R^\mathsf T\delta\phi(\omega)=R^\mathsf T\mathcal L(\omega)^{-1}\delta S(\omega)=\chi(\omega)^\mathsf T\delta S(\omega),\qquad \mathcal L(\omega)^\mathsf T\chi=R .$$

- The adjoint run uses the transposed library at the **same** $+\omega$, with $R$ as a real external source. Nothing is conjugated in the run.
- Users of the sesquilinear product $\langle a,b\rangle=\sum a^*b$ get $\psi^\dagger=\chi^*$. Conjugation happens only in post-processing.
- This is the bilinear convention of Section 1.4. No part of this plan runs the adjoint at $-\omega$; the $-\omega$ form there is the sesquilinear alternative.

**Change.** The transport code does not change. The adjoint is the forward noise solver run on a transposed library:

- **Scattering.** Built by the ce_mg_triga transposition, copied verbatim into `noise_dev/common/adjoint_lib.py` with a sync test.
  - The nu-scatter matrix is transposed and keeps each transfer's μ histogram.
  - Absorption becomes $\Sigma_a'=\max(0,\Sigma_t-\text{in-scatter})$.
  - The row multiplicity is $c_g=\text{in-scatter}_g/(\Sigma_t-\Sigma_a')_g$.
- **Fission and delayed data.** Added by `noise_dev/common/adjoint_noise_lib.py` as transposed production matrices (`prompt-nu-fission` $[g_\text{in},g_\text{out}]$, `delayed-nu-fission` $[k,g_\text{in},g_\text{out}]$).
  - OpenMC's existing MG matrix reader and `sample_fission_energy` then do what the original Stage 7 text derived: the family is chosen in proportion to the transposed row sums, and the outgoing group in proportion to $\nu\Sigma_{f,g}$.
  - The delayed factor $\lambda_k/(\lambda_k+i\omega)$ is a scalar per family and is applied unchanged.
  - Inverse velocity, decay rates and $\Sigma_t$ are copied.
- **Normalisation.** The adjoint source follows the ce_mg_triga convention $S_i=V_i\Sigma_{R,i}$, so that the columns of $G(\omega)$ compare with the static $G$ at $\omega\to0$.
- **Variance.** A transposed library with $c_g>1$ can have an infinite second moment: `adjoint_noise_lib.second_moment_radius` $\ge1$. Example N3x in `stage7/stage7_prep_log.md`: radius 2.05, and the estimators disagree. Test materials must have a radius below 1. Larger radii, as in the TRIGA library, need weight windows, which do not yet work in noise mode (N3c).

**Tests.**

- **7a, infinite medium, full operator.**
  - The adjoint from the N2 library against `ref_inf.py`'s adjoint (`adjoint=True`), with delayed neutrons on.
  - The 3-group material of `stage7/run_stage7_prep.py` (transposed multiplicity 2.4) and the 2-group CORTEX fuel.
  - Frequencies 0.1, 1, 10, 30 and 100 Hz.
  - Amplitude and phase. This is N3b promoted to a gate.
- **7b, slab reciprocity, external source.** A Stage 4 slab: 4b, two groups, so that the transpose is not trivial.
  - $\delta S$ is a known external source in one cell and $R$ is a detector region.
  - Compare the forward $\delta d=R^\mathsf T\delta\phi$ with the adjoint $\chi^\mathsf T\delta S$, where $\chi$ comes from an adjoint run with source $R$. Use the two-sample test on replicas.
  - Frequencies below the slab's point-kinetic break, at it and above it.
  - Low-frequency sanity check: as $\omega\to0$, with no delayed neutrons and on a subcritical fixed-source version of the slab, the adjoint noise solution approaches a static adjoint run with the same transposed library.
- **7c, reciprocity with a sampled source.** Stage 6's perturbation, $\delta S=-\delta\mathcal L\,\phi_0$.
  - **Absorber first** (capture and total only): the scalar static flux in the perturbed cells is enough for the contraction.
  - **Then a scattering perturbation.** It needs the P0 and P1 angular moments of $\phi_0$ and $\chi$. Copy the contraction verbatim from `ce_mg_triga/step5e_pool_pt.py`, with a header and the source commit, as for N1. $\chi$ is now complex; keep the same sign convention and the factor 3 in the P1 term.
- **7d (optional), $S_N$ cross-check.** If `ref_sn1d.py` solves the transposed problem, use it as an independent adjoint reference for the 7b slab.
- **7e.** Kept from the original plan: the adjoint at $-\omega$ returns the complex conjugate of the adjoint at $+\omega$, within statistics.
- **7f.** Kept from the original plan: case 7b with the histogram-anisotropic library of 5e.

**Pass.** These are the existing criteria:
- STD against the references (7a, 7d, 7e);
- the two-sample test between forward and adjoint (7b, 7c, 7f);
- STD-C where the tallies are correlated.

**Log and plot.** $\delta\phi^\dagger_g(x)$ amplitude and phase, which is the transfer map $G(\omega)$ for that detector; FOM of forward against adjoint for the detector response; penetration statistics (particles reaching $V_d$ per source particle).

  - *Amendment (2026-11-20b, Stage 7): 7e dropped; 7c contraction and pairing.*
    - **7e is not run.** The code rejects negative frequencies (`settings.cpp`, "Noise frequency must be non-negative"), and under the bilinear convention above no run needs $-\omega$. For real cross sections $\mathcal L(-\omega)=\overline{\mathcal L(\omega)}$ exactly, so the sesquilinear $\psi^\dagger=\bar\chi$ is a post-processing step; `ref_inf.py`'s self-test checks the identity. Running at $-\omega$ would need code changes (sign of the copy weight, $|\omega|$ in $\Sigma_\omega$) that this plan does not otherwise need. The Pass list's STD for 7e is void.
    - **7c uses P0 + P1 for both perturbations.** The P1 term enters through $\delta\Sigma_t$ as well as through the P1 scattering moment, so the scalar flux alone is not exact for an absorber either. On the 7b slab with the perturbed cell at $x\in[30,31]$ cm, `ref_sn1d` gives the P0-only error as 4e-4 for a thermal absorber but 25 % for a scattering density change ($\delta\Sigma_t=\delta\Sigma_s$, P0 parts nearly cancel), and P0 + P1 within 3e-4 for both, with one 1 cm bin. Both cases therefore use `pt_lib.pt_terms` (copied from ce_mg_triga at 540ab07), applied to the complex $\chi$ by linearity (`adjoint_noise_lib.pt_terms_complex`).
    - **7c is a paired test.** The contraction uses the static moments of $\phi_0$ from the forward run itself (same normalisation per static particle), so forward replica $r$ and its contraction are correlated. Each forward replica is paired with an independent adjoint replica; the differences are i.i.d. and `harness.paired_correlated` replaces the two-sample test.
    - **Penetration statistics are not logged** separately: in the 7b slab every source particle can reach the detector, and the forward/adjoint FOM ratio carries the information. They return with the out-of-core detectors of Stage 9b.

### Stage 9b preparation — weight windows with complex weights (N3c) and the shared TRIGA library (N5)

*Amendment (2026-11-22). Order: this section runs after Stage 7 and **before Stage 8**. Stage 9b does not depend on vibration, but it needs windows to reach the out-of-core detectors (the TRIGA adjoints were 79 % low from weight tails without them), and Stage 7 showed the adjoint game 10–25 times less efficient than the forward one on the N3 and CORTEX materials. N5 (checks of the regenerated ce_mg_triga library) also belongs here and runs when that library arrives. Stage 8 follows.*

**Change (N3c).** Weight windows in noise mode, which until now refuse to run (guard of 2026-11-20). In `apply_weight_window`, when `settings::noise_mode`:

- **Absolute windows.** Bounds are scaled by the noise weight scale only: 1 for an external source of unit-modulus sites, or the batch's mean source modulus for a sampled source. They are *not* renormalised by each particle's birth weight (OpenMC's `wgt_born / wgt_ww_born`): every fission site is a new particle, so per-birth renormalisation would reset the windows each generation and never control the weight spread across generations. The user supplies windows in the units of a unit source, typically lower bounds $C/|\psi^\dagger|$ from an importance estimate (CADIS).
- **Modulus.** Every comparison (weight cutoff, lower and upper bound, `max_lb_ratio`, survival weight) uses $|w|$.
- **Phase-preserving split.** Splitting into $n$ gives each copy $(\mathrm{Re}\,w/n,\ \mathrm{Im}\,w/n)$, and the survivor is scaled by $1/n$. Roulette uses the existing norm-based `russian_roulette`, which keeps the phase.
- **Fission sites** are checked against the windows when their track starts, as source particles are.
- **Roulette.** Inside the window mesh the windows replace the noise roulette. Otherwise its fixed cutoff would kill the small-weight particles that windows create near a detector. Outside the mesh the noise roulette is unchanged, and `apply_russian_roulette` (survival biasing) is not used in noise mode.
- **Unchanged.**
  - Non-noise runs follow the same code path with $\mathrm{Im}\,w=0$; their results are bit-identical.
  - Weight-window *generators* stay refused in noise mode: they tally a real flux.
- **Cancellation** acts on the fission bank between generations and windows within a generation, so the two do not interact.

**Tests.**

- **N3c-u (unit, `tests/unit_tests/test_noise.py`).** Noise mode accepts windows. A tight window that forces splitting and roulette gives the same answer as no windows within statistics. Generators are refused.
- **N3c-a, deep penetration.** A fission-free 2-group water slab (CORTEX water data), fast source at one face and a thermal detector at the other, at least 30 mean free paths apart, at two frequencies where the detector response has a large phase lag (one where Re changes sign along the slab).
  - Windows: CADIS lower bounds $C/|\chi|$ on a 1 cm × 2-group mesh, from the `ref_sn1d` transposed solution for the detector, with $C$ placing a unit source particle at its window centre; upper/lower ratio 5.
  - Runs: forward with and without windows; reference `ref_sn1d`.
- **N3c-b, adjoint efficiency.** The 7a infinite-medium adjoint cases (N3 material, CORTEX fuel) at 1, 30 and 100 Hz with group windows $C/|\psi_g|$. Here $\psi$ is the forward solution for the adjoint source $R$, i.e. the importance for $\langle\chi,R\rangle$, and $C$ places the source group's unit weight at its window centre. Runs: with and without windows; forward controls at the same particle counts (N3f for N3, a CORTEX forward run).

**Pass.**

- N3c-u: the unit tests pass.
- N3c-a:
  - STD of the windowed detector response against `ref_sn1d`;
  - the two-sample test between windows on and off, where the off run scores;
  - FOM with windows above FOM without (95 % bootstrap lower bound of the ratio above 1).
- N3c-b: STD against `ref_inf`'s adjoint with windows on. The FOM ratio with/without windows, and adjoint-with-windows against forward, are logged and not gated.

### Stage 8 — Vibrating component

**Precondition.** Stage 0 by-eye answers recorded; Stage 5 passed; cancellation on.

**Change.** A planar-interface vibration source: for harmonic $n$ at $\omega=n\omega_0$, at static collisions inside the swept layer $|x-x_0|\le\varepsilon$, bank source sites by the Section 2.4 rules with $\delta\Sigma_x=\Delta\Sigma_x\,c_n(x)$, $\Delta=$ (negative-side material) − (positive-side material) from the two materials' group data. Collisions are sampled in the material actually present, and the weights use the ratio of $\delta\Sigma_x$ to that material's $\Sigma_x$. Overlapping vibrations sum.

**Tests.**

- **8a.** A 1-D slab with a fuel–water interface inside the core region of case 7c, $\varepsilon=0.2$ cm, $\omega_0=2\pi$ rad/s, harmonics $n=1,2$. Reference: `ref_sn1d.py` with the source built from the same $c_n(x)$ on a mesh resolving the layer (at least 40 cells across $2\varepsilon$).
- **8b.** $\varepsilon\to0$ check: the $n=1$ response scales linearly with $\varepsilon$ for $\varepsilon=0.05,0.1,0.2$ cm.
- **8c, published benchmark.** Vibrating pin in the CORTEX assembly at 1 Hz, harmonics 1 and 2, cancellation mesh 170×170 and 340×340. Reference: Abeille `noise_vibration.yaml`, in its own orientation (pin at $(+2.52,-2.52)$ moving along $y$; amendment 2026-11-22b). Compare with the published statement of agreement with APOLLO3 within error bars.

**Pass.**

- 8a: STD.
- 8b: linear within statistics.
- 8c: STD against Abeille on cells above 10 % of maximum modulus; results at the two cancellation meshes mutually consistent.
- The real part of the first harmonic must not show the seed-dependent bias reported for runs without cancellation.

**Log and plot.** Source tally (Re and Im) across the layer against $c_n$; maps per harmonic; FOM with cancellation; statement of which $c_n$ form was used.

  - *Amendment (2026-11-22b, Stage 8): pin orientation and the second estimator.*
    - **Orientation resolved by symmetry.**
      - The paper text has the pin at $(-2.52,-2.52)$ cm (two cells left of and two below centre) moving along $x$.
      - Abeille's `noise_vibration.yaml` has the pin at $x\in[2.1543,2.8857]$, $y\in[-2.8857,-2.1543]$ (centre $(+2.52,-2.52)$) moving along $+y$: two flat-vibration layers of half-width 0.2 cm on the pin's $y$ faces, positive side water above and fuel below.
      - The rotation $(x,y)\to(-y,x)$ maps the first onto the second, including the sense of motion. The CORTEX assembly is invariant under it: the 17×17 lattice with 25 water holes (`stage6/cortex.py`) is unchanged by 90° rotation and by every reflection, and both positions are fuel pins. The pins, water blade, reflective boundaries and data are identical.
      - So the two set-ups give the same field up to that rotation. 8c uses Abeille's configuration, which carries the reference numbers.
    - **Statistics.** The vibration source is two layers of nearly equal and opposite weight (fuel ↔ water). Stage 7c showed that a forward sampled source of this cancelling kind can have a σ of 40–150 %. Stage 8 starts by measuring the forward σ of each gated quantity on a short run and sizing the gate runs from it, with cancellation on as the precondition requires.
    - **Second estimator in the gate.** Each gated detector response is also computed as the adjoint contraction $\langle\chi,\delta S\rangle$.
      - $\chi$ is from an adjoint run (transposed library, Stage 7) with the detector as source.
      - $\delta S=-\delta\mathcal L\,\phi_0$ is built from tallies on the swept layers: the static $\phi_0$ moments and the material difference times $c_n$, contracted with `pt_lib.pt_terms` (P0 + P1) as in 7c.
      - Gate: the forward estimate and the contraction agree (paired or two-sample test as in 7c), in addition to the 8a–8c criteria.

### Stage 9 — Many-group benchmark and TRIGA dry run

**Change.** None, except input tooling for 70 groups and a cancellation mesh aligned with fuel elements.

**Tests.**

- **9a.** C5G7 2-D noise (7 groups, 8 families, 5 % capture oscillation in one MOX pin, 1 Hz) against FEMFFUSION SP3 (`examples/2D_C5G7_NN`) and the published NOISE-SN description. Informational on amplitude near the source, where diffusion-like solvers differ; gated on the global pattern and phase.
- **9b.** TRIGA model with the CASMO-70 histogram library. Compute cut-offs from the IFP values. Select $\eta$ per frequency: 1 up to 30 Hz, 3 at 100 Hz. A central absorber of varying strength; forward detector responses and one adjoint run per detector; F7.

**Pass.**

- 9a: phase within 2° and amplitude within 5 % of FEMFFUSION SP3 away from the perturbed assembly; differences near the source recorded.
- 9b: forward–adjoint reciprocity within statistics, using a source region and spectrum for which the adjoint tally is exact, as in 7c. $\sum|w|$ bounded. Results independent of cancellation mesh at two sizes. On the plateau $\delta d/d_0$ within about 10 % of $G_0\,\delta\rho$ for a central absorber, as a sanity check and not a gate.

**Log and plot.** Bode plots of $\delta d/d_0$ per detector, forward and adjoint; $G(\omega)$ maps; cost per frequency; memory high-water mark.

**Changes from the suggested order.**

- A preflight stage is added.
- The generation loop and guards enter with delayed neutrons (Stage 3), because near-critical chains need them before cancellation exists.
- Norm-based roulette and covariance are in Stage 1.
- The absorber source (Stage 6) is separated from vibration (Stage 8), which now follows cancellation.
- The adjoint (Stage 7) precedes vibration because it does not depend on it and serves your main goal.
- The published benchmarks are attached to the stages they verify (6c, 8c), with a final many-group stage.

---

## 6. Risks and open questions

### 6.1 Variance near critical and at low frequency

- The extended generation radius is about $1-\beta$ on the plateau and tends to 1 as $\omega\to0$. The game is nearly critical and "very hard to converge" at low frequency ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)). At 0.1 Hz your band is above $\omega_{LF}=\lambda_\text{eff}$ but close to the bottom of the plateau. The MGMC thesis excluded frequencies below 0.1 Hz for convergence reasons ([Ferrero thesis](https://webthesis.biblio.polito.it/41032/1/tesi.pdf)).
- No FOM-against-frequency data exists in any source read, and none below 1 Hz for a heterogeneous system.
- The top of the band is conditionally convergent at $\eta=1$ (Section 2.3). The cut-off values are point-kinetics estimates. No convergence analysis exists for a realistic multigroup heterogeneous core, and Fauvel et al. list this as open.
- Sensitivity to the static $k_\text{eff}$ grows as $|G_0|$ (Stage 4 note).
- Batch error bars are optimistic; replicas are required.

### 6.2 Memory and runtime

- Published runs are large: $10^4$–$10^5$ batches or replicas of $10^4$–$10^6$ particles for a 2-group 2-D assembly at 1 Hz. Absolute runtimes are almost never given; the only figure is a cap of two days on 512 cores ([Belanger et al.](https://ar5iv.labs.arxiv.org/html/2207.09822)).
- Each source particle drives of order 140 fission generations on your plateau, and each generation needs a global rendezvous.
- Seven frequencies, several detectors and one adjoint run per detector per frequency multiply the cost.
- The mesh-averaging `MPI_Allreduce` design avoids the master gather but is unpublished.
- Multigroup fixed-source OpenMC has no secondary-bank guard today. The guard in Stage 3 is required, not optional.
- A 70-group library with 32-bin histograms per transfer is static memory you already carry. The noise additions are two doubles per particle and site, plus the cancellation mesh.

### 6.3 Histogram angular scattering and weight cancellation

- **Mesh-averaged cancellation of fission-born particles is unaffected** by the scattering representation. Only the position, and possibly the outgoing group, of the banked site matters (inference from the published algorithm). Tests 5e and 7e check this.
- **Exact regional cancellation is affected.** It needs the scattering-cosine density evaluated toward every point of the region and strictly positive there ([Belanger et al. 2022](https://ar5iv.labs.arxiv.org/html/2205.12062)). The following points are inferences from the notes:
  - OpenMC's order (outgoing group, then $\mu$ from that transfer's histogram) matches the required factorisation.
  - Any zero-probability histogram bin makes $\zeta=0$ on a cone, so those particles must take $\beta=0$.
  - The histogram is discontinuous, so the minimum of $\zeta$ is not at a region corner. Use the sampled $\langle\zeta\rangle$, $\langle1/\zeta\rangle$ rule.
  - OpenMC has only surface tracking, so each $\zeta$ needs an optical-path ray trace.
  - Fissions induced on a source particle's first flight, or after a copy, cannot take part.
- No demonstration of exact cancellation with anisotropic scattering exists in the sources read.
- If $\chi_d\ne\chi_p$, position-only binning requires banking sites before the outgoing group is sampled. The published papers do not discuss this choice.
- Cancelling scattered particles would need direction- and group-aware bins and "was not investigated" ([Fauvel et al. 2025](https://doi.org/10.1051/epjn/2025062)).

### 6.4 What makes continuous energy much harder

- The delayed term needs per-nuclide, per-precursor yields and spectra rather than $\beta_k$ and $\chi_{d,k,g}$.
- Vibration sources need the fictitious union material and nuclide-wise weights (Belanger Eqs. 30–34); in multigroup $\Delta\Sigma$ is a matrix difference.
- Exact cancellation is "straightforward in multi-group codes but harder in continuous-energy codes", which lack lab-frame conditional scattering distributions ([Belanger et al. 2022](https://ar5iv.labs.arxiv.org/html/2205.12062)).
- Fission in OpenMC's continuous-energy path has its own weight assumptions, photon coupling and thermal scattering. None of these were mapped for noise.
- A continuous-energy adjoint does not exist in OpenMC. Your transposed-library route is multigroup only.

### 6.5 Out-of-core detectors, the adjoint, and items to check first

**Deep penetration.** Cancellation reduces the variance of the in-core complex fission source but does nothing for penetration statistics through graphite and water. No source was found on weight windows or CADIS-type methods for complex weights, or on out-of-core noise detectors with Monte Carlo. A modulus-based window that preserves phase is a design hypothesis. OpenMC's weight-window generators need a real flux, and FW-CADIS runs only with random ray. This is the main argument for the adjoint route, where the detector is the source.

**Adjoint.** No published adjoint noise Monte Carlo was found, so Stage 7 has no external reference beyond reciprocity and the $S_N$ script. Two specific open points:

- The fission branching rule in the transposed game is a derivation.
- Applying $G(\omega)=\delta\psi^\dagger$ to a real perturbation needs $\langle\delta\psi^\dagger,-\delta\Sigma_t\psi_0\rangle$, a product of **angular** fluxes. Scalar-flux tallies give it exactly only where one factor is isotropic. Otherwise it needs angular moments or bins, with an error that must be quantified against a forward run. No source addresses this.

**Model limits.** The orthodox linearization fails for mechanical vibrations at realistic amplitudes, over-predicting higher harmonics. A reference that shares the linearization verifies the code, not the physics.

**Items to check by eye or by search before coding.**

| Item | Why it matters | Stage |
|---|---|---|
| Belanger Eq. 26: arcsin or arccos | Determines $c_n$; transcription and Abeille code disagree | 0, 8 |
| Belanger Eq. 17: sign between the deltas | Phase of the oscillation source at $-\omega_0$ | 0, 6 |
| Vibrating pin position and axis | Paper text and shipped Abeille input disagree; equivalent by symmetry, Abeille's used (amendment 2026-11-22b) | 0, 8 |
| OpenMC issues, open pull requests, Discussions, Discourse | Not searched; prior or parallel work is unknown | 0 |
| GPU fork data-layout constraints | Unverified; affects only the weight-representation argument | before upstreaming |
| FEMFFUSION 1-D group-2 column values | Printed numbers look inconsistent | 0, 4 |
| How MGMC and Abeille cull the bank after averaging; prompt against delayed in the cancellation bank | Not established by the notes | 5 |
| Counting and zero-power noise formulas | No primary source read | before interpreting measured spectra |
| `read_source_bank` tolerance of a missing `wgt_im` | Backward compatibility of source files | 1 |

---

## References

**Read (full text or large parts, some through HTML extraction or a summarising tool)**

1. H. Belanger, D. Mancusi, A. Rouchon, A. Zoia, "Variance reduction and noise source sampling techniques for Monte Carlo simulations of neutron noise induced by mechanical vibrations", arXiv:2207.09822. https://arxiv.org/abs/2207.09822 (HTML: https://ar5iv.labs.arxiv.org/html/2207.09822). Journal reference as recorded in two notes: Nucl. Sci. Eng. 197(4), 534–557 (2023), doi:10.1080/00295639.2022.2126719.
2. A. Fauvel, A. Rouchon, D. Mancusi, A. Zoia, "Convergence of Monte Carlo algorithms for power reactor noise", EPJ Nuclear Sci. Technol. 11, 72 (2025). https://doi.org/10.1051/epjn/2025062 (about two thirds read; Eq. 68 details and Table 2 not retrieved by one note).
3. H. Belanger, D. Mancusi, A. Zoia, "Exact weight cancellation in Monte Carlo eigenvalue transport problems", Phys. Rev. E 104, 015306 (2021), doi:10.1103/PhysRevE.104.015306. https://arxiv.org/abs/2103.13891
4. H. Belanger, D. Mancusi, A. Rouchon, A. Zoia, "Unbiasedness and optimization of regional weight cancellation", Phys. Rev. E 106, 025302 (2022) (journal reference as given in the task, not confirmed on the arXiv page). https://arxiv.org/abs/2205.12062
5. A. Rouchon, M. Vadée Le Brun, A. Zoia, "Analysis and comparison of APOLLO3 and TRIPOLI-4 neutron noise solvers", PHYSOR 2020, EPJ Web Conf. 247, 21002 (2021). https://doi.org/10.1051/epjconf/202124721002
6. A. Rouchon, R. Sanchez, I. Zmijarevic, A. Zoia, CEA slides on deterministic and stochastic noise solvers (2017). https://www.cea.fr/energies/tripoli-4/Documents/8_rouchon_manu.pdf (several equations garbled).
7. P. Vinai et al., "On the simulation of neutron noise induced by vibrations of fuel pins in a fuel assembly", Ann. Nucl. Energy 181, 109521 (2023), doi:10.1016/j.anucene.2022.109521. https://research.chalmers.se/publication/533337/file/533337_Fulltext.pdf
8. O. Ferrero, "Non-perturbative computation of neutron noise using the Monte Carlo method", MSc thesis, Politecnico di Torino with CEA (2025/2026). https://webthesis.biblio.polito.it/41032/1/tesi.pdf
9. C. Demazière, A. Tatidis, P. Vinai, "Monte Carlo-based dynamic calculations of stationary perturbations", PHYSOR 2020, EPJ Web Conf. 247, 21003, doi:10.1051/epjconf/202124721003. https://research.chalmers.se/publication/524729/file/524729_Fulltext.pdf
10. X. Yang, L. Cao, Q. Zheng, Q. He, H. Wu, "A new Monte-Carlo-based iterative method for neutron noise calculation", RPHA 2023. https://www.kns.org/files/int_paper/paper/RPHA-2023_2023_2/RPHA-2023-XuranYang.pdf
11. Hansen and Hansson, "On the feasibility of neutron noise diagnostics in fast Gen-IV reactors … using linear multi-group diffusion theory", MSc thesis, Chalmers (2025). https://odr.chalmers.se/server/api/core/bitstreams/c88f0408-122b-40a7-b838-be90e2d0a4ec/content
12. Pázsit group, "Research and development program in reactor diagnostics and monitoring with neutron noise methods", report to the Swedish regulator (stage and year not confirmed). https://www.stralsakerhetsmyndigheten.se/contentassets/b1a81ffee6da4b9d970b54b302109cfb/0028-research-and-development-program-in-reactor-diagnostics-and-monitoring-with-neutron-noise-methods
13. Salih, "Response of an in-core neutron detector to a vibrating absorber based on the detector adjoint function", MSc thesis, Iowa State University (1980). https://behost.lib.iastate.edu/DR/Salih_ISU-1980-Sa33.pdf (scanned; equations partly garbled).
14. C. Demazière, "Core diagnostics using noise analysis: from proof-of-principle to industrial demonstration", M&C 2017. https://www.kns.org/files/int_paper/paper/MC2017_2017_6/P304S06-10DemaziereC.pdf
15. C. Demazière, "Open source neutron noise modelling", IAEA Technical Meeting, Milano (2022). https://research.chalmers.se/publication/530864/file/530864_Fulltext.pdf
16. H. Yi, P. Vinai, C. Demazière, "On the simulation of neutron noise using a discrete ordinates method". https://research.chalmers.se/publication/525456/file/525456_Fulltext.pdf
17. Coupled CORE SIM user's manual. https://publications.lib.chalmers.se/records/fulltext/204034/204034.pdf
18. T. Yamamoto, KAKENHI final report 17K07011 (publication list and summary only). https://kaken.nii.ac.jp/file/KAKENHI-PROJECT-17K07011/17K07011seika.pdf
19. J. Tramm et al., "Performance portable Monte Carlo particle transport on Intel, NVIDIA, and AMD GPUs", arXiv:2403.12345. https://arxiv.org/html/2403.12345v1
20. Zywiec, Ingegno, Heinrichs, Forget, "Rossi-alpha benchmark validation of a static alpha eigenvalue capability in OpenMC", arXiv:2606.00907 (abstract page). https://arxiv.org/abs/2606.00907

**Code and data read directly**

21. OpenMC, `openmc-dev/openmc`, branch `develop`, commit `013cfbb999ed42da13ec10f1948198c5b9db0765`. https://github.com/openmc-dev/openmc/tree/013cfbb999ed42da13ec10f1948198c5b9db0765 ; documentation https://docs.openmc.org/en/latest/usersguide/kinetics.html and https://docs.openmc.org/en/latest/methods/random_ray.html
22. Abeille, commit `c6aaa5fd5a53c076a4d8db6d412e90768fa16b87`, GPL-3.0. https://github.com/HunterBelanger/abeille
23. MGMC, commit `1231c27`, CeCILL-2.1, Zenodo doi:10.5281/zenodo.4585368 (README and inputs only). https://github.com/HunterBelanger/mgmc
24. FEMFFUSION, commit `1611f90`, GPL-3.0; reference paper doi:10.1016/j.anucene.2023.110239 (not read). https://github.com/Zonni/FEMFFUSION ; https://femffusion.webs.upv.es/
25. SCONE, commit `a4a02bf` (grep only). https://github.com/CambridgeNuclear/SCONE
26. CEA TRIPOLI-4 distribution page. https://www.cea.fr/energies/tripoli-4/distribution
27. OECD/NEA C5G7-TD benchmark page. https://oecd-nea.org/jcms/pl_32145/deterministic-time-dependent-neutron-transport-benchmark-without-spatial-homogenisation-c5g7-td
28. COLIBRI first campaign paper, Zenodo doi:10.5281/zenodo.3989875. https://zenodo.org/records/3989875

**Abstract or record only**

29. A. Zoia, A. Rouchon, B. Gasse, C. Demazière, P. Vinai, Ann. Nucl. Energy 154, 108061 (2021), doi:10.1016/j.anucene.2020.108061. https://research.chalmers.se/en/publication/522444
30. I. Pázsit, V. Dykin, "The dynamic adjoint as a Green's function", Ann. Nucl. Energy 86, 29–34 (2015), doi:10.1016/j.anucene.2015.01.029. https://research.chalmers.se/en/publication/229494
31. A. Jonsson, H. N. Tran, V. Dykin, I. Pázsit, Kerntechnik 77(5), 371–380 (2012), doi:10.3139/124.110258. https://research.chalmers.se/en/publication/168848
32. CORE SIM report record. https://research.chalmers.se/en/publication/154380 ; CORE SIM+ record https://research.chalmers.se/en/publication/523604
33. T. Booth, J. Gubernatis, Nucl. Sci. Eng. 165(3), 283–291 (2010), doi:10.13182/NSE09-62. https://ans.org/pubs/journals/nse/article-10762
34. Grablevsky, H. Belanger, "Distributed memory algorithms for weight cancellation in Monte Carlo particle transport simulations", EPJ Web Conf. 302, 09007 (2024), doi:10.1051/epjconf/202430209007.
35. Saliba et al., "Development of an event tracking feature in OpenMC for neutron noise analysis", EPJ Web Conf. 338, 04001 (2025), doi:10.1051/epjconf/202533804001.
36. P. Vinai et al., M&C 2021 multi-solver comparison. https://research.chalmers.se/publication/547937 ; slides doi:10.5281/zenodo.5575801

**Not directly read (cited second-hand; contents unverified)**

37. T. Yamamoto, "Monte Carlo method with complex-valued weights for frequency domain analyses of neutron noise", Ann. Nucl. Energy 58, 72–79 (2013), doi:10.1016/j.anucene.2013.03.002.
38. T. Yamamoto, Ann. Nucl. Energy 115, 494–501 (2018), doi:10.1016/j.anucene.2018.02.008.
39. T. Yamamoto, H. Sakamoto, Prog. Nucl. Energy 117, 103098 (2019), doi:10.1016/j.pnucene.2019.103098; Prog. Nucl. Energy 145, 104130 (2022), doi:10.1016/j.pnucene.2022.104130; J. Nucl. Sci. Technol. 58, 190–200 (2021).
40. A. Rouchon, A. Zoia, R. Sanchez, "A new Monte Carlo method for neutron noise calculations in the frequency domain", Ann. Nucl. Energy 102, 465–475 (2017), doi:10.1016/j.anucene.2016.11.035 (HAL cea-02509654).
41. A. Rouchon, PhD thesis, Université Paris-Saclay (2016).
42. I. Pázsit, C. Demazière, "Noise techniques in nuclear systems", Handbook of Nuclear Engineering (2010).
43. I. Pázsit, G. Analytis, Ann. Nucl. Energy 7, 171 (1980).
44. C. Demazière, CORE SIM model description, Ann. Nucl. Energy 38, 2698 (2011); report CTH-NT-241.
45. T. Booth, Nucl. Sci. Eng. 143, 291–300 (2003), doi:10.13182/NSE02-10TN.
46. Zhang, Lee, Lee, J. Comput. Phys. 305, 387–402 (2016), doi:10.1016/j.jcp.2015.10.042.
47. A. Fauvel et al., "Convergence of Monte Carlo methods for neutron noise", EPJ Web Conf. 302, 08004 (2024), doi:10.1051/epjconf/202430208004.
48. M. M. R. Williams, *Random Processes in Nuclear Reactors*; G. Bell, S. Glasstone, *Nuclear Reactor Theory*; J. Thie, *Power Reactor Noise*; I. Pázsit, L. Pál, *Neutron Fluctuations*.
49. NEA/NSC/DOC(2003)16, C5G7 benchmark specification. https://www.oecd-nea.org/science/docs/2003/nsc-doc2003-16.pdf (link taken from the FEMFFUSION README; not opened).
