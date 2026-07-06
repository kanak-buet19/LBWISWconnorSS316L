# Modelling Wire / Powder Directed Energy Deposition by Extending laserbeamFoam

**A CFD/multiphase-VOF literature knowledge base and code-modification workflow**

*Scope.* This document is a working knowledge base for adapting the open-source OpenFOAM solver **laserbeamFoam** ([Flint et al. 2023](https://doi.org/10.1016/j.softx.2022.101299); v2.0 [Flint et al. 2024](https://doi.org/10.1016/j.softx.2023.101612)) from its native laser powder-bed / melting configuration into a **directed energy deposition (DED)** solver — either **wire-fed (wire-DED / laser hot-wire)** or **blown-powder DED**. It has three parts: (I) a map of the existing solver's governing equations and source-code structure, so you know exactly what you are modifying; (II) a synthesis of the modelling literature — what interface-capturing schemes, heat-source models, and mass-addition strategies the field uses, each anchored to a primary paper; and (III) a concrete, file-by-file modification workflow with a reasoned recommendation on which process (wire vs powder) to attempt first.

The central technical problem is stated once and returned to throughout: **laserbeamFoam is a fixed-mass solver, and DED is fundamentally a mass-addition process.** Every architectural decision below flows from how you choose to add feedstock mass, momentum, and enthalpy to the computational domain.

---

## Part I — What laserbeamFoam already is

### I.1 Physical model and governing equations

laserbeamFoam is an **incompressible, two-phase, one-fluid Volume-of-Fluid (VOF)** solver. Two immiscible phases are tracked: the **metal** (substrate/melt) and the **shielding gas / atmosphere**, distinguished by a phase fraction $\alpha_1 \in [0,1]$ (`alpha.metal`). A second scalar, the **liquid fraction** $\varepsilon_1 \in [0,1]$ (`epsilon1`), tracks the melted state of the metal via an enthalpy-porosity formulation. The solver couples five physics blocks, all visible in the source:

**Phase advection (VOF).** The metal fraction is transported by either algebraic MULES or geometric isoAdvector:
$$\frac{\partial \alpha_1}{\partial t} + \nabla\cdot(\mathbf{u}\,\alpha_1) + \nabla\cdot\big(\mathbf{u}_r\,\alpha_1(1-\alpha_1)\big) = S_\alpha .$$
The compression term (third on the LHS) and the interface-sharpening machinery are the classical VOF apparatus. In the current code $S_\alpha=0$ — there is no mass source — but note that the MULES routine `alphaEqn.H` already carries an `Su`/`Sp` source slot (`Su + fvm::Sp(Sp + divU, alpha1)`). **This slot is the single most important entry point for DED mass injection.**

**Momentum (single-fluid, shared velocity).** One velocity field $\mathbf{u}$ is shared by both phases:
$$\frac{\partial(\rho\mathbf{u})}{\partial t} + \nabla\cdot(\rho\mathbf{u}\mathbf{u}) = -\nabla p_{rgh} - (\mathbf{g}\cdot\mathbf{x})\nabla\rho + \nabla\cdot\boldsymbol{\tau} + \mathbf{f}_{\sigma} + \mathbf{f}_{Ma} + \mathbf{f}_{rec} - K_{D}\,\mathbf{u} + \mathbf{f}_{L} .$$
Reading `UEqn.H`, the right-hand side assembles, term by term:
- **Surface tension** $\mathbf{f}_\sigma$ via the continuum-surface-force (CSF) model, `mixture.surfaceTensionForce()`.
- **Marangoni (thermocapillary) stress** $\mathbf{f}_{Ma} = \varepsilon_1\,\dfrac{d\gamma}{dT}\,\big(\nabla T - \hat{\mathbf{n}}(\hat{\mathbf{n}}\cdot\nabla T)\big)\,|\nabla\alpha_1|$, active only where $T>T_{liquidus}$ (the `Marangoni_Constant` term).
- **Recoil pressure** $\mathbf{f}_{rec} = p_{vap}\,\nabla\alpha_1$, with $p_{vap} = 0.54\,p_0\exp\!\big(\tfrac{L_v M_m (T-T_{vap})}{R\,T\,T_{vap}}\big)$ (the `pVap` field) — the vapour-driven surface depression that produces keyholing.
- **Darcy mushy-zone drag** $-K_D\mathbf{u}$ with the Carman–Kozeny form $K_D = C\dfrac{(1-\varepsilon_1)^2}{\varepsilon_1^3 + \epsilon_0}$ (`DarcyConstantlarge`, `DarcyConstantsmall`) — this is what freezes solid metal in place.
- **Lorentz force** $\mathbf{f}_L = \mathbf{J}\times\mathbf{B}$, optional, supplied by the `MTHD` (magneto-thermo-hydrodynamic) module for arc/EM-assisted processes.
- A **buoyancy** density factor $\rho_k = 1 - \varepsilon_1\beta(T-T_{solidus})$ (Boussinesq-type thermal expansion).

**Energy (enthalpy-porosity).** From `TEqn.H`, temperature is solved with latent heat handled by an iterative liquid-fraction update:
$$\frac{\partial(\rho c_p T)}{\partial t} + \nabla\cdot(\rho c_p \mathbf{u} T) = \nabla\cdot(\kappa\nabla T) + \dot{q}_{laser} + \dot{q}_{Joule} - L_f\Big(\frac{\partial(\rho\varepsilon_1)}{\partial t}+\nabla\cdot(\rho\mathbf{u}\varepsilon_1)\Big) + \dot{q}_{rad} - \dot{q}_{vap},$$
where the latent-heat sink is the `TRHS` term, $\dot q_{rad}$ is linearised Stefan–Boltzmann radiation, and $\dot q_{vap}$ is the evaporative cooling flux $Q_v = 0.82\,\dfrac{L_v M_m p_0}{\sqrt{2\pi M_m R T}}\exp\!\big(\tfrac{L_v M_m (T-T_{vap})}{R T T_{vap}}\big)$. After each solve, the liquid fraction is relaxed toward the local solidus/liquidus band: $\varepsilon_1 \leftarrow \varepsilon_1 + \dfrac{\varepsilon_{rel} c_p}{L_f}(T - T_{corr})$, clipped to $[0,1]$.

**Laser heat source (ray tracing).** $\dot q_{laser}$ (`laser.deposition()`) is not a surface flux prescribed by hand — it is computed by a **Lagrangian ray-tracing** module (`laserHeatSource`, `laserRayParticle`). Rays are launched from the beam, reflected/absorbed at the metal–gas interface using a Fresnel absorptivity model, and deposit energy along their path (multiple bounces resolve keyhole absorption). The module exposes a `powderSim()` switch that, when on, treats a **pre-packed DEM powder bed** (individual particles resolved on the mesh) rather than a flat substrate — this is why the repository can already render ray tracing through a powder bed.

### I.2 Source-code map (what lives where)

| Path | Role | DED-relevant? |
|---|---|---|
| `applications/solvers/laserbeamFoam/laserbeamFoam.C` | Main PIMPLE loop; phase→props→laser→U→T→p ordering | Add feed logic here |
| `.../UEqn.H` | Momentum assembly (Marangoni, recoil, Darcy, Lorentz) | Add feed momentum |
| `.../TEqn.H` | Energy + enthalpy-porosity latent-heat iteration | Add feed enthalpy |
| `.../updateProps.H` | Mixture properties, `alpha_filtered`, `epsilon1mask` (powder masking) | Extend for feed material |
| `.../MULES/alphaEqn.H` | **VOF transport with `Su`/`Sp` source slot** | **Primary mass-injection point** |
| `.../createFields.H` | Field & material-property declarations | Declare feed fields |
| `src/laserHeatSource/` | Ray-tracing heat source; `powderSim_` switch; `laserRayParticle` | Powder-DED beam attenuation |
| `src/MTHD/` | Magneto-thermo-hydrodynamics (arc/EM Lorentz force) | Wire-arc (WAAM) variants |
| `src/geometricVoF/` | isoAdvector geometric interface advection | Keep; sharp free surface |
| `tutorials/` | `laserbeamFoam`, `compressiblelaserbeamFoam`, `laserMeltFoam`, `multiComponentlaserbeamFoam` | Templates to fork |

The multi-component solver (`multiComponentlaserbeamFoam`) and its `laserThreePhaseMixture` are especially relevant: they already generalise the VOF machinery beyond two phases, which is the natural home for adding a *distinct feedstock phase* if you go the Eulerian route.

---

## Part II — Literature review

The modelling literature for DED and its cousins (WAAM, laser welding, powder-bed fusion) is large, but for the purpose of extending a VOF melt-pool solver it organises into five questions: what numerical building blocks the melt-pool physics rests on; how the free surface is captured; how the laser–feedstock energy coupling is modelled for powder vs wire; how feedstock **mass** is injected; and what the laserbeamFoam lineage itself has already demonstrated. The through-line is that mass addition, not heat input, is the hard part — the heat-source and melt-pool machinery in laserbeamFoam is already competitive with the published state of the art.

### II.1 Numerical foundations the melt-pool model rests on

Four papers underpin essentially every solver in this space, laserbeamFoam included, and are worth reading first because the modification work reuses their formulations verbatim. The one-fluid VOF method itself originates with [Hirt & Nichols 1981](https://doi.org/10.1016/0021-9991(81)90145-5), whose scalar advection of a fractional volume function remains the backbone of free-surface AM simulation. Surface tension enters as a volumetric body force through the continuum-surface-force model of [Brackbill et al. 1992](https://doi.org/10.1016/0021-9991(92)90240-Y), which is exactly the `surfaceTensionForce()` term in `UEqn.H` and the reason curvature is computed from $\nabla\alpha$. Melting and solidification without an explicit front are handled by the fixed-grid **enthalpy-porosity** method of [Voller & Prakash 1987](https://doi.org/10.1016/0017-9310(87)90317-6), whose mushy-zone momentum sink $K_D=C(1-\varepsilon)^2/(\varepsilon^3+\epsilon_0)$ is reproduced term-for-term in the solver's Darcy drag; the choice of the constant $C$ (`DarcyConstantlarge`) is a known sensitivity, quantified recently for keyhole AM by the Voller–Prakash-constant study surfaced in the laserbeamFoam citation network. Finally, the sharpness of the captured interface — decisive for a slender wire bead or a resolved powder particle — is set by the advection scheme; laserbeamFoam offers the geometric **isoAdvector** method of [Roenby et al. 2016](https://doi.org/10.1098/rsos.160405), which advects a reconstructed sharp interface (PLIC-like isosurface) and gives markedly less numerical smearing than algebraic MULES on the same mesh. For DED, where the deposited-bead geometry *is* the answer you are after, that interface fidelity is not a luxury.

### II.2 One-fluid vs two-fluid, and the interface-capturing choice

Before committing to laserbeamFoam's one-fluid VOF formulation it is worth understanding what it trades away. [Zenz et al. 2025](https://doi.org/10.1007/s43939-025-00434-0) provide a direct, controlled comparison of one-fluid and two-fluid approaches for laser-induced melt-pool formation and vaporisation, and the practical conclusion is that the shared-velocity one-fluid model — which laserbeamFoam uses — is adequate and far cheaper for melt-pool dynamics up to moderate vaporisation, while a two-fluid (separate metal/vapour momentum) treatment becomes necessary only when the vapour jet and its back-pressure on the melt are themselves the object of study. The same group's compressible **Mass-of-Fluid** solver ([Zenz et al. 2023](https://doi.org/10.1016/j.compfluid.2023.106109)) shows the route to a fully compressible vapour phase if keyhole/vapour dynamics later become central; laserbeamFoam's own `compressibleLaserbeamFoam` variant is the in-family analogue. For a first DED extension the recommendation from this literature is unambiguous: stay one-fluid and incompressible, and spend the modelling budget on mass addition rather than on vapour compressibility.

### II.3 Blown-powder DED: beam attenuation, catchment, and particle dynamics

Powder-DED modelling divides on how the powder stream is represented. The high-fidelity end resolves individual particles and their interaction with the beam. [Aggarwal et al. 2024](https://doi.org/10.1016/j.addma.2024.104344) couple high-fidelity modelling with high-speed imaging to quantify **laser-beam shadowing** and in-flight powder-particle dynamics — the attenuation of the beam by the powder cloud and the resulting redistribution of energy — which is precisely the physics laserbeamFoam's ray tracer must reproduce for a blown-powder configuration. [Khairallah et al. 2023](https://doi.org/10.1016/j.addma.2023.103684) present a high-fidelity DED model of the laser–powder–melt-pool interaction and show how the incident beam profile propagates through to solidification microstructure, establishing that beam-shape and powder-catchment fidelity are first-order for track geometry. At the continuum end, [Zhang et al. 2021](https://doi.org/10.1016/j.ijthermalsci.2021.106954) model heat transfer in the melt pool and **clad generation** in DED of stainless steel by treating powder addition as a distributed mass/energy source on the free surface rather than resolving particles — computationally the cheaper path and a good template for a first Eulerian implementation. Single-track geometry prediction with experimental validation for 316L-Si by [Biyikli et al. 2022](https://doi.org/10.1007/s12540-022-01243-3) demonstrates the multi-physics-plus-regression workflow used to tie such a model to measured bead dimensions, and [Zhu et al. 2025](https://doi.org/10.1016/j.optlastec.2024.111732) extend the flow-and-temperature-field simulation to a hollow (annular) laser beam, relevant if you later model coaxial-nozzle beam shapes. Process-mapping and defect prediction with a full multi-physics melt-pool model — the endgame for a validated solver — is exemplified for refractory C103 by [Prabhune et al. 2025](https://doi.org/10.1016/j.jmapro.2025.09.082).

The key architectural insight from this body of work: laserbeamFoam's existing `powderSim()` mode (DEM-packed particles + ray tracing through the bed) is already a *particle-resolved* representation. Converting it to blown-powder DED means (a) making particles arrive continuously from a nozzle rather than pre-packing a static bed, and (b) melting-and-merging them into the pool — both of which are mass-addition problems, not heat-source problems.

### II.4 Wire-DED and laser hot-wire: the primary target

Wire feeding is geometrically simpler than a powder cloud — a single continuous solid entering the pool — but the melt-transfer physics (continuous bridge vs dripping) is subtle, and this is the most directly relevant literature for the intended modification. The essential reference is [Wei et al. 2018](https://doi.org/10.1016/j.ijheatmasstransfer.2018.04.164), a **comprehensive model of transport phenomena in laser hot-wire deposition** that resolves the coupled wire heating, melting, liquid-bridge formation, and pool dynamics — it is the closest published analogue to what an extended laserbeamFoam should produce, and its treatment of the wire as a heated, melting solid region feeding the pool is the template for the momentum/enthalpy source terms. [Hu et al. 2019](https://doi.org/10.1080/13621718.2019.1591039) focus specifically on the **thermal-fluid dynamics of liquid-bridge transfer** in laser wire deposition, distinguishing the stable-bridge regime (smooth deposition) from dripping and stubbing — the transfer-mode map any wire-DED solver must be able to reproduce. The experimentally-validated multiphysics wire-DED model of [Nagaraja et al. 2024](https://doi.org/10.1007/s00170-024-14905-w) is the most recent full CFD treatment with validation and is the natural benchmark case; its reference list is effectively a curated bibliography for this exact problem. Earlier, [Nie et al. 2016](https://doi.org/10.1016/j.jmatprotec.2016.04.006) combined experiment and modelling for H13-steel laser hot-wire AM, and [Demir 2018](https://doi.org/10.1016/j.optlaseng.2017.07.003) demonstrated micro laser metal wire deposition, both useful for material properties and process windows. On the process-comparison side, [Li et al. 2022](https://doi.org/10.1016/j.ijmachtools.2022.103942) contrast high-deposition-rate powder- and wire-based laser DED directly, and [Bambach et al. 2021](https://doi.org/10.1016/j.addma.2021.102269) compare Inconel-718 deposition using powder, cold wire, and hot wire — jointly they frame the powder-vs-wire decision in measured terms. A thermal-fluid-mechanical model coupling melt flow to residual-stress evolution for wire DED is given by [Liang et al. 2021](https://doi.org/10.1016/j.jmapro.2021.08.008), pointing to the eventual thermomechanical extension.

### II.5 Arc-based deposition (WAAM) as the mass-injection exemplar

Wire-arc AM solvers solved the *droplet/continuous mass-injection into a VOF pool* problem years ago, and their numerics transfer directly even though the heat source differs (arc, not laser). The 3D heat-transfer/fluid-flow/electromagnetic CMT-WAAM model of [Cadiou et al. 2020a](https://doi.org/10.1016/j.addma.2020.101541) and the companion droplet-generation-and-pool model [Cadiou et al. 2020b](https://doi.org/10.1016/j.ijheatmasstransfer.2019.119102) show how molten-metal droplets are introduced as a VOF phase with prescribed mass, momentum, and enthalpy and then merged into the pool — this is the exact `Su` source-term pattern you will implement, and the electromagnetic (Lorentz) coupling maps onto laserbeamFoam's `MTHD` module. Layer-by-layer heat-and-mass-transfer WAAM modelling by [Li et al. 2024](https://doi.org/10.1016/j.amf.2024.200159) demonstrates the multi-layer build-up bookkeeping (adding material layer on layer) that a DED solver needs for anything beyond a single track. The metal-transfer physics itself — droplet detachment frequency and momentum — is grounded in the GMAW modelling literature: [Xu et al. 2009](https://doi.org/10.1016/j.ijheatmasstransfer.2008.09.018) model arc plasma and metal transfer in three dimensions, [Wang et al. 2003](https://doi.org/10.1088/0022-3727/36/9/313) analyse metal-transfer modes, and [Zhou et al. 2022](https://doi.org/10.1016/j.jmapro.2022.07.063) characterise CMT-based transfer behaviour for WAAM. A broad orientation to the wire-arc DED field is given by the state-of-the-art review of [Costello et al. 2023](https://doi.org/10.1080/0951192x.2022.2162597). For coaxial wire-powder hybrid feeding — an advanced configuration — see [Zhou et al. 2023](https://doi.org/10.1115/1.4062216) and the recent multiphysics-plus-microstructure model of [Kong et al. 2025](https://doi.org/10.1016/j.mtla.2025.102461). A CFD feasibility study of laser-wire DED for on-orbit (microgravity) manufacturing by [Noori Rahim Abadi et al. 2022](https://doi.org/10.3389/frspt.2022.880012) is a clean, self-contained OpenFOAM-adjacent wire-DED CFD case worth reading for its boundary-condition setup.

### II.6 The laserbeamFoam lineage and its demonstrated extensions

Finally, the solver's own development trajectory shows what is straightforward to build on it. The original release ([Flint et al. 2023](https://doi.org/10.1016/j.softx.2022.101299)) established the ray-tracing-plus-VOF-plus-state-transition architecture; version 2.0 ([Flint et al. 2024](https://doi.org/10.1016/j.softx.2023.101612)) added thermally-induced state transitions and the powder-bed capability. Downstream, the framework has been extended to **dissimilar-material chemical mixing and molten-pool shape** — the `alpha.Mn` species-transport results in the repository come from this line of work — demonstrating that adding transported scalar species and multiple metal phases is a proven modification pattern, not speculative. Related high-fidelity multiphysics melt-pool models that inform parameter choices include the single/multiple-track defect-mechanism study of [Yan et al. 2017](https://doi.org/10.1016/j.actamat.2017.05.061), the high-temperature laser-absorption model of [Yang et al. 2021](https://doi.org/10.1002/adem.202100137) (relevant to the ray-tracer's Fresnel absorptivity), and the internal-flow-behaviour study of [Ebrahimi et al. 2022](https://doi.org/10.1016/j.matdes.2022.110385).

---

## Part III — Modification workflow

### III.1 The one decision that governs everything: how to add mass

DED adds feedstock; laserbeamFoam conserves mass. There are three published strategies for injecting feedstock into a VOF melt-pool solver, in increasing order of fidelity and effort. Choosing among them is the first architectural decision, and it should be made before any code is touched.

**Strategy A — Eulerian volumetric source ("virtual mould / source-term" method).** Add a source $S_\alpha$ to the VOF equation in a small "feed region" above the substrate, converting gas cells to metal cells at a rate set by the feed velocity and cross-section, with matching momentum and enthalpy sources in `UEqn.H`/`TEqn.H`. This is the [Zhang et al. 2021](https://doi.org/10.1016/j.ijthermalsci.2021.106954) / continuum-clad approach and the [Cadiou et al. 2020a](https://doi.org/10.1016/j.addma.2020.101541) droplet-source pattern. It reuses the existing `Su` slot in `alphaEqn.H`, needs no new mesh motion, and is by far the fastest to a working single-track simulation. Its limitation is that the incoming feedstock geometry is imposed, not resolved.

**Strategy B — Resolved feedstock phase (geometric VOF).** Introduce the wire (or arriving droplets/particles) as an actual region of $\alpha_1=1$ that is advected into the domain with a prescribed inlet velocity, letting isoAdvector resolve the liquid bridge / droplet coalescence. This is the fidelity of [Wei et al. 2018](https://doi.org/10.1016/j.ijheatmasstransfer.2018.04.164) and [Hu et al. 2019](https://doi.org/10.1080/13621718.2019.1591039) for wire, and is the only way to capture transfer-mode transitions (bridge → dripping → stubbing). More expensive and requires careful inlet boundary conditions on a patch representing the wire tip.

**Strategy C — Lagrangian particle feed (powder only).** Extend the existing DEM/`laserRayParticle` machinery so powder particles are injected continuously from a nozzle, heated in flight by the ray tracer, and deposited/melted on impact — the [Aggarwal et al. 2024](https://doi.org/10.1016/j.addma.2024.104344) beam-shadowing picture. Highest fidelity for blown powder, but couples the Lagrangian cloud to the Eulerian pool, which is the most invasive change to the codebase.

The recommendation, developed below, is **wire-DED via Strategy A first, then B**. Wire is a single continuous inlet rather than a stochastic particle cloud, so the mass bookkeeping is deterministic and one-dimensional; Strategy A gets you a running, inspectable single-track deposition on the existing mesh with the smallest diff; and the wire-DED validation literature (Nagaraja, Wei, Hu) gives concrete benchmark geometries. Powder-DED via Strategy C is the more spectacular demonstration but multiplies the moving parts (nozzle model + in-flight heating + catchment efficiency + cloud–pool coupling), and its mass injection is stochastic — a harder first target.

### III.2 Governing-equation additions

Whichever strategy, the conservation equations gain matched source terms so that mass, momentum, and energy are added *consistently* (adding metal without its enthalpy will spuriously chill the pool; adding it without momentum will stall the bead). For a feed of mass rate $\dot m$ over a feed-region volume $V_f$, distributed by a normalised spatial kernel $G(\mathbf{x})$ (support = feed region, $\int G\,dV = 1$):

$$S_\alpha = \frac{\dot m}{\rho_1}\,G(\mathbf{x}) \qquad\text{(VOF: gas}\to\text{metal conversion)}$$
$$\mathbf{S}_U = \dot m\,\mathbf{u}_{feed}\,G(\mathbf{x}) \qquad\text{(momentum of incoming feedstock, }\mathbf{u}_{feed}\text{ = wire/nozzle velocity)}$$
$$S_T = \dot m\big[c_p(T_{feed}-T_{ref})\big]G(\mathbf{x}) \qquad\text{(sensible enthalpy of feedstock at its delivery temperature)}$$

For wire, $\dot m = \rho_1 A_{wire} v_{wire}$ with $A_{wire}=\pi d_{wire}^2/4$; $\mathbf u_{feed}$ points along the wire axis at the wire-feed speed; $T_{feed}$ is ambient for cold wire or the preheat temperature for hot wire (the [Wei 2018](https://doi.org/10.1016/j.ijheatmasstransfer.2018.04.164) hot-wire case adds a resistive $\dot q_{Joule}$, which laserbeamFoam's `MTHD` Joule term can supply). The enthalpy-porosity machinery then melts the added metal exactly as it melts the substrate — no special melting logic is needed for the feedstock, which is the elegance of reusing this solver.

### III.3 File-by-file change list (Strategy A, wire-DED)

1. **`createFields.H`** — declare the feed fields and read feed parameters from a new `DEDProperties` dict: `vector feedVelocity`, `scalar feedMassRate` (or `wireDiameter` + `wireFeedSpeed`), `scalar feedTemperature`, and a `feedRegion` (a `topoSet`/`cellSet` or an analytic mask marking where feedstock enters). Construct a `volScalarField Sfeed` (the kernel $G$) and, if hot-wire, a Joule-source scalar.
2. **`readControls.H`** — read/refresh any time-varying feed parameters (feed on/off, ramp).
3. **`MULES/alphaEqn.H`** (and the isoAdvector `alphaEqnSubCycle`) — populate the existing `Su` term with $S_\alpha$ over the feed region. This is the mass injection. Verify the MULES limiter still bounds $\alpha_1\in[0,1]$ with a positive source.
4. **`UEqn.H`** — add $\mathbf S_U$ as an explicit source (or `fvOptions` momentum source) so injected metal carries its feed momentum; keep it inside the feed region only.
5. **`TEqn.H`** — add $S_T$ so injected metal carries its enthalpy; for hot-wire, route the resistive preheat through the existing `Joulesource` hook.
6. **`updateProps.H`** — ensure the feed region's material properties (density, $c_p$, $\kappa$, melting range) are set; if depositing a *dissimilar* alloy, follow the existing `alpha.Mn` species-mixing pattern from the multi-component solver.
7. **`laserbeamFoam.C`** — move the beam (and feed region) along the scan path each timestep; the beam is already time-addressable via `timeVsLaserPower_`, so add a coupled feed-position update. For multi-layer builds, add the layer-increment bookkeeping from [Li et al. 2024](https://doi.org/10.1016/j.amf.2024.200159).
8. **`Make/options`, `Make/files`** — no new libraries needed for Strategy A (only if you add a bespoke `fvOption` or a new mixture class). Copy the solver to `wireDEDFoam` before editing so the upstream solver stays intact.
9. **Tutorial case** — fork `tutorials/laserbeamFoam` into `tutorials/wireDEDFoam`: add `constant/DEDProperties`, define the substrate + feed region in `blockMeshDict`/`topoSetDict`, set the scan path, and (for Strategy B later) refine the mesh at the wire tip. Validate single-track bead height/width against a Nagaraja or Wei benchmark.

### III.4 Verification ladder (do these in order)

1. **Mass conservation check.** With the laser off and feed on, integrate $\int\rho\alpha_1\,dV$ over time; it must increase at exactly $\dot m$. This isolates the `Su` term before any thermal coupling.
2. **Static melting.** Feed off, laser on a bare substrate — reproduce a known laserbeamFoam melt pool to confirm you have not broken the baseline.
3. **Single-track deposition.** Feed + laser + scan on a 2D or thin-3D slice; check that a bead forms, that its cross-section is stable, and that Marangoni/recoil are behaving. Compare bead width/height to [Biyikli et al. 2022](https://doi.org/10.1007/s12540-022-01243-3) or [Nagaraja et al. 2024](https://doi.org/10.1007/s00170-024-14905-w).
4. **Transfer mode (Strategy B).** Only once A works, resolve the wire tip and look for the bridge/dripping transition of [Hu et al. 2019](https://doi.org/10.1080/13621718.2019.1591039).
5. **Multi-track / multi-layer.** Add layer bookkeeping; check inter-track fusion and residual heat accumulation.

### III.5 Pitfalls flagged by the literature and the source

- **Enthalpy consistency.** The dominant beginner error is adding $S_\alpha$ without $S_T$: the pool then cools unphysically as cold "metal" appears. Always add all three sources together.
- **Darcy constant.** The mushy-zone constant `DarcyConstantlarge` strongly affects pool depth and bead shape; the laserbeamFoam citation network's Voller–Prakash-constant sensitivity study shows it must be tuned, not left at a default.
- **Interface smearing.** Prefer isoAdvector over MULES for the deposited bead — the [Roenby et al. 2016](https://doi.org/10.1098/rsos.160405) geometric scheme keeps the free surface sharp, which matters when the bead geometry is the deliverable.
- **One-fluid limits.** Per [Zenz et al. 2025](https://doi.org/10.1007/s43939-025-00434-0), do not over-interpret vapour-jet or spatter behaviour in the incompressible one-fluid model; if that becomes central, migrate to `compressibleLaserbeamFoam`.
- **Timestep / Courant.** Recoil pressure and surface tension impose severe capillary timestep limits already; adding a fast feed inlet can tighten them further. Use the solver's adaptive `setDeltaT` and expect small timesteps.
- **Powder catchment (if you go Strategy C).** Catchment efficiency is well below 100% and beam shadowing is real ([Aggarwal et al. 2024](https://doi.org/10.1016/j.addma.2024.104344)); a naive "all powder lands and melts" assumption will over-predict deposition.

### III.6 Suggested milestone sequence

1. Fork solver → `wireDEDFoam`; compile unchanged to confirm toolchain.
2. Implement Strategy A mass source in `alphaEqn.H`; pass the mass-conservation check (III.4.1).
3. Add momentum + enthalpy sources; pass single-track deposition (III.4.3).
4. Fork and validate a tutorial against a published wire-DED bead.
5. (Optional) Strategy B resolved wire tip for transfer-mode physics.
6. (Optional) Species transport for dissimilar-alloy DED, reusing the `alpha.Mn` pattern.
7. (Later) Thermomechanical coupling for residual stress, per [Liang et al. 2021](https://doi.org/10.1016/j.jmapro.2021.08.008).

---

## Reference list

All DOIs below were resolved and checked against CrossRef.

**Numerical foundations** — Hirt & Nichols 1981, VOF method, [10.1016/0021-9991(81)90145-5](https://doi.org/10.1016/0021-9991(81)90145-5) · Brackbill et al. 1992, continuum surface force, [10.1016/0021-9991(92)90240-Y](https://doi.org/10.1016/0021-9991(92)90240-Y) · Voller & Prakash 1987, enthalpy-porosity, [10.1016/0017-9310(87)90317-6](https://doi.org/10.1016/0017-9310(87)90317-6) · Roenby et al. 2016, isoAdvector, [10.1098/rsos.160405](https://doi.org/10.1098/rsos.160405)

**laserbeamFoam lineage** — Flint et al. 2023, laserbeamFoam, [10.1016/j.softx.2022.101299](https://doi.org/10.1016/j.softx.2022.101299) · Flint et al. 2024, v2.0, [10.1016/j.softx.2023.101612](https://doi.org/10.1016/j.softx.2023.101612)

**One- vs two-fluid / compressible VOF** — Zenz et al. 2025, one- vs two-fluid, [10.1007/s43939-025-00434-0](https://doi.org/10.1007/s43939-025-00434-0) · Zenz et al. 2023, Mass-of-Fluid, [10.1016/j.compfluid.2023.106109](https://doi.org/10.1016/j.compfluid.2023.106109)

**Powder-DED** — Aggarwal et al. 2024, beam shadowing/particle dynamics, [10.1016/j.addma.2024.104344](https://doi.org/10.1016/j.addma.2024.104344) · Khairallah et al. 2023, high-fidelity laser-powder-pool, [10.1016/j.addma.2023.103684](https://doi.org/10.1016/j.addma.2023.103684) · Zhang et al. 2021, clad generation, [10.1016/j.ijthermalsci.2021.106954](https://doi.org/10.1016/j.ijthermalsci.2021.106954) · Biyikli et al. 2022, 316L-Si single track, [10.1007/s12540-022-01243-3](https://doi.org/10.1007/s12540-022-01243-3) · Zhu et al. 2025, hollow-laser DED, [10.1016/j.optlastec.2024.111732](https://doi.org/10.1016/j.optlastec.2024.111732) · Prabhune et al. 2025, C103 process mapping, [10.1016/j.jmapro.2025.09.082](https://doi.org/10.1016/j.jmapro.2025.09.082)

**Wire-DED / laser hot-wire** — Wei et al. 2018, transport in laser hot-wire, [10.1016/j.ijheatmasstransfer.2018.04.164](https://doi.org/10.1016/j.ijheatmasstransfer.2018.04.164) · Hu et al. 2019, liquid-bridge transfer, [10.1080/13621718.2019.1591039](https://doi.org/10.1080/13621718.2019.1591039) · Nagaraja et al. 2024, validated wire-DED multiphysics, [10.1007/s00170-024-14905-w](https://doi.org/10.1007/s00170-024-14905-w) · Nie et al. 2016, H13 laser hot-wire, [10.1016/j.jmatprotec.2016.04.006](https://doi.org/10.1016/j.jmatprotec.2016.04.006) · Demir 2018, micro laser wire deposition, [10.1016/j.optlaseng.2017.07.003](https://doi.org/10.1016/j.optlaseng.2017.07.003) · Li et al. 2022, powder vs wire high-rate DED, [10.1016/j.ijmachtools.2022.103942](https://doi.org/10.1016/j.ijmachtools.2022.103942) · Bambach et al. 2021, IN718 powder/cold/hot wire, [10.1016/j.addma.2021.102269](https://doi.org/10.1016/j.addma.2021.102269) · Liang et al. 2021, thermal-fluid-mechanical stress, [10.1016/j.jmapro.2021.08.008](https://doi.org/10.1016/j.jmapro.2021.08.008) · Noori Rahim Abadi et al. 2022, on-orbit laser-wire CFD, [10.3389/frspt.2022.880012](https://doi.org/10.3389/frspt.2022.880012)

**WAAM / arc mass-injection & metal transfer** — Cadiou et al. 2020a, CMT-WAAM 3D model, [10.1016/j.addma.2020.101541](https://doi.org/10.1016/j.addma.2020.101541) · Cadiou et al. 2020b, droplet generation & pool, [10.1016/j.ijheatmasstransfer.2019.119102](https://doi.org/10.1016/j.ijheatmasstransfer.2019.119102) · Li et al. 2024, layer-by-layer WAAM, [10.1016/j.amf.2024.200159](https://doi.org/10.1016/j.amf.2024.200159) · Xu et al. 2009, arc plasma & metal transfer, [10.1016/j.ijheatmasstransfer.2008.09.018](https://doi.org/10.1016/j.ijheatmasstransfer.2008.09.018) · Wang et al. 2003, metal-transfer modes, [10.1088/0022-3727/36/9/313](https://doi.org/10.1088/0022-3727/36/9/313) · Zhou et al. 2022, CMT transfer behaviour, [10.1016/j.jmapro.2022.07.063](https://doi.org/10.1016/j.jmapro.2022.07.063) · Costello et al. 2023, wire-arc DED review, [10.1080/0951192x.2022.2162597](https://doi.org/10.1080/0951192x.2022.2162597)

**Coaxial wire-powder hybrid** — Zhou et al. 2023, coaxial wire-powder feeding, [10.1115/1.4062216](https://doi.org/10.1115/1.4062216) · Kong et al. 2025, multiphysics + microstructure, [10.1016/j.mtla.2025.102461](https://doi.org/10.1016/j.mtla.2025.102461)

**Melt-pool physics & absorption** — Yan et al. 2017, single/multi-track defects, [10.1016/j.actamat.2017.05.061](https://doi.org/10.1016/j.actamat.2017.05.061) · Yang et al. 2021, high-temperature laser absorption, [10.1002/adem.202100137](https://doi.org/10.1002/adem.202100137) · Ebrahimi et al. 2022, internal melt-pool flow, [10.1016/j.matdes.2022.110385](https://doi.org/10.1016/j.matdes.2022.110385)




