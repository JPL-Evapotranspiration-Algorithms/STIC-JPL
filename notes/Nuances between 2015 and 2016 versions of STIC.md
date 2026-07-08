# Comparative Analysis of Mallick et al. (2015) and Mallick et al. (2016)

The papers **Mallick et al. (2015)** (published in *Water Resources Research*) and **Mallick et al. (2016)** (published in *Hydrology and Earth System Sciences*) represent successive milestones in the development of the **Surface Temperature Initiated Closure (STIC)** framework. This framework physically integrates thermal infrared remote sensing—specifically radiometric surface temperature ($T_R$)—into the Penman-Monteith (PM) equation to estimate land surface fluxes and conductances analytically without relying on empirical submodels. 

While both papers share the foundational goal of removing exogenous parameterizations for aerodynamic and surface conductances, they differ significantly in model architecture, flux granularity, geographic scope, and ecohydrological objectives.

---

## 1. Core Nuances and Methodological Evolution

### A. Flux Granularity: Lumped vs. Partitioned Components
* **Mallick et al. (2015):** This study operates on a single-source "big-leaf" perspective to calculate the **total, lumped** terrestrial latent heat flux ($\lambda E$) and sensible heat flux ($H$). The conductances derived are bulk, aggregated terms: aerodynamic conductance ($g_B$) and surface conductance ($g_s$), where $g_s$ represents the combined soil and canopy pathways under partial cover conditions.
* **Mallick et al. (2016):** This paper advances the framework to **explicitly partition** the total latent heat flux into its sub-components: **Transpiration ($\lambda E_T$)** and **Evaporation ($\lambda E_E$)** (comprising soil and canopy interception evaporation), fulfilling an explicit partitioning necessity noted in the earlier framework. To achieve this, it establishes an analytical structure where moisture availability ($M$) is leveraged as a direct "evaporation efficiency index" to isolate the components. Correspondingly, the conductances are specifically isolated as canopy-scale aerodynamic conductance ($g_A$) and canopy surface conductance ($g_C$).

### B. Underlying Model Architecture and State Equations
* **Mallick et al. (2015):** The model relies on the simultaneous solution of **four state equations** ($g_B, g_s, T_0, \Lambda$) by combining the PM equation with the Advection-Aridity (AA) complementary hypothesis. It iteratively updates the Priestley-Taylor parameter ($\alpha$) to find a stable value under limiting ecohydrological conditions.
* **Mallick et al. (2016):** Designated as **STIC1.2**, this version introduces a structural breakthrough by integrating the Penman-Monteith model with the **Shuttleworth-Wallace (SW) dual-source framework**. Rather than approximating in-canopy states simplistically, STIC1.2 establishes a recursive feedback loop to analytically estimate the precise vapor pressure ($e_0$) and saturation vapor pressure ($e_0^*$) *at the canopy source/sink height* (roughness length stream), directly resolving previous challenges in retrieving aerodynamic vapor pressure at the $T_0$ level. Furthermore, it derives a fully physical analytical equation for a dynamically variable $\alpha$ based on the coupled PM-SW equations.

### C. Conceptualization and Treatment of Hysteresis
* **Mallick et al. (2015):** Hysteresis ($\lambda E - D_A - T_R$) is identified as a **mathematical error source** causing the original STIC model to systematically overestimate $\lambda E$ under dry conditions. The paper resolves this by introducing a dual-regime retrieval for moisture availability ($M$): switching from a standard surface wetness equation to a root-zone moisture availability equation derived from Granger and Gray (1989) when an afternoon decoupling between net radiation ($R_N$), $T_R$, and vapor pressure deficit ($D_A$) is detected.
* **Mallick et al. (2016):** Hysteresis is examined as an **eco-physiological trait**. Utilizing the decoupling coefficient ($\Omega$), the paper evaluates how diurnal hysteresis loops between $g_C$ and $\lambda E_T$ reflect biological adaptations to water supply-demand constraints across different ecosystems, mapping how these loops widen or narrow based on soil water reserves and rooting depths.

---

## 2. Geographic Focus and Ecohydrological Application

The shift between these two publications also marks a transition from wide-scale model validation to a highly targeted regional ecological diagnostic study:

| Feature | Mallick et al. (2015) (WRR) | Mallick et al. (2016) (HESS) |
| :--- | :--- | :--- |
| **Geographic Scope** | Global / Multi-Experimental | Regional Ecosystem (Amazon Basin) |
| **Data Sources** | 4 distinct international field campaigns: SMEX02 (Iowa), BEAREX08 (Texas), FIFE (Kansas), SAFARI2000 (Zambia/Botswana). | 6 long-term flux tower sites from the Large-scale Biosphere-Atmosphere Experiment in Amazonia (LBA). |
| **Ecosystems Covered** | Intensive agro-ecosystems (corn, soybeans, cotton), prairie grasslands, and open savannas. | Diverse tropical forest classes: Tropical Rainforest (TRF), Tropical Moist Forest (TMF), Tropical Dry Forest (TDF), and converted Pastures (PAS). |
| **Core Scientific Question** | Can a non-parametric model utilizing thermal remote sensing match or exceed the accuracy of heavily parameterized surface energy balance models across diverse climates? | How do climate anomalies (e.g., the 2005 Amazon drought) and deforestation (forest-to-pasture conversion) alter canopy-atmosphere coupling, water-use efficiency, and $\lambda E_E / \lambda E_T$ partitioning? |

---

## 3. Key Findings and Diagnostic Insights

### Mallick et al. (2015) Insights
* **Wind and Roughness Independence:** Proved that omitting explicit wind speed data introduces a relatively minor error (~4.5%), because the mechanical effects of wind are implicitly captured via the lower boundary signature embedded within $T_R$, $R_N$, and $G$ measurements.
* **Atypical Boundary Layers:** Uncovered that the primary failure point of the model occurred under large-scale horizontal advective conditions (e.g., irrigated fields adjacent to hot, dry contrasting surfaces in BEAREX08), where standard energy closures cannot fully resolve non-local boundary layer developments.

### Mallick et al. (2016) Insights
* **Regime Shifts in Tropical Forests:** Demonstrated that during the wet season, tropical forests are entirely radiation-driven ($R_N$ dictates 75–80% of transpiration variance), exhibiting weak biophysical coupling. In contrast, during dry seasons and severe drought phases, control shifts drastically to a soil-moisture-driven biophysical regime (explaining 50–65% of transpiration variance).
* **The Deforestation Paradox:** Revealed that despite massive structural differences in aerodynamic roughness ($g_A$) between towering rainforests and short pastures, both ecosystems exhibit remarkably similar canopy-atmosphere coupling during the dry season. This occurs because the pasture experiences severe soil-moisture depletion due to shallow root systems, forcing a steep biological reduction in $g_C$ that effectively mirrors the coupled state of the forest.