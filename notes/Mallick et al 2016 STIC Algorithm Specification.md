# STIC (Surface Temperature Initiated Closure) Algorithm Specification
**Mallick et al. HESS 2016 (Canopy-Scale Biophysical Controls of ET Components)**

## 1. Introduction
This document specifies the STIC framework as used in Mallick et al. (2016, Hydrology and Earth System Sciences) to analyze canopy-scale controls of transpiration and evaporation in the Amazon Basin.

The 2016 study uses the STIC analytical closure lineage introduced in Mallick et al. (2014) and extended with the radiometric iterative update framework in Mallick et al. (2015). In that sense, the 2016 contribution is not a new algebraic closure, but a physically based application and interpretation of STIC conductances and ET partitioning across plant functional types, seasons, and drought states.

---

## 2. Physical Variables and Constants

### 2.1 Core State Variables
- $R_N$: net radiation (W m$^{-2}$)
- $G$: ground heat flux (W m$^{-2}$)
- $\phi = R_N - G$: available energy (W m$^{-2}$)
- $T_R$: radiometric surface temperature (degC)
- $T_A$: air temperature (degC)
- $RH$: relative humidity (fraction)
- $e_A$: actual vapor pressure (hPa)
- $D_A$: atmospheric vapor pressure deficit (hPa)
- $s$: slope of saturation vapor pressure curve (hPa K$^{-1}$)

### 2.2 Closure Unknowns Retrieved by STIC
- $g_A$: aerodynamic conductance (m s$^{-1}$)
- $g_C$: canopy or aggregated surface conductance (m s$^{-1}$)
- $T_0$: aerodynamic source or sink temperature (degC)
- $e_0$: vapor pressure at source or sink height (hPa)
- $e_0^*$: saturation vapor pressure at source or sink height (hPa)
- $M$: moisture availability factor (0 to 1)

### 2.3 Constants
- $\rho$: air density
- $c_p$: specific heat of air
- $\gamma$: psychrometric constant
- $\alpha$: Priestley-Taylor type coefficient used in the STIC formulation

---

## 3. Fundamental Penman-Monteith Framework
STIC is posed as a closure of the Penman-Monteith equation:

$$
\lambda E = \frac{s\phi + \rho c_p g_A D_A}{s + \gamma\left(1 + \frac{g_A}{g_C}\right)}
$$

where $\lambda E$ is latent heat flux and $g_A$, $g_C$ are solved analytically from coupled state equations rather than prescribed by empirical parameterizations.

---

## 4. STIC State Equations (2016 Theoretical Core)
Using radiometric temperature information with surface energy balance constraints, STIC uses the following coupled equations:

$$
g_A = \frac{\phi}{\rho c_p \left[(T_0 - T_A) + \frac{(e_0 - e_A)}{\gamma}\right]}
$$

$$
g_C = g_A \frac{(e_0 - e_A)}{(e_0^* - e_0)}
$$

$$
T_0 = T_A + \left(\frac{e_0 - e_A}{\gamma}\right)\left(\frac{1 - \Lambda}{\Lambda}\right)
$$

$$
\Lambda = \frac{2\alpha s}{2s + 2\gamma + \gamma\frac{g_A}{g_C}(1+M)}
$$

These equations provide an internally consistent analytical retrieval of conductances and latent heat flux by embedding thermal-radiometric information directly in the closure.

---

## 5. Moisture Availability and Vapor-State Constraints
The STIC framework defines an effective moisture-availability control $M$ that links thermal state, vapor-pressure gradients, and conductance partitioning.

Conceptually:

1. $M$ regulates the source-air vapor gradient and therefore $g_C$.
2. Under moist conditions, larger effective $M$ supports larger transpiration and weaker canopy-atmosphere coupling constraints.
3. Under dry conditions, reduced $M$ drives stronger stomatal or surface limitation, greater role of atmospheric demand, and stronger coupling behavior.

In the 2016 analysis, this moisture-conductance interaction is central to interpreting wet-season versus dry-season behavior and drought response.

---

## 6. Iterative Radiometric STIC Update (STIC1.2 Behavior Used in 2016 Context)
The HESS 2016 application uses the radiometric STIC logic in which source-air vapor state and conductances are iteratively updated until latent heat flux convergence.

Representative iterative sequence:

1. Initialize $\alpha$ and source vapor states from thermal and meteorological constraints.
2. Compute closure variables ($g_A$, $g_C$, $T_0$ and related vapor terms).
3. Recompute $\lambda E$ from Penman-Monteith form.
4. Update source vapor deficit and moisture state.
5. Repeat until convergence in $\lambda E$.

This iterative feedback allows STIC to represent dynamic transitions between radiation control and moisture or demand control, which is a core feature interpreted in the 2016 study.

---

## 7. ET Partitioning and Conductance Diagnostics

### 7.1 Flux Components
From converged closure states, total latent heat flux is partitioned into:

- transpiration component ($\lambda E_T$)
- evaporation component ($\lambda E_E$)

The 2016 analysis is focused on physically consistent partitioning through conductance-moisture controls rather than empirical vegetation-only scaling.

### 7.2 Biophysical Control Metrics
The primary diagnostics interpreted in the 2016 paper are:

1. $g_C$ response to moisture status and atmospheric demand.
2. $g_A$ modulation of evaporation transfer.
3. Conductance-ratio behavior ($g_C/g_A$) as an indicator of coupling regime.
4. Seasonal and drought hysteresis in the $\lambda E_T$-$g_C$ relationship.

---

## 8. Physical Interpretation Emphasized in Mallick et al. (2016)
The article reports the following behavior patterns captured with STIC:

1. Wet-season transpiration is predominantly radiation controlled in many Amazon sites.
2. Dry-season and drought periods show stronger biophysical control, with increased moisture limitation through $g_C$.
3. Similar canopy-atmosphere coupling can emerge across different plant functional types when soil-moisture stress reshapes $g_C/g_A$.
4. Evaporation remains strongly influenced by aerodynamic transfer, while transpiration shows stronger sensitivity to moisture-driven canopy conductance changes.

---

## 9. Assumptions and Scope
This 2016 specification reflects the article-level scientific formulation:

1. STIC is treated as a physically based single-source closure framework.
2. Conductances are retrieved analytically from thermal-radiometric and meteorological state variables.
3. Interpretation centers on canopy-scale controls and ET partitioning dynamics across hydrologic regimes.

This is intentionally distinct from software-specific runtime options or implementation defaults.

---

## 10. Bibliography
- Mallick, K., Trebs, I., Boegh, E., Giustarini, L., Schlerf, M., Drewry, D. T., et al. (2016). Canopy-scale biophysical controls of transpiration and evaporation in the Amazon Basin. Hydrology and Earth System Sciences, 20, 4237-4264.
- Mallick, K., Jarvis, A., Fisher, J. B., Tu, K. P., Boegh, E., and Niyogi, D. (2014). A Surface Temperature Initiated Closure (STIC) for surface energy balance fluxes. Remote Sensing of Environment, 141, 243-261.
- Mallick, K., Boegh, E., Trebs, I., Alfieri, J. G., Kustas, W. P., Prueger, J. H., et al. (2015). Reintroducing radiometric surface temperature into the Penman-Monteith formulation. Water Resources Research, 51(8), 6214-6243.
