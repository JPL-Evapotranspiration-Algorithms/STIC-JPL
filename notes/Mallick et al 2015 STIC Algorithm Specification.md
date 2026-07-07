# STIC (Surface Temperature Initiated Closure) Algorithm Specification
**Mallick et al. 2015 Configuration (Radiometric STIC 1.2 Formulation)**

## 1. Introduction
This document specifies the STIC framework described in Mallick et al. (2015, Water Resources Research), where radiometric surface temperature ($T_R$) is reintroduced into the Penman-Monteith (PM) formulation to retrieve latent and sensible heat fluxes without exogenous parameterizations of aerodynamic and surface conductances.

The 2015 work extends the 2014 STIC formulation by explicitly introducing:

1. improved moisture-availability retrieval under hysteretic conditions,
2. a moisture-constrained evaporative-fraction state equation,
3. iterative estimation of the Priestley-Taylor coefficient ($\alpha$), and
4. iterative closure updates until stable flux states are obtained.

---

## 2. Core Variables and Constants

### 2.1 Observed and Derived State Variables
- $R_N$: net radiation (W m$^{-2}$)
- $G$: ground heat flux (W m$^{-2}$)
- $\phi = R_N - G$: available energy (W m$^{-2}$)
- $T_R$: radiometric surface temperature (degC)
- $T_A$: air temperature at reference height (degC)
- $RH$: relative humidity
- $e_A$: atmospheric vapor pressure (hPa)
- $D_A = e_A^* - e_A$: atmospheric vapor pressure deficit (hPa)
- $s$: slope of saturation vapor pressure curve at $T_A$

### 2.2 Retrieved STIC State Variables
- $g_B$: aerodynamic conductance (m s$^{-1}$)
- $g_s$: surface conductance (m s$^{-1}$; canopy or aggregated canopy-soil)
- $T_0$: source/sink aerodynamic temperature (degC)
- $e_0$: vapor pressure at source/sink height (hPa)
- $e_0^*$: saturation vapor pressure at source/sink height (hPa)
- $\Lambda$: evaporative fraction
- $M$: moisture availability factor
- $\alpha$: Priestley-Taylor coefficient (retrieved iteratively)

### 2.3 Constants
- $\rho$: air density
- $c_p$: specific heat of air
- $\gamma$: psychrometric constant

---

## 3. Penman-Monteith Basis
STIC solves closure for the PM equation:

$$
\lambda E = \frac{s\phi + \rho c_p g_B D_A}{s + \gamma\left(1 + \frac{g_B}{g_s}\right)}
$$

with unknown internal states $g_B$ and $g_s$ retrieved analytically from coupled state equations using $T_R$, radiation, and meteorology.

---

## 4. STIC State Equations
Using surface energy balance and aerodynamic transfer relations:

$$
\phi = \lambda E + H
$$

$$
H = \rho c_p g_B (T_0 - T_A)
$$

$$
\lambda E = \frac{\rho c_p}{\gamma} g_B (e_0 - e_A) = \frac{\rho c_p}{\gamma} g_s (e_0^* - e_0)
$$

Analytical conductance forms:

$$
g_B = \frac{\phi}{\rho c_p\left[(T_0 - T_A) + \frac{(e_0 - e_A)}{\gamma}\right]}
$$

$$
g_s = g_B\frac{(e_0 - e_A)}{(e_0^* - e_0)}
$$

Source temperature equation:

$$
T_0 = T_A + \left(\frac{e_0 - e_A}{\gamma}\right)\left(\frac{1-\Lambda}{\Lambda}\right)
$$

These equations require additional closure for $\Lambda$ and $e_0$.

---

## 5. Moisture-Constrained Evaporative Fraction ($\Lambda$)
A key 2015 contribution is deriving $\Lambda$ using a moisture-constrained advection-aridity formulation. Final expression:

$$
\Lambda = \frac{2\alpha s}{2s + 2\gamma + \gamma\frac{g_B}{g_s}(1+M)}
$$

This introduces explicit dependence of evaporative fraction on both conductance ratio and moisture availability.

---

## 6. Moisture Availability ($M$), Hysteresis, and Vapor States

### 6.1 Surface Moisture Role
The surface vapor state is represented as:

$$
e_s = e_A(1-M) + M e_s^*
$$

where $M \in [0,1]$ transitions from dry to saturated conditions.

### 6.2 Hysteresis-Aware $M$ Retrieval
Mallick et al. (2015) identifies that single-form $M$ retrieval can overestimate evaporation under dry-down hysteresis between $\lambda E$, $D_A$, and $T_R$.

To account for this, two moisture formulations are used depending on hysteretic regime:

1. a near-surface wetness indicator (2014-style form),
2. a root-zone or hysteresis-sensitive form with explicit $D_A$ feedback.

Representative hysteresis-sensitive form:

$$
M = \frac{(e_s - e_A)\gamma}{(e_s^* - e_s)s + (e_A^* - e_A)\gamma}
$$

and temperature-form expression:

$$
M = \frac{\gamma s_1 (T_{SD} - T_D)}{s_3 (T_R - T_{SD})s + \gamma s_4 (T_A - T_D)}
$$

where $T_{SD}$ is surface dewpoint temperature and $T_D$ is air dewpoint temperature.

### 6.3 Source Vapor Pressure State
The source/sink vapor pressure is expressed as:

$$
e_0 = e_A(1-M) + M e_0^*
$$

This relation links source vapor state to the retrieved moisture condition.

---

## 7. Iterative Retrieval of $\alpha$ and Closure Variables
The 2015 framework derives an analytical $\alpha$ expression from PM decomposition under limiting conditions:

$$
\alpha = \frac{s+\gamma}{s+\gamma\left(1+\frac{g_B}{g_s}\right)} +
\frac{\rho c_p g_B D_A (s+\gamma)}{s\phi\left\{s+\gamma\left(1+\frac{g_B}{g_s}\right)\right\}}
$$

Algorithm sequence:

1. initialize with $\alpha = 1.26$,
2. retrieve initial $M$, $e_0$, $g_B$, $g_s$, $\Lambda$, and $T_0$,
3. update $\alpha$ from retrieved conductance state,
4. recompute closure variables,
5. iterate until stable $\alpha$ and fluxes are obtained.

The paper reports stable solutions typically within approximately 10-12 iterations for $\alpha$ and around 25 iterations for full flux stability in STIC1.2 contexts.

---

## 8. Flux Computation and Partitioning Context
After convergence:

1. compute latent heat flux ($\lambda E$) from PM closure,
2. compute sensible heat flux ($H = \phi - \lambda E$),
3. analyze resulting conductance states and evaporative fraction.

The paper emphasizes physically consistent retrieval of $g_B$ and $g_s$ and improved robustness under dry, advective, and hysteretic conditions, rather than empirical canopy parameter tuning.

---

## 9. Distinction from 2014 STIC
Relative to Mallick et al. (2014), the 2015 framework adds:

1. radiometric iterative feedback with dynamic state updating,
2. moisture-constrained $\Lambda$ formulation,
3. hysteresis-aware moisture retrieval strategy,
4. analytical iterative retrieval of variable $\alpha$,
5. improved behavior under arid and strongly water-limited conditions.

---

## 10. Bibliography
- Mallick, K., Boegh, E., Trebs, I., Alfieri, J. G., Kustas, W. P., Prueger, J. H., et al. (2015). Reintroducing radiometric surface temperature into the Penman-Monteith formulation. Water Resources Research, 51(8), 6214-6243.
- Mallick, K., Jarvis, A., Fisher, J. B., Tu, K. P., Boegh, E., and Niyogi, D. (2014). A Surface Temperature Initiated Closure (STIC) for surface energy balance fluxes. Remote Sensing of Environment, 141, 243-261.
- Brutsaert, W., and Stricker, H. (1979). An advection-aridity approach to estimate actual regional evapotranspiration. Water Resources Research, 15(2), 443-450.
- Priestley, C. H. B., and Taylor, R. J. (1972). On the assessment of surface heat flux and evaporation using large-scale parameters. Monthly Weather Review, 100(2), 81-92.
