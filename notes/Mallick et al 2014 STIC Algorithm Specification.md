# STIC (Surface Temperature Initiated Closure) Algorithm Specification
**Mallick et al. 2014 Configuration (Original STIC Formulation)**

## 1. Introduction
This document specifies the original STIC formulation presented in Mallick et al. (2014, Remote Sensing of Environment), where radiometric surface temperature is physically integrated into the Penman-Monteith framework to retrieve latent and sensible heat fluxes.

The central objective of the 2014 formulation is to avoid exogenous parameterization of aerodynamic and surface conductances by treating them as internal state variables recovered through analytical closure.

---

## 2. Inputs, States, and Constants

### 2.1 Primary Inputs
- $R_N$: net radiation (W m$^{-2}$)
- $G$: ground heat flux (W m$^{-2}$)
- $\Phi = R_N - G$: available energy (W m$^{-2}$)
- $T_S$: radiometric surface temperature (degC)
- $T_A$: air temperature (degC)
- $RH$ (or $e_A$): relative humidity (or atmospheric vapor pressure)

### 2.2 Retrieved Internal States
- $g_B$: aerodynamic (boundary layer) conductance
- $g_S$: surface conductance (canopy for full vegetation; aggregated canopy-soil under partial cover)
- $\Delta T$: aerodynamic temperature difference between source/sink and air
- $\Lambda$: evaporative fraction
- $M$: near-surface moisture availability
- $e_S$: effective vapor pressure at the evaporating front

### 2.3 Constants
- $\rho$: air density
- $c_P$: specific heat of air
- $\gamma$: psychrometric constant
- $s$: slope of saturation vapor pressure curve at $T_A$
- $\alpha$: Priestley-Taylor coefficient (nominally 1.26 for wet, weak-advection conditions)

---

## 3. Penman-Monteith Basis
The STIC closure is built on the Penman-Monteith equation:

$$
\lambda E = \frac{s\Phi + \rho c_P g_B D_A}{s + \gamma\left(1 + \frac{g_B}{g_S}\right)}
$$

where $D_A$ is atmospheric vapor pressure deficit.

The original PM equation contains unknown $g_B$ and $g_S$. STIC derives these analytically using coupled state equations constrained by $T_S$, meteorology, and radiation.

---

## 4. Moisture Availability and Effective Surface Vapor Pressure
STIC introduces near-surface moisture availability $M$ and defines effective evaporating-front vapor pressure:

$$
e_S = e_A(1-M) + M e_S^*
$$

with $M \in [0,1]$, where:

- $M \to 0$: dry surface, minimal latent heat transfer,
- $M \to 1$: saturated evaporating front.

Equivalent vapor-gradient form:

$$
M = \frac{e_S - e_A}{e_S^* - e_A}
$$

### 4.1 Dewpoint-Based Retrieval Logic
The 2014 paper derives $M$ through a notional surface dewpoint temperature ($T_{SD}$), linking psychrometric gradients between:

- air dewpoint and surface dewpoint,
- surface dewpoint and radiometric surface temperature.

This yields a physically constrained estimate of $M$ from $T_S$, $T_A$, and humidity information, using local linearization of the saturation vapor pressure curve.

---

## 5. Evaporative Fraction from Advection-Aridity Theory
To close the system, STIC derives an explicit conductance-dependent expression for evaporative fraction $\Lambda$ by combining:

1. Penman potential evaporation,
2. Priestley-Taylor potential evaporation,
3. Brutsaert-Stricker advection-aridity complementarity.

Final 2014 expression:

$$
\Lambda = \frac{2\alpha s}{2s + \gamma\left(2 + \frac{g_B}{g_S}\right)}
$$

This relation is a key component enabling analytical closure of the PM framework.

---

## 6. Original STIC Closure Equations (Four Unknowns, Four Equations)
The 2014 STIC system solves for $g_B$, $g_S$, $\Delta T$, and $\Lambda$ with:

$$
g_B = \frac{\Phi}{\rho c_P\left(\Delta T + \frac{e_S - e_A}{\gamma}\right)}
$$

$$
g_S = g_B\frac{e_S - e_A}{e_S^* - e_S}
$$

$$
\Delta T = \left(\frac{e_S - e_A}{\gamma}\right)\left(\frac{1-\Lambda}{\Lambda}\right)
$$

$$
\Lambda = \frac{2\alpha s}{2s + \gamma\left(2 + \frac{g_B}{g_S}\right)}
$$

Solving this system provides analytical retrievals of the unobserved conductance and thermodynamic state variables.

---

## 7. Flux Computation
With closure variables recovered:

1. latent heat flux from PM equation,
2. sensible heat flux from residual energy balance:

$$
H = \Phi - \lambda E
$$

These are the core flux outputs evaluated against eddy-covariance measurements in the 2014 study.

---

## 8. Partitioning of Latent Heat Components
The 2014 framework uses $M$ to partition latent heat into evaporation and transpiration components:

$$
\lambda E = \lambda E_E + \lambda E_T = M\lambda E^* + (1-M)\lambda E_T
$$

where:

- $\lambda E^*$ is potential evaporation,
- $\lambda E_E$ is evaporation-dominated component,
- $\lambda E_T$ is transpiration component.

This partitioning enables interpretation of moisture controls on canopy versus soil contributions.

---

## 9. Practical Scope of the 2014 Formulation
The paper emphasizes:

1. analytical closure without empirical conductance submodels,
2. explicit use of radiometric temperature as a state constraint,
3. physically based retrieval of moisture availability,
4. broad biome applicability with strong relevance to water-limited and semi-arid regimes.

The 2014 framework is a single-pass analytical closure formulation; later iterative enhancements are introduced in subsequent STIC versions.

---

## 10. Bibliography
- Mallick, K., Jarvis, A., Fisher, J. B., Tu, K. P., Boegh, E., and Niyogi, D. (2014). A Surface Temperature Initiated Closure (STIC) for surface energy balance fluxes. Remote Sensing of Environment, 141, 243-261.
- Brutsaert, W., and Stricker, H. (1979). An advection-aridity approach to estimate actual regional evapotranspiration. Water Resources Research, 15(2), 443-450.
- Priestley, C. H. B., and Taylor, R. J. (1972). On the assessment of surface heat flux and evaporation using large-scale parameters. Monthly Weather Review, 100(2), 81-92.
- Penman, H. L. (1948). Natural evaporation from open water, bare soil and grass. Proceedings of the Royal Society A, 193(1032), 120-145.
- Monteith, J. L. (1965). Evaporation and environment. In The State and Movement of Water in Living Organisms (19th Symposia of the Society for Experimental Biology), 205-234.
