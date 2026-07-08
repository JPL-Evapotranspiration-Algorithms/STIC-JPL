# STIC (Surface Temperature Initiated Closure) Algorithm Specification
**ECOSTRESS Collection 3 Configuration (Current STIC-JPL Implementation)**

## 1. Introduction
This document describes the current STIC-JPL configuration used for Collection 3-oriented processing in this repository.

It is intentionally implementation-centric: equations, defaults, branching logic, and constraints are documented as they are currently encoded in the Python package, rather than as a strict historical C2 reconstruction.

---

## 2. Global Constants and Runtime Defaults

### 2.1 Physical Constants
The current implementation uses:

- $\sigma = 5.67 \times 10^{-8}$ W m$^{-2}$ K$^{-4}$ (Stefan-Boltzmann constant)
- $\rho = 1.2$ kg m$^{-3}$ (air density)
- $c_p = 1013$ J kg$^{-1}$ K$^{-1}$ (specific heat of air)
- $\gamma = 0.67$ hPa K$^{-1}$ (psychrometric constant)
- $\alpha_{PT} = 1.26$ (Priestley-Taylor coefficient)

### 2.2 Model Configuration Defaults
Current package defaults are:

- `RESAMPLING = "cubic"`
- `DEFAULT_G_METHOD = "santanello"`
- `LE_CONVERGENCE_TARGET_WM2 = 2.0`
- `MAX_ITERATIONS = 30`
- `USE_VARIABLE_ALPHA = True`
- `CONSTRAIN_NEGATIVE_LE = False`
- `CONSTRAIN_PET = False`
- `UPSCALE_TO_DAYLIGHT = False`
- `SHOW_DISTRIBUTIONS = True`

These defaults create a more permissive and longer-iteration configuration than Collection 2 operational profiles that often emphasized tighter iteration limits.

---

## 3. Inputs and Pathway Decision

### 3.1 Primary Inputs
Core state and forcing inputs include:

- $R_n$: net radiation (W m$^{-2}$)
- $T_A$: air temperature (degC)
- $T_S$: surface temperature (degC)
- $RH$: relative humidity (fraction 0-1)
- $\epsilon$: surface emissivity
- $\alpha$: albedo
- $NDVI$: normalized difference vegetation index
- optional $R_g$ via `SWin_Wm2` (incoming shortwave)

### 3.2 Forcing Fallbacks
If `Ta_C` or `RH` are omitted and `offline_mode=False`, they are fetched from GEOS-5 FP.

If `offline_mode=True`, missing `Ta_C` or `RH` raises `MissingOfflineParameter`.

### 3.3 Branching Behavior
The model still supports two pathways:

1. **Solar branch** (`SWin_Wm2` provided): STIC 1.2/1.3-style iterative branch.
2. **No-solar branch** (`SWin_Wm2` missing): fallback legacy-style initialization and iteration.

So, unlike strict C2 pipeline assumptions where SWin was always provided, the current implementation remains dual-path.

---

## 4. Psychrometric Initialization

The current implementation computes:

$$e_{s,A} = 6.13753 \cdot \exp\left(\frac{17.27 T_A}{T_A + 237.3}\right)$$
$$e_a = e_{s,A} \cdot RH$$
$$D_A = e_{s,A} - e_a$$
$$\Delta = \frac{4098 e_{s,A}}{(T_A + 237.3)^2}$$
$$T_D = T_A - \frac{100 - 100 RH}{5.0}$$
$$dT_S = T_S - T_A$$
$$e_s^* = 6.13753 \cdot \exp\left(\frac{17.27 T_S}{T_S + 237.3}\right)$$

If not provided, $\Delta$ is derived from $T_A$ and $e_{s,A}$ as above.

---

## 5. Solar Branch (Collection 3 Current Configuration)

When `SWin_Wm2` is provided, initialization uses:

### 5.1 Net Longwave
Atmospheric emissivity and net longwave are computed as:

$$\epsilon_a = 1.24 \left(\frac{e_a}{T_A + 273.15}\right)^{1/7}$$
$$LWin = \sigma \epsilon_a (T_A + 273.15)^4$$
$$LWout = \sigma \epsilon (T_S + 273.15)^4$$

Default behavior uses:

$$L_{net} = LWin - LWout$$

An optional switch can apply surface emissivity to incoming longwave:

$$L_{net} = \epsilon \cdot LWin - LWout$$

### 5.2 Soil Moisture Initialization
The initialization computes:

- $M_{surf}$ from surface thermodynamic constraints
- $M_{rz}$ from root-zone moisture closure
- hysteresis-conditioned composite $M$

The implementation includes aridity and low-cover conditional logic using:

- fractional vegetation cover thresholds (`FVC <= 0.25`)
- vapor pressure comparisons (e.g., $D_{surf} > D_A$)
- thermal and radiative predicates (e.g., $T_D < 0$, $LWnet < -125$)

### 5.3 Soil Heat Flux in Current Code
Although the interface exposes `G_method`, the current code path computes soil heat flux through `calculate_SEBAL_soil_heat_flux` in both initialization and iteration. Available energy is:

$$\Phi = R_n - G$$

### 5.4 Initial Surface Vapor Pressure State
The source vapor pressure term is initialized with conditional root-zone versus surface wetness selection:

$$e_s = e_a + M( e_s^* - e_a )$$

with $M$ selected by branch conditions tied to dryness/advection logic.

---

## 6. No-Solar Branch (Fallback Path)

When `SWin_Wm2` is missing:

1. If `G_Wm2` is missing, $G$ is estimated using SEBAL soil heat flux.
2. $\Phi = R_n - G$.
3. Moisture state is initialized without explicit incoming shortwave forcing.

The no-solar iterative update still computes canopy-air stream vapor states, updates moisture hysteresis, and updates $\alpha_N$ when variable alpha is enabled.

---

## 7. STIC Analytical Closure (Current)

The closure solves for $g_B$, $g_S$, and $dT$ from the expanded STIC algebraic system, then applies hard bounds:

- $g_B \in [0.0001, 0.2]$ m s$^{-1}$
- $g_S \in [0.0001, 0.2]$ m s$^{-1}$
- $dT \in [-10, 50]$ degC
- $EF \in [0, 1]$

The symbolic form remains consistent with Mallick-style STIC closure expressions using $\Phi$, $\Delta$, $\gamma$, $e_s$, $e_a$, $e_s^*$, and $M$.

---

## 8. Iterative Loop and Convergence (Current Defaults)

## 8.1 Loop Condition
The loop iterates while:

$$\max(|LE_{old} - LE_{new}|) \ge LE_{target} \quad \text{and} \quad iteration \le max\_iterations$$

with current defaults:

- $LE_{target} = 2.0$ W m$^{-2}$
- `max_iterations = 30`

## 8.2 Iteration Updates
Each iteration performs:

1. Canopy-air stream update for $e_0^*$, $D_0$, and $e_0$.
2. Soil moisture update with iterative hysteresis logic.
3. Optional variable-$\alpha$ update:

$$\alpha \leftarrow \alpha_N$$

when `use_variable_alpha=True`.

4. Recompute STIC closure for updated $g_B$, $g_S$, and $dT$.
5. Recompute latent heat flux:

$$LE = \frac{\Delta \Phi + \rho c_p g_B D_A}{\Delta + \gamma\left(1 + \frac{g_B}{g_S}\right)}$$

6. Enforce energy cap:

$$LE \le \Phi$$

7. If `constrain_negative_LE=True`, clamp negative LE to zero.

## 8.3 Dryness Adjustment
The implementation includes additional dryness logic for initialization and stability in arid conditions, using combinations of:

- source versus atmospheric deficit
- potential evaporation versus available energy
- thermal indicators (e.g., $dT_S > 0$, $T_D \le 0$)

---

## 9. Final Flux Partitioning and Diagnostics

After convergence (or max-iteration stop), outputs include:

- `LE_Wm2`
- `LE_change`
- `LE_soil_Wm2`
- `LE_canopy_Wm2`
- `potential_transpiration_Wm2`
- `PET_Wm2`
- `G_Wm2`
- `iteration`
- `LE_max_change`

Partitioning follows:

$$LE_{soil} = \text{clip}(M \cdot PET, 0, \infty)$$
$$LE_{canopy} = \text{clip}(LE - LE_{soil}, 0, \infty)$$

Potential evaporation is Penman-type:

$$PET = \frac{\Delta \Phi + \rho c_p g_B D_A}{\Delta + \gamma}$$

Potential transpiration follows the STIC moisture-weighted conductance form.

---

## 10. Optional Daily Upscaling

Daily upscaling is available but disabled by default (`UPSCALE_TO_DAYLIGHT=False`).

When enabled and `time_UTC` is present, the model calls daylight ET upscaling and merges daily products into the result dictionary.

---

## 11. Collection 3 Operational Notes for This Codebase

This current implementation differs from a strict C2 profile in several practical ways:

1. Longer and looser convergence defaults (`30` iterations, `2.0` W m$^{-2}$ target).
2. Dual-path support (with and without SWin) remains active.
3. PET and negative LE constraints are configurable but default off.
4. Daily upscaling is optional and default off.
5. Soil heat flux is currently computed via the SEBAL helper in active code paths.

These settings make the implementation flexible for research and cross-sensor workflows while retaining STIC 1.2/1.3 closure and iterative physics.

---

## 12. Bibliography

- Mallick, K., et al. (2014). A Surface Temperature Initiated Closure (STIC) for surface energy balance fluxes. Remote Sensing of Environment, 141, 243-261.
- Mallick, K., et al. (2015). Reintroducing radiometric surface temperature into the Penman-Monteith formulation. Water Resources Research, 51(8), 6214-6243.
- Mallick, K., et al. (2018). Bridging Thermal Infrared Sensing and Physically-Based Evapotranspiration Modeling across an aridity gradient. Water Resources Research, 54(5), 3409-3443.
- Mallick, K., et al. (2022). Insights into aerodynamic versus radiometric surface temperature in thermal-based evaporation modeling. Geophysical Research Letters, 49(15), e2021GL097568.
