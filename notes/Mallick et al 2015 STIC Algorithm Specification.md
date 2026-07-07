# STIC (Surface Temperature Initiated Closure) Algorithm Specification
**Mallick et al. 2015 Configuration (Radiometric STIC Update / Current STIC-JPL Lineage)**

## 1. Introduction
This document describes the Mallick et al. (2015) STIC formulation, which reintroduces radiometric surface temperature into the Penman-Monteith framework and is the closest paper-level match to the current STIC-JPL implementation.

The focus is implementation-centric. The equations are written to match the current package logic, while highlighting the major updates relative to the original 2014 STIC formulation: radiometric temperature feedback, iterative closure, variable Priestley-Taylor alpha, and the more complete handling of surface-moisture hysteresis.

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
Current package defaults used by this lineage are:

- `RESAMPLING = "cubic"`
- `DEFAULT_G_METHOD = "santanello"`
- `LE_CONVERGENCE_TARGET_WM2 = 2.0`
- `MAX_ITERATIONS = 30`
- `USE_VARIABLE_ALPHA = True`
- `CONSTRAIN_NEGATIVE_LE = False`
- `CONSTRAIN_PET = False`
- `UPSCALE_TO_DAYLIGHT = False`
- `SHOW_DISTRIBUTIONS = True`
- `APPLY_SURFACE_EMISSIVITY_TO_LWin = False`

These defaults describe the current research/operational package line rather than the strictest interpretation of the original paper.

---

## 3. Inputs and Pathway Decision

### 3.1 Primary Inputs
Core state and forcing inputs include:

- $R_n$: net radiation (W m$^{-2}$)
- $G$: ground heat flux (W m$^{-2}$)
- $T_S$: surface temperature (degC)
- $T_A$: air temperature (degC)
- $RH$: relative humidity (fraction 0-1)
- $\epsilon$: surface emissivity
- $\alpha$: albedo
- $NDVI$: normalized difference vegetation index
- optional $R_g$ via `SWin_Wm2` (incoming shortwave)

### 3.2 Forcing Fallbacks
If `Ta_C` or `RH` are omitted and `offline_mode=False`, they are fetched from GEOS-5 FP.

If `offline_mode=True`, missing `Ta_C` or `RH` raises `MissingOfflineParameter`.

### 3.3 Branching Behavior
The current implementation supports two pathways:

1. **Solar branch** (`SWin_Wm2` provided): the main radiometric STIC path.
2. **No-solar branch** (`SWin_Wm2` missing): fallback path that still updates moisture and closure states.

The 2015-style formulation is most naturally associated with the solar branch and the iterative closure loop, but the package preserves the no-solar branch for incomplete forcing sets.

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

The current code path can optionally use the Buck (1981) dewpoint formulation to preserve the 2014 mode behavior, but the default 2015-style path uses the simpler temperature-RH approximation.

---

## 5. Radiometric STIC Update

### 5.1 Net Longwave Radiation
The solar branch computes atmospheric emissivity and net longwave radiation as:

$$\epsilon_a = 1.24 \left(\frac{e_a}{T_A + 273.15}\right)^{1/7}$$
$$LWin = \sigma \epsilon_a (T_A + 273.15)^4$$
$$LWout = \sigma \epsilon (T_S + 273.15)^4$$

By default,

$$L_{net} = LWin - LWout$$

An optional switch can instead apply surface emissivity to incoming longwave:

$$L_{net} = \epsilon \cdot LWin - LWout$$

### 5.2 Soil Moisture Initialization
The initialization computes:

- $M_{surf}$ from surface thermodynamic constraints
- $M_{rz}$ from root-zone moisture closure
- hysteresis-conditioned composite $M$

The implementation uses conditional logic that responds to:

- fractional vegetation cover thresholds
- vapor pressure comparisons
- thermal and radiative predicates

Representative conditions include $FVC \le 0.25$, $D_{surf} > D_A$, $T_D < 0$, and $LWnet < -125$.

### 5.3 Soil Heat Flux
In active code paths, soil heat flux is computed by a helper rather than by a static closed form:

$$\Phi = R_n - G$$

The default current implementation routes soil heat flux through the Santanello-style helper when available, with SEBAL as the explicit fallback method. This is one of the practical differences between the paper-level description and the package runtime.

### 5.4 Initial Surface Vapor Pressure State
The source vapor pressure is initialized with moisture-weighted mixing:

$$e_s = e_a + M (e_s^* - e_a)$$

This is the same core moisture mixing relation as the 2014 formulation, but in the 2015 lineage it is embedded in a richer iterative state update.

---

## 6. STIC Analytical Closure

The closure solves for $g_B$, $g_S$, $dT$, and $EF$ from the expanded STIC system and then bounds the results physically:

- $g_B \in [0.0001, 0.2]$ m s$^{-1}$
- $g_S \in [0.0001, 0.2]$ m s$^{-1}$
- $dT \in [-10, 50]$ degC
- $EF \in [0, 1]$

The symbolic form remains consistent with Mallick-style STIC closure expressions using $\Phi$, $\Delta$, $\gamma$, $e_s$, $e_a$, $e_s^*$, and $M$.

In operational form, the current implementation uses the 2014-style closure equations as the algebraic core while updating the moisture state and Priestley-Taylor coefficient iteratively.

---

## 7. Iterative Loop and Convergence

### 7.1 Loop Condition
The iterative loop runs while:

$$\max(|LE_{old} - LE_{new}|) \ge LE_{target} \quad \text{and} \quad iteration \le max\_iterations$$

with current defaults:

- $LE_{target} = 2.0$ W m$^{-2}$
- `max_iterations = 30`

### 7.2 Iteration Updates
Each iteration performs:

1. Canopy-air stream update for $e_0^*$, $D_0$, and $e_0$.
2. Soil moisture update with iterative hysteresis logic.
3. Optional variable-$\alpha$ update:

$$\alpha \leftarrow \alpha_N$$

when `use_variable_alpha=True`.

4. Recompute the STIC closure for updated $g_B$, $g_S$, and $dT$.
5. Recompute latent heat flux:

$$LE = \frac{\Delta \Phi + \rho c_p g_B D_A}{\Delta + \gamma\left(1 + \frac{g_B}{g_S}\right)}$$

6. Enforce the energy cap:

$$LE \le \Phi$$

7. If `constrain_negative_LE=True`, clamp negative LE to zero.

### 7.3 Dryness Adjustment
The current implementation includes additional dryness logic that becomes important in arid scenes and strong advection cases. The logic is driven by combinations of:

- source versus atmospheric deficit
- potential evaporation versus available energy
- thermal indicators such as $dT_S > 0$ and $T_D \le 0$

This is the main place where the 2015 implementation extends the simpler 2014 logic.

---

## 8. Final Flux Partitioning and Diagnostics

After convergence or max-iteration stop, outputs include:

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

## 9. Optional Daily Upscaling

Daily upscaling is available but disabled by default (`UPSCALE_TO_DAYLIGHT=False`).

When enabled and `time_UTC` is present, the model calls daylight ET upscaling and merges the resulting daily products into the output dictionary.

---

## 10. Collection 3 / Implementation Notes

This 2015-lineage implementation differs from the original 2014 paper in several practical ways:

1. It uses an iterative convergence loop by default.
2. Variable $\alpha$ is enabled by default.
3. The soil-moisture state includes hysteresis logic and root-zone updates.
4. The package supports both solar and no-solar pathways.
5. Optional physical constraints on LE and PET remain available but are off by default.
6. Daily upscaling is optional and off by default.

These settings make the implementation flexible for research and operational workflows while preserving the STIC closure physics introduced in the Mallick et al. line of work.

---

## 11. Bibliography

- Mallick, K., Boegh, E., Trebs, I., Alfieri, J. G., Kustas, W. P., Prueger, J. H., ... & Jarvis, A. J. (2015). Reintroducing radiometric surface temperature into the Penman-Monteith formulation. Water Resources Research, 51(8), 6214-6243.
- Mallick, K., Toivonen, E., Trebs, I., Boegh, E., Cleverly, J., Eamus, D., ... & Garcia, M. (2018). Bridging Thermal Infrared Sensing and Physically-Based Evapotranspiration Modeling: From Theoretical Implementation to Validation Across an Aridity Gradient in Australian Ecosystems. Water Resources Research, 54(5), 3409-3435.
- Mallick, K., Baldocchi, D., Jarvis, A., Hu, T., Trebs, I., Sulis, M., ... & Kustas, W. P. (2022). Insights Into the Aerodynamic Versus Radiometric Surface Temperature Debate in Thermal-Based Evaporation Modeling. Geophysical Research Letters, 49(15), e2021GL097568.
