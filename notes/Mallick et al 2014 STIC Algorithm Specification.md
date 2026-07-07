# STIC (Surface Temperature Initiated Closure) Algorithm Specification
**Mallick et al. 2014 Configuration (Original STIC Formulation)**

## 1. Introduction
This document describes the original STIC formulation introduced by Mallick et al. (2014) and the corresponding STIC-JPL package mode used to represent that paper in this repository.

The emphasis is implementation-facing rather than purely historical: the equations are written in the same style as the current code path, while noting where the original 2014 formulation differs from the later radiometric-temperature update and from the current default runtime configuration.

---

## 2. Global Constants and Runtime Defaults

### 2.1 Physical Constants
The current implementation uses the following constants for this mode:

- $\sigma = 5.67 \times 10^{-8}$ W m$^{-2}$ K$^{-4}$ (Stefan-Boltzmann constant)
- $\rho = 1.2$ kg m$^{-3}$ (air density)
- $c_p = 1013$ J kg$^{-1}$ K$^{-1}$ (specific heat of air)
- $\gamma = 0.67$ hPa K$^{-1}$ (psychrometric constant)
- $\alpha_{PT} = 1.26$ (Priestley-Taylor coefficient)

### 2.2 Package Mode Mapping
The repository exposes a dedicated `MALLICK2014` mode with the following behavior:

- `G_method = "sebal"`
- `RUN_ITERATIVE_CONVERGENCE = False`
- `USE_BUCK_DEWPOINT = True`
- `USE_VARIABLE_ALPHA = False`
- `CONSTRAIN_NEGATIVE_LE = False`
- `CONSTRAIN_PET = False`
- `CONSTRAIN_LE_TO_AVAILABLE_ENERGY = False`
- `APPLY_SURFACE_EMISSIVITY_TO_LWin = False`

This mode is intentionally conservative and closer to the original paper’s analytical presentation than to the more iterative current default configuration.

---

## 3. Inputs and Pathway Decision

### 3.1 Primary Inputs
The 2014 STIC formulation uses the following core inputs:

- $R_n$: net radiation (W m$^{-2}$)
- $G$: ground heat flux (W m$^{-2}$)
- $T_S$: surface temperature (degC)
- $T_A$: air temperature (degC)
- $RH$: relative humidity (fraction 0-1)
- $\epsilon$: surface emissivity
- $\alpha$: albedo
- $NDVI$: normalized difference vegetation index

The current package also accepts optional incoming shortwave radiation through `SWin_Wm2`, but the 2014 mode does not require it.

### 3.2 Forcing Fallbacks
If `Ta_C` or `RH` are omitted and `offline_mode=False`, they are fetched from GEOS-5 FP.

If `offline_mode=True`, missing `Ta_C` or `RH` raises `MissingOfflineParameter`.

### 3.3 Branching Behavior
The package still supports both pathways:

1. **With-solar branch** when `SWin_Wm2` is provided.
2. **No-solar branch** when `SWin_Wm2` is missing.

For the 2014 mode, the no-solar branch is the closest representation of the original paper’s direct STIC closure logic, because the paper itself does not depend on a later iterative radiometric-temperature update.

---

## 4. Psychrometric Initialization

The standard initialization computes:

$$e_{s,A} = 6.13753 \cdot \exp\left(\frac{17.27 T_A}{T_A + 237.3}\right)$$
$$e_a = e_{s,A} \cdot RH$$
$$D_A = e_{s,A} - e_a$$
$$\Delta = \frac{4098 e_{s,A}}{(T_A + 237.3)^2}$$
$$T_D = T_A - \frac{100 - 100 RH}{5.0}$$
$$dT_S = T_S - T_A$$
$$e_s^* = 6.13753 \cdot \exp\left(\frac{17.27 T_S}{T_S + 237.3}\right)$$

The 2014 package mode uses the Buck-style dewpoint formulation for consistency with the original-mode handling in this repository.

---

## 5. Original STIC Formulation

### 5.1 Moisture Availability
The original paper defines surface moisture availability through the effective evaporating-front vapor pressure:

$$e_s = e_a(1 - M) + M e_s^*$$

which is equivalent to:

$$M = \frac{e_s - e_a}{e_s^* - e_a}$$

The paper derives $M$ from the surface dewpoint construction $T_{SD}$ and the slopes of the saturation vapor pressure curve. In practice, the current package mode uses the same wetness logic and hysteresis structure that leads to $M$ and root-zone wetness states $M_{rz}$.

### 5.2 STIC Closure Equations
The 2014 analytical closure solves for $g_B$, $g_S$, $dT$, and $EF$ from the following system:

$$g_B = \frac{\Phi}{\rho c_p \left(\Delta T + \frac{e_s - e_a}{\gamma}\right)}$$
$$g_S = g_B \frac{(e_s - e_a)}{(e_s^* - e_s)}$$
$$\Delta T = \left(\frac{e_s - e_a}{\gamma}\right)\left(\frac{1 - \Lambda}{\Lambda}\right)$$
$$\Lambda = \frac{2\alpha s}{2s + \gamma \left(2 + \frac{g_B}{g_S}\right)}$$

where $\Phi = R_n - G$.

The closure is then used to compute latent heat flux from the Penman-Monteith form:

$$LE = \frac{\Delta \Phi + \rho c_p g_B D_A}{\Delta + \gamma\left(1 + \frac{g_B}{g_S}\right)}$$

and the corresponding sensible heat flux is:

$$H = \Phi - LE$$

### 5.3 Flux Partitioning
The 2014 formulation partitions latent heat into soil and canopy components using moisture availability:

$$LE_{soil} = \text{clip}(M \cdot PET, 0, \infty)$$
$$LE_{canopy} = \text{clip}(LE - LE_{soil}, 0, \infty)$$

Potential evaporation is computed in Penman form:

$$PET = \frac{\Delta \Phi + \rho c_p g_B D_A}{\Delta + \gamma}$$

Potential transpiration is computed from the same conductance framework but weighted by the moisture state.

---

## 6. Soil Heat Flux and Available Energy

In the current package mode, soil heat flux is computed with the SEBAL helper:

$$\Phi = R_n - G$$

The `MALLICK2014` mode does not route the active closure through the later Santanello-based iterative flux path.

---

## 7. Iteration and Convergence Behavior

### 7.1 Loop Condition
The original 2014 mode is exposed as a non-iterative mode in the current package:

- `run_iterative_convergence = False`
- `max_iterations = 0`

This matches the paper’s emphasis on a direct analytical closure rather than a long iterative refinement loop.

### 7.2 Operational Behavior in the Package
When the 2014 mode is used, the package computes the initial closure and flux partitioning directly, without the later variable-$\alpha$ update loop used in the radiometric-temperature update path.

---

## 8. Output Products
The current package returns the standard STIC outputs for this mode:

- `LE_Wm2`
- `LE_change`
- `LE_soil_Wm2`
- `LE_canopy_Wm2`
- `potential_transpiration_Wm2`
- `PET_Wm2`
- `G_Wm2`
- `iteration`
- `LE_max_change`

For the original 2014 configuration, the emphasis is on the closure variables and moisture partitioning rather than on a lengthy convergence history.

---

## 9. Collection 3 Comparison Notes
Relative to the later STIC-JPL default runtime, the 2014 mode differs mainly in implementation posture rather than in the core physics:

1. It is non-iterative in the current package.
2. It does not enable variable $\alpha$ updates.
3. It uses Buck dewpoint handling.
4. It does not impose additional LE or PET constraints by default.
5. It keeps the analytical STIC closure as the central algorithmic object.

---

## 10. Bibliography

- Mallick, K., Jarvis, A., Fisher, J. B., Tu, K. P., Boegh, E., & Niyogi, D. (2014). A Surface Temperature Initiated Closure (STIC) for surface energy balance fluxes. Remote Sensing of Environment, 141, 243-261.
- Mallick, K., Boegh, E., Trebs, I., Alfieri, J. G., Kustas, W. P., Prueger, J. H., ... & Jarvis, A. J. (2015). Reintroducing radiometric surface temperature into the Penman-Monteith formulation. Water Resources Research, 51(8), 6214-6243.
