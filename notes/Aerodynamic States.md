# Comprehensive Technical Report: Aerodynamic State Dynamics and Iterative Regularization in the STIC-JPL Evapotranspiration Algorithm

---

## Executive Summary

This report documents the technical findings, algorithmic code audit, and physical reasoning surrounding the **Surface Temperature Initiated Closure (STIC)** evapotranspiration model. Specifically, it examines the aerodynamic source vapor pressure variables **`e0`** ($e_0$) and **`e0star`** ($e_0^*$), the **`update_aerodynamic_states`** configuration toggle, the differences between historical **ECOSTRESS Collection 2** and **Collection 3** implementations, and the scientific justification for modifying the algorithm's iterative feedback loop.

---

## 1. Theoretical Definitions of `e0` and `e0star`

In the STIC model, transport phenomena occur across three distinct physical levels:
1. **Reference Height ($z_a$):** Characterized by ambient air temperature ($T_a$) and actual vapor pressure ($e_a$ / `Ea_hPa`).
2. **Physical Skin Surface ($z_s$):** Characterized by radiometric skin surface temperature ($T_s$ / `ST_C`).
3. **Aerodynamic Source Height ($z_0$):** The effective scalar exchange plane inside the vegetation canopy or soil interface.

### `e0` ($e_0$) — Aerodynamic Source Vapor Pressure
`e0` represents the actual vapor pressure of the air at the aerodynamic source height $z_0$ (in $\text{hPa}$). It serves as the intermediate vapor pressure between internal stomatal cavities / soil pores and the free atmosphere.

### `e0star` ($e_0^*$) — Saturation Vapor Pressure at Aerodynamic Temperature
`e0star` represents the saturation vapor pressure corresponding to the aerodynamic surface temperature $T_0$ (in $\text{hPa}$), where $e_0^* = f(T_0)$. Unlike `Estar_hPa` ($e_s^*$), which is evaluated at the radiometric *skin* temperature $T_s$, $e_0^*$ is evaluated at the effective *aerodynamic* temperature $T_0$.

---

## 2. Code Audit: Historical Collection 2 (`STIC.py`) vs. Refactored Pipeline (`model.py`)

A detailed code comparison was conducted between the legacy ECOSTRESS Collection 2 standalone script (`STIC.py`)[cite: 1] and the refactored, multi-configuration pipeline (`model.py`)[cite: 2].

```
+-------------------------------------------------------------------------------+
|                                  STIC Solver                                  |
+-------------------------------------------------------------------------------+
                                        |
       +--------------------------------+--------------------------------+
       |                                                                 |
       v                                                                 v
+---------------------------------------+   +---------------------------------------+
|          Historical STIC.py           |   |            Refactored Pipeline        |
|             (Collection 2)            |   |               (model.py)              |
+---------------------------------------+   +---------------------------------------+
| Hardcoded dynamic refreshes of:       |   | Isolated behind toggle:               |
|  - e0 (Source vapor pressure)[cite: 1]|  |  - update_aerodynamic_states[cite: 2] |
|  - e0star (Sat. vapor pressure)[cite: 1]| |                                       |
|  - phi (Available energy)    [cite: 1]|   | ECOv003 Default: False                |
+---------------------------------------+   +---------------------------------------+
```

### Key Findings from Code Audit

1. **Collection 2 Implementation (`STIC.py`):**
   * During every solver iteration, `STIC.py` calculated source-level vapor pressures `e0star` and `e0` using the latent heat flux estimate `LE_new`[cite: 1].
   * It recalculated soil heat flux $G$ and available energy $\phi = R_n - G$ using `f_G_PHI_actualsurface`[cite: 1].
   * It explicitly passed `e0`, `e0star`, and updated `phi` back into `STIC_closure` on each iteration step[cite: 1]:
     ```python
     [gB, gS, dT, EF] = STIC_closure(delta, phi, e0, Ea_hPa, e0star, M, RHO, CP, GAMMA, alphaN)
     ```

2. **Refactored Implementation (`model.py`):**
   * The refactored `STIC_JPL` function isolated this dynamic update behind the boolean parameter `update_aerodynamic_states`[cite: 2]:
     ```python
     if update_aerodynamic_states:
         Es_hPa = e0
         Estar_hPa = e0star
         phi_Wm2 = Rn_Wm2 - G_Wm2
     ```
   * Initially, `_resolve_mode_defaults("ECOv002")` did not set `"update_aerodynamic_states": True`, causing Collection 2 backward-compatibility runs to default to `False`[cite: 2].

3. **Code Pipeline Resolution:**
   To guarantee exact consistency with historical Collection 2 processing[cite: 1], `_resolve_mode_defaults` was updated to explicitly set `"update_aerodynamic_states": True` for configurations `"ECOv002"`, `"MALLICK2014"`, `"MALLICK2015"`, and `"MALLICK2016"`, while retaining `"update_aerodynamic_states": False` for `"ECOv003"`[cite: 2]. Furthermore, `update_aerodynamic_states` in `STIC_JPL` was converted to `Optional[bool] = None` to dynamically inherit configuration defaults when an explicit override is not passed[cite: 2].

---

## 3. Scientific Justification for Bypassing Updates in Collection 3

In ECOSTRESS Collection 3, disabling `update_aerodynamic_states` (`update_aerodynamic_states = False`) yielded significant improvements in validation accuracy against Eddy Covariance (EC) flux towers[cite: 2].

### Why Dynamic Updates Failed in Satellite Implementations

* **Thermal Noise Exponentiation:** Satellite Land Surface Temperature ($T_s$) carries inherent retrieval uncertainty ($\pm 0.5\text{ to }1.5\text{ K}$). Saturation vapor pressure ($e^*$) depends exponentially on temperature via the Clausius-Clapeyron equation. Updating $e_0^*$ dynamically on every iteration creates a positive feedback loop that amplifies small thermal errors into artificial, non-physical swings in source vapor pressure deficit ($D_0$) and conductances ($g_b, g_s$).
* **Radiometric ($T_s$) vs. Aerodynamic ($T_0$) Disconnect:** STIC's analytical equations assume $T_s \approx T_0$. Over sparse vegetation or dry soils, skin temperature $T_s$ can exceed aerodynamic temperature $T_0$ by $5\text{–}15\text{ K}$. Updating $e_0^*$ using a radiometric skin temperature baseline forces an unadjusted $T_s$ into an aerodynamic framework, leading to severe overestimation of sensible heat flux ($H$) and underestimation of latent heat flux ($LE$).
* **Reanalysis Meteorology Mismatch:** Operational processing ingests coarse gridded reanalysis meteorology ($0.25^\circ\text{–}0.5^\circ$). Coupling coarse atmospheric actual vapor pressure ($e_a$) with high-resolution ($70\text{ m}$) iteratively updated $e_0$ creates unrealistic microclimatic gradients.

### The Operational Solution: Physically-Guided Regularization

Freezing $e_0$ and $e_0^*$ at their initial boundary conditions acts as a numerical damper. It converts STIC from an unconstrained theoretical microclimate solver into an **optimized diagnostic translation engine** that anchors its calculations to the satellite-observed land surface temperature ($T_s$) as an established physical truth.

---

## 4. Summary Matrix of Algorithmic Configurations

| Parameter / Option | Mallick Literature (2014–2016) | ECOSTRESS Collection 2 (`ECOv002`) | ECOSTRESS Collection 3 (`ECOv003`) |
| :--- | :--- | :--- | :--- |
| **`update_aerodynamic_states`** | `True` (Theoretical necessity) | `True` (Hardcoded in legacy `STIC.py`)[cite: 1] | **`False` (Regularized operational default)**[cite: 2] |
| **Closure Version** | 2014 / 2015 | `"2015"`[cite: 2] | `"2015"`[cite: 2] |
| **Dewpoint Method** | Buck (1981) or Linear | Linear (`use_buck_dewpoint = False`)[cite: 1, 2] | Linear (`use_buck_dewpoint = False`)[cite: 2] |
| **Variable $\alpha$** | Enabled | Enabled (`use_variable_alpha = True`)[cite: 1, 2] | Enabled (`use_variable_alpha = True`)[cite: 2] |
| **Grid Stability Focus** | Tower-scale precision | Prone to non-convergence on noisy pixels | **Optimized for global 70m satellite grids** |

---

## 5. Next-Generation Middle-Ground Architectures

While bypassing aerodynamic updates was the correct operational decision for Collection 3, future satellite missions (e.g., NASA SBG, ESA LSTM) can adopt middle-ground regularization strategies:

1. **Under-Relaxed Feedback Loop:** Apply a damping factor $\beta \in (0.15, 0.25)$ to aerodynamic updates:
   $$e_0^{(k+1)} = (1 - \beta) e_0^{(k)} + \beta e_{0, \text{calculated}}^{(k)}$$
2. **Explicit $T_s \to T_0$ Translation:** Integrate excess resistance ($kB^{-1}$) formulations to correct skin temperature $T_s$ to aerodynamic temperature $T_0$ prior to evaluating Clausius-Clapeyron equations.
3. **Regime-Dependent Coupling:** Dynamically enable full iterative updates over dense, low-stress vegetation ($\text{NDVI} > 0.6, T_s - T_a < 2\text{ K}$) while freezing states over sparse or water-stressed pixels.
4. **Bounded Source Deficit ($D_0$) Envelopes:** Enforce strict physical bounds ($e_a \le e_0 \le e_0^*$) during iterations to prevent mathematical divergence.