# ECOSTRESS STIC-JPL: Collection 2 vs Collection 3 Differences

## 1. Scope
This document explains the differences between:

- **Collection 2 (C2):** the ECOSTRESS Collection 2 algorithm specification in this repository's notes.
- **Collection 3 (C3):** the current STIC-JPL implementation and defaults in code.

It focuses on practical configuration and behavior differences, not just theoretical equation forms.

---

## 2. At-a-Glance Summary

1. C3 uses **looser convergence defaults** (more iterations, higher LE threshold).
2. C3 keeps **dual-path operation** (with and without SWin), while C2 operations are centered on SWin-driven STIC 1.2/1.3.
3. C3 currently computes soil heat flux through the **SEBAL helper in active paths**, even where Santanello-style framing is referenced.
4. C3 has several **optional physical constraints disabled by default** (negative LE and PET constraints).
5. C3 introduces a few **implementation toggles and convenience behaviors** (GEOS forcing fallback, optional daylight upscaling, longwave emissivity option).

---

## 3. Detailed Differences

## 3.1 Convergence Target and Iteration Limit
**C2:**
- Iterative loop target is tighter: stop when $\max(LE_{change}) < 1.0$ W m$^{-2}$.
- Maximum loop count is 3 iterations.

**C3 (current code):**
- Default convergence target is $2.0$ W m$^{-2}$.
- Default `max_iterations` is 30.

**Practical effect:**
- C3 can run significantly longer and accept a larger final LE residual by default.
- C3 is more tolerant for challenging scenes, but less aligned with strict C2 operational stopping rules.

---

## 3.2 Pathway Selection (SWin vs No-SWin)
**C2:**
- Operational framing emphasizes the SWin ($R_g$) branch; JET pipeline passes SWin and runs STIC 1.2/1.3 loop.

**C3 (current code):**
- Still supports two pathways:
1. with `SWin_Wm2` (solar branch)
2. without `SWin_Wm2` (fallback/legacy-style path)

**Practical effect:**
- C3 remains broadly usable for datasets where incoming shortwave is missing.
- Behavior is less pipeline-fixed than C2 operations.

---

## 3.3 Soil Heat Flux Implementation
**C2 spec narrative:**
- Describes Santanello and Friedl style diurnal $G$ treatment in the solar branch.

**C3 (current code):**
- Active initialization and iteration paths call `calculate_SEBAL_soil_heat_flux`.
- `G_method` appears in interfaces but is not currently driving alternate computational branches in the main flow.

**Practical effect:**
- Current C3 runtime behavior for $G$ is effectively tied to the SEBAL helper in the active code paths.
- If strict Santanello behavior is required, code-level method routing should be hardened.

---

## 3.4 Longwave Radiation Handling Toggle
**C2:**
- Spec equation form applies surface emissivity to incoming longwave term.

**C3 (current code):**
- Net longwave function supports two forms:
1. default: do not apply emissivity multiplier to incoming longwave
2. optional: apply emissivity multiplier (`apply_surface_emissivity_to_LWin=True`)
- Solar initialization currently uses the default (option off).

**Practical effect:**
- C3 defaults may produce slightly different longwave partitioning than strict C2 equation form, depending on emissivity and temperature regime.

---

## 3.5 Physical Constraints as Defaults
**C2 intent:**
- The algorithmic narrative is strongly physically constrained, especially in arid-regime handling and bounded states.

**C3 (current defaults):**
- `CONSTRAIN_NEGATIVE_LE = False`
- `CONSTRAIN_PET = False`
- Constraint logic exists but must be enabled.

**Practical effect:**
- C3 defaults are more permissive and may allow physically marginal PET/LE outcomes unless constraints are explicitly activated.

---

## 3.6 Forcing and Data-Access Behavior
**C2 operations:**
- In production pipelines, meteorological forcing is generally explicit and controlled.

**C3 (current code):**
- If `Ta_C` or `RH` are not supplied and `offline_mode=False`, values are fetched from GEOS-5 FP.
- `offline_mode=True` enforces user-provided `Ta_C` and `RH`.
- Default resampling is `"cubic"`.

**Practical effect:**
- C3 supports both operational and research modes more flexibly.
- Reproducibility depends on explicitly pinning forcing and mode settings.

---

## 3.7 Daily Upscaling Behavior
**C2 context:**
- Product framing is strongly tied to ET production pipelines.

**C3 (current code):**
- Daylight upscaling is available, but default is off (`UPSCALE_TO_DAYLIGHT=False`).

**Practical effect:**
- C3 instantaneous output is the default; daily products require explicit opt-in.

---

## 3.8 Diagnostics and Runtime Instrumentation
**C2:**
- Specification is algorithmic and product-oriented.

**C3 (current code):**
- Includes runtime diagnostics and distribution checks (`check_distribution`) and richer intermediate outputs.

**Practical effect:**
- C3 is more analysis/debug friendly for development and calibration workflows.

---

## 4. Areas That Are Largely Unchanged

1. Core psychrometric setup equations are consistent (SVP, VPD, dewpoint formulation).
2. STIC closure structure for $g_B$, $g_S$, and $dT$ is retained.
3. Closure bounds on conductances and aerodynamic temperature are retained.
4. Iterative $\alpha$ update remains central (variable alpha enabled by default in C3).

---

## 5. Recommended “C2-Like” Settings in C3 Code
To approximate C2 behavior more closely when using current code:

1. Set `max_iterations=3`.
2. Set `LE_convergence_target=1.0`.
3. Ensure `SWin_Wm2` is provided for all scenes.
4. Decide and pin longwave emissivity handling (`apply_surface_emissivity_to_LWin`) if strict equation matching is required.
5. Consider enabling `constrain_negative_LE=True` and `constrain_PET=True` for stricter physical control.
6. Pin forcing inputs (`Ta_C`, `RH`) and run with `offline_mode=True` for reproducible studies.

---

## 6. Bottom Line
Collection 3, as currently implemented, keeps the same STIC scientific backbone but is configured as a **more flexible, research-capable runtime** than a strict Collection 2 operations profile. Most differences are in defaults, method routing, and optional constraints rather than a wholesale change in core closure physics.
