# Implementing STIC-RCEEP in STIC-JPL

## 1. Purpose
This document defines what needs to be added to this repository to implement STIC-RCEEP, not just STIC ET.

Current repository status:
- STIC ET physics is implemented and returned as energy fluxes (W m^-2).
- ET partition outputs are available (`LE_soil_Wm2`, `LE_canopy_Wm2`, `potential_transpiration_Wm2`).
- The RCEEP GPP coupling layer is not implemented.

Target outcome:
- Add a production-ready STIC-RCEEP pathway that computes daily GPP from STIC-derived ET and ET_T using the Bai et al. (2022) formulation.

## 2. Scope of New Capability
Implement the following model equations from the STIC-RCEEP reference:

- GPP coupling

$$
GPP = (uWUE \cdot \sqrt{b}) \times \sqrt{\frac{1.3 \cdot EVI \cdot ET}{P_a \cdot ET_{max}}} \times ET_T
$$

- Transpiration from ET

$$
ET_T = ET \cdot F_t
$$

- Transpiration fraction

$$
F_t = \max\left(\frac{f_g \cdot f_T \cdot f_c}{f_{SM} + f_c \left(f_g \cdot f_T - f_{SM}\right)}, f_c\right)
$$

- Daily ET upscaling method expected by paper

$$
ET_{daily} = ETrF \cdot ET_{0,daily}, \quad ETrF = \frac{ET_{inst}}{ET_{0,inst}}
$$

Notes:
- STIC ET remains non-parametric.
- STIC-RCEEP introduces calibrated parameter(s), notably $uWUE \cdot \sqrt{b}$.

## 3. High-Level Design
Add a new layer above STIC ET, not inside core STIC closure math.

Recommended architecture:
1. Keep `STIC_JPL/model.py` focused on ET physics.
2. Add new RCEEP module(s) to transform STIC ET outputs into GPP.
3. Expose a clean public API for STIC-RCEEP, with optional calibration utilities.

Suggested new files:
- `STIC_JPL/stic_rceep.py`
- `STIC_JPL/stic_rceep_factors.py`
- `STIC_JPL/stic_rceep_units.py`
- `STIC_JPL/stic_rceep_calibration.py`
- `tests/test_stic_rceep_equations.py`
- `tests/test_stic_rceep_units.py`
- `tests/test_stic_rceep_end_to_end.py`

## 4. Detailed Work Items

### 4.1 API and Data Contract
Implement a top-level function (name can vary, but keep explicit semantics), for example:

`run_STIC_RCEEP(...) -> dict or DataFrame`

Required inputs:
- STIC ET terms: either instantaneous LE plus upscaling inputs, or daily ET directly.
- Meteorology: pressure (`P_a`) or enough state to derive it.
- Vegetation terms: EVI, NDVI/f_c as needed.
- Terms needed to build `f_g`, `f_T`, and `f_SM`.
- Parameters: `uWUE_sqrt_b`, `ET_max` (default per paper), optional crop-type logic.

Required outputs:
- `ET_daily` (if computed)
- `ET_T_daily`
- `F_t`
- `GPP_daily`
- diagnostics for each scalar term in Eq. (3-5)

### 4.2 Unit Harmonization (Critical)
Current STIC outputs are primarily W m^-2 and sometimes kg/m^2/day equivalents via daylight helper.
RCEEP equations require ET and ET_T in mol m^-2 s^-1.

Implement explicit, tested conversion utilities:
1. `LE_Wm2 -> ET_kg_m2_s`
2. `ET_kg_m2_s -> ET_mol_m2_s`
3. Daily aggregation/disaggregation conventions

Define and document constants and assumptions:
- latent heat of vaporization formulation
- molar mass of water
- time basis and daylight vs full-day interpretation

Do not mix units implicitly in RCEEP equations.

### 4.3 ET Daily Upscaling Alignment
Decide one of two implementation tracks:
1. Strict paper parity track:
- Implement explicit `ETrF * ET0_daily` flow (with ET0_inst and ET0_daily).
2. Operational compatibility track:
- Reuse existing daylight upscaling outputs from `daylight_evapotranspiration`, then convert to RCEEP-required units.

Recommendation:
- Implement both behind a mode flag:
  - `upscale_method="etrf"`
  - `upscale_method="daylight_helper"`
- Add tests to quantify differences.

### 4.4 Factor Computation (`f_g`, `f_T`, `f_SM`, `f_c`)
Implement factor computations in one dedicated module with transparent equations and clipping behavior.

Requirements:
- Every factor must have:
  - equation citation in docstring
  - valid range constraints
  - NaN handling policy
- Provide per-factor diagnostics in output to support scientific debugging.

Open issue to resolve before final freeze:
- Confirm exact operational formulas for `f_g`, `f_T`, and `f_SM` to mirror Bai et al. (2021/2022) implementation details.

### 4.5 GPP Coupling Core
Implement Eq. (3-5) in a numerically robust function:
- clip invalid inputs (negative ET, non-physical EVI, non-positive pressure)
- guard all square-root and division operations
- preserve NaN masks from upstream ET retrievals

Return intermediate terms:
- coupling coefficient (`uWUE_sqrt_b`)
- sqrt term in Eq. (3)
- ET_T term

### 4.6 Calibration Workflow
Add optional calibration utility for `uWUE_sqrt_b`:
- objective: minimize RMSE against observed daily GPP
- support global fit and stratified fit (for example C3/C4 crop classes)
- save/load calibrated parameter set (YAML or JSON)

Do not alter STIC ET parameters in this calibration path.

### 4.7 Batch Processing Integration
Extend `process_STIC_table.py` with an opt-in RCEEP mode:
- `run_rceep=True/False`
- include required RCEEP input columns validation
- append output columns (`F_t`, `ET_T`, `GPP` etc.)

Keep backward compatibility for ET-only users.

### 4.8 Public Package Interface
Update exports in `STIC_JPL/__init__.py` and `STIC_JPL/STIC_JPL.py`:
- expose STIC-RCEEP runner
- preserve existing imports and behavior

## 5. Proposed Development Sequence
1. Implement and test unit conversion utilities.
2. Implement factor computations (`f_g`, `f_T`, `f_SM`, `f_c`) with range tests.
3. Implement GPP coupling core Eq. (3-5) with synthetic tests.
4. Integrate with ET upscaling pathways.
5. Add calibration utilities and parameter persistence.
6. Integrate table/batch processing.
7. Add docs and notebooks demonstrating ET-only vs STIC-RCEEP outputs.

## 6. Testing and Validation Requirements

### 6.1 Unit Tests
- Equation-level checks for Eq. (3), Eq. (4), Eq. (5).
- Unit conversion round-trip checks.
- Edge-case tests:
  - zero ET
  - very low/high pressure
  - EVI out of expected range
  - missing factors / NaN propagation

### 6.2 Integration Tests
- End-to-end run from STIC ET outputs to daily GPP.
- Compare `upscale_method="etrf"` vs `upscale_method="daylight_helper"` on sample scenes.

### 6.3 Scientific Regression Tests
Using tower benchmark subsets:
- Compute RMSE, bias, and R^2 for GPP.
- Stratify performance by dryness class.
- Confirm expected dry-condition skill gains relative to baseline LUE options used in this repo.

## 7. Documentation Deliverables
Add/update the following:
- New algorithm note for STIC-RCEEP equations and implementation assumptions.
- README section: ET-only vs STIC-RCEEP capabilities.
- Notebook example showing:
  - STIC ET production
  - factor calculation
  - GPP output
  - calibration workflow

## 8. Risks and Open Decisions
1. Exact operational definitions of `f_g`, `f_T`, `f_SM` must be pinned to a canonical source implementation.
2. Upscaling pathway choice can shift absolute GPP; this must be versioned and documented.
3. Unit mismatches are the highest risk; conversion tests are mandatory before scientific evaluation.
4. Calibration transferability across crop types and climates should be treated as a model configuration, not hard-coded constants.

## 9. Definition of Done
Implementation is complete when all conditions are met:
1. New STIC-RCEEP API returns daily `GPP`, `ET_T`, and `F_t` with diagnostics.
2. Equation and unit tests pass in CI.
3. Batch table processing supports RCEEP mode without breaking ET-only workflows.
4. Calibration utility can fit and persist `uWUE_sqrt_b`.
5. Documentation and one executable notebook demonstrate full workflow.
6. Validation report includes overall and dryness-stratified metrics.
