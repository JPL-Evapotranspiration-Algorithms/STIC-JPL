# Nuances Between the 2014 and 2015 Versions of STIC

The Surface Temperature Initiated Closure (STIC) model physically integrates radiometric land surface temperature into the Penman-Monteith (PM) formulation to estimate terrestrial sensible and latent heat fluxes ($H$ and $\lambda E$) without requiring exogenous conductance parameterizations. However, the 2015 version introduces vital structural modifications to correct the 2014 version's tendency to overestimate evapotranspiration ($E$) under extremely dry land surface conditions.

---

## Key Architectural Differences

### 1. Priestley-Taylor Parameter (alpha) Estimation
* **2014 Version:** Treated the Priestley-Taylor parameter (alpha) as a fixed, constant exogenous input.
* **2015 Version:** Introduces an analytical, physical expression decomposed from the PM equation to calculate alpha under limited environmental and ecohydrological conditions.
* **Iteration Process:** The 2015 version dynamically updates alpha through numerical iterations starting from an initial value of 1.26 until a stable value is achieved.

### 2. Evaporative Fraction (Lambda) State Equation
* **2014 Version:** Formulated the state equation for the evaporative fraction (Lambda) without an explicit land-surface moisture availability constraint built into the advection-aridity equation.
* **2015 Version:** Implements a direct moisture availability constraint ($M$) into the original advection-aridity hypothesis to derive the Lambda state equation.
* **Modified Equation:** The upgraded formulation adjusts the energy partitioning using the following moisture-constrained state equation:
  `Lambda = (2 * alpha * s) / (2*s + 2*gamma + gamma * (g_B / g_s) * (1 + M))`

### 3. Moisture Availability (M) and Diurnal Hysteresis
* **2014 Version:** Neglected the diurnal hysteresis effects occurring between evapotranspiration, atmospheric vapor pressure deficit, and land surface temperature.
* **2015 Version:** Employs a two-equation system for retrieving $M$ depending on the real-time detection of clear-day hysteretic limbs.
* **Hysteresis Handling:** When a shift from energy limitation to water limitation occurs, it triggers an alternative retrieval equation based on Granger and Gray (1989) to incorporate root-zone wetness controls and vapor pressure deficit feedback.

### 4. Temporal Scale and Evaluation Datasets
* **2014 Version:** Interpreted and evaluated using coarse temporal resolution data aggregated into 8-day blocks from satellite and tower averages.
* **2015 Version:** Validated and fine-tuned for high-temporal frequencies at half-hourly and hourly scales.
* **Field Validation:** Extensively tested against continuous micrometeorological and high-frequency temperature measurements from intensive field campaigns, including SMEX02, BEAREX08, FIFE, and SAFARI2000.