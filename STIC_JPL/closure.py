from typing import Union, Tuple
import numpy as np

import rasters as rt

from rasters import Raster

from .constants import *

def STIC_closure(
        delta_hPa: Union[Raster, np.ndarray],  # Slope of the saturation vapor pressure-temperature curve (hPa/°C)
        phi_Wm2: Union[Raster, np.ndarray],  # available energy (W/m^2)
        Es_hPa: Union[Raster, np.ndarray],  # Vapor pressure at the evaporating front (hPa)
        Ea_hPa: Union[Raster, np.ndarray],  # Actual vapor pressure (hPa)
        Estar_hPa: Union[Raster, np.ndarray],  # Saturation vapor pressure at surface temperature (hPa)
        SM: Union[Raster, np.ndarray],  # Soil moisture / moisture availability M (0-1)
        gamma_hPa: float = GAMMA_HPA,  # Psychrometric constant (hPa/°C)
        rho_kgm3: float = RHO_KGM3,  # Air density (kg/m³)
        Cp_Jkg: float = CP_JKG,  # Specific heat capacity of air (J/kg/°C)
        alpha: float = PT_ALPHA,  # Priestley-Taylor alpha
        EF_moisture_factor: Union[Raster, np.ndarray, float, None] = None,
    ) -> Tuple[Union[Raster, np.ndarray]]:
    """
    STIC analytical closure for conductances, aerodynamic temperature, and evaporative fraction.

    All four outputs are derived from a single EF expression, parameterized by
    EF_moisture_factor to support both the 2014 and 2015/2016 formulations:

        EF = 2α δ (Es-Ea) / [ (Es-Ea)(2δ+2γ) + γ (1 + EF_moisture_factor)(Estar-Es) ]

    - EF_moisture_factor = 0   → strict 2014 formula (Mallick et al. 2014, Eq. 17)
    - EF_moisture_factor = SM  → moisture-constrained 2015/2016 formula (Mallick et al. 2015, Eq. 5)
    - EF_moisture_factor = None → defaults to SM (preserves 2015-style behavior)

    From EF, gB, gS, and dT follow from the same state equations for both versions:
        gB = Φ · EF · γ / (ρ cp (Es-Ea))        [Eq. 14 / Eq. 4 in 2014/2015]
        gS = Φ · EF · γ / (ρ cp (Estar-Es))      [Eq. 15 / Eq. 4]
        dT = (Es-Ea)/γ · (1-EF)/EF               [Eq. 16 / Eq. 4]

    Parameters
    ----------
    EF_moisture_factor : array-like or float or None
        Moisture factor in the EF denominator. Pass 0 for strict 2014 closure,
        SM (or None, which defaults to SM) for 2015/2016 closure.

    Returns
    -------
    gB, gS, dT, EF
    """

    epsilon = 1e-8

    if EF_moisture_factor is None:
        EF_moisture_factor = SM

    # Intermediate vapor gradients reused across all four outputs.
    vapor_gradient = Es_hPa - Ea_hPa        # Es - Ea  (> 0 under typical evaporation)
    vapor_deficit_sfc = Estar_hPa - Es_hPa  # Estar - Es  (> 0 when surface sub-saturated)

    # Shared denominator for EF, dT, gB, and gS.
    # 2014: EF_moisture_factor = 0  → (1+0)  = 1
    # 2015: EF_moisture_factor = SM → (1+SM)
    shared_den = (
        vapor_gradient * (2 * delta_hPa + 2 * gamma_hPa)
        + gamma_hPa * (1 + EF_moisture_factor) * vapor_deficit_sfc
    )

    # Evaporative fraction (Lambda) — version-parameterized formula.
    EF = rt.clip(
        2 * alpha * delta_hPa * vapor_gradient / (shared_den + epsilon),
        0, 1,
    )

    # Aerodynamic-source to air temperature difference: dT = (Es-Ea)/γ · (1-EF)/EF
    # Equivalent to (shared_den - EF_num) / (2·α·δ·γ) after substitution.
    dT_num = (2 * delta_hPa + 2 * gamma_hPa - 2 * alpha * delta_hPa) * vapor_gradient + gamma_hPa * (1 + EF_moisture_factor) * vapor_deficit_sfc
    dT = rt.clip(dT_num / (2 * alpha * delta_hPa * gamma_hPa + epsilon), -10, 50)

    # Boundary layer conductance: gB = 2 α δ γ Φ / (ρ cp · shared_den)
    # Equivalent to Φ · EF · γ / (ρ cp (Es-Ea)) but epsilon placement matches
    # the original polynomial form to preserve numerical behavior.
    gB_denom = rho_kgm3 * Cp_Jkg * shared_den
    gB = rt.clip(
        2 * phi_Wm2 * alpha * delta_hPa * gamma_hPa / (gB_denom + epsilon),
        0.0001, 0.2,
    )

    # Surface/stomatal conductance: gS = gB · (Es-Ea) / (Estar-Es)
    gS_denom = rho_kgm3 * Cp_Jkg * vapor_deficit_sfc * shared_den
    gS = rt.clip(
        2 * phi_Wm2 * alpha * delta_hPa * gamma_hPa * vapor_gradient / (gS_denom + epsilon),
        0.0001, 0.2,
    )

    return gB, gS, dT, EF
