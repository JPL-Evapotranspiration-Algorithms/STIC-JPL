from typing import Union, Tuple
import numpy as np

import rasters as rt

from rasters import Raster

from .constants import *

def STIC_closure(
        delta_hPa: Union[Raster, np.ndarray],  # Slope of the saturation vapor pressure-temperature curve (hPa/°C)
        phi_Wm2: Union[Raster, np.ndarray],  # available energy (W/m^2)
        Es_hPa: Union[Raster, np.ndarray],  # Vapor pressure at the reference height (hPa)
        Ea_hPa: Union[Raster, np.ndarray],  # Actual vapor pressure (hPa)
        Estar_hPa: Union[Raster, np.ndarray],  # Saturation vapor pressure at the reference height (hPa)
        SM: Union[Raster, np.ndarray],  # Soil moisture (m³/m³)
        gamma_hPa: float = GAMMA_HPA,  # Psychrometric constant (hPa/°C)
        rho_kgm3: float = RHO_KGM3,  # Air density (kg/m³)
        Cp_Jkg: float = CP_JKG,  # Specific heat capacity of air (J/kg/°C)
        alpha: float = PT_ALPHA  # Priestley-Taylor alpha
    ) -> Tuple[Union[Raster, np.ndarray]]:
    """
    STIC closure equations with modified Priestley Taylor and Penman Monteith
    (Mallick et al., 2015, Water Resources research)

    Parameters:
    delta (np.ndarray): Slope of the saturation vapor pressure-temperature curve (hPa/°C)
    phi (np.ndarray): available energy (W/m^2)
    Es (np.ndarray): Vapor pressure at the reference height (hPa)
    Ea (np.ndarray): Actual vapor pressure (hPa)
    Estar (np.ndarray): Saturation vapor pressure at the reference height (hPa)
    SM (np.ndarray): Soil moisture (m³/m³)
    rho (float, optional): Air density (kg/m³). Defaults to RHO.
    cp (float, optional): Specific heat capacity of air (J/kg/°C). Defaults to CP.
    gamma (float, optional): Psychrometric constant (hPa/°C). Defaults to GAMMA.
    alpha (float, optional): Stability correction factor for conductance. Defaults to ALPHA.

    Returns:
    gB (np.ndarray): boundary layer conductance (m s^-1)
    gS (np.ndarray): stomatal conductance (m s^-1)
    dT (np.ndarray): difference between surface and air temperature
    EF (np.ndarray): evaporative fraction
    """

    epsilon = 1e-8  # Small value to prevent division by zero

    # boundary layer conductance (m s^-1)
    gB_den = (2 * Cp_Jkg * delta_hPa * Es_hPa * rho_kgm3 - 2 * Cp_Jkg * delta_hPa * Ea_hPa * rho_kgm3 - 2 * Cp_Jkg * Ea_hPa * gamma_hPa * rho_kgm3 + Cp_Jkg * Es_hPa * gamma_hPa * rho_kgm3 + Cp_Jkg * Estar_hPa * gamma_hPa * rho_kgm3 - Cp_Jkg * SM * Es_hPa * gamma_hPa * rho_kgm3 + Cp_Jkg * SM * Estar_hPa * gamma_hPa * rho_kgm3)
    gB = ((2 * phi_Wm2 * alpha * delta_hPa * gamma_hPa) / (gB_den + epsilon))
    gB = rt.clip(gB, 0.0001, 0.2)

    # stomatal conductance (m s^-1)
    gS_den = (Cp_Jkg * Estar_hPa ** 2 * gamma_hPa * rho_kgm3 - Cp_Jkg * Es_hPa ** 2 * gamma_hPa * rho_kgm3 - 2 * Cp_Jkg * delta_hPa * Es_hPa ** 2 * rho_kgm3 + 2 * Cp_Jkg * delta_hPa * Ea_hPa * Es_hPa * rho_kgm3 - 2 * Cp_Jkg * delta_hPa * Ea_hPa * Estar_hPa * rho_kgm3 + 2 * Cp_Jkg * delta_hPa * Es_hPa * Estar_hPa * rho_kgm3 + 2 * Cp_Jkg * Ea_hPa * Es_hPa * gamma_hPa * rho_kgm3 - 2 * Cp_Jkg * Ea_hPa * Estar_hPa * gamma_hPa * rho_kgm3 + Cp_Jkg * SM * Es_hPa ** 2 * gamma_hPa * rho_kgm3 + Cp_Jkg * SM * Estar_hPa ** 2 * gamma_hPa * rho_kgm3 - 2 * Cp_Jkg * SM * Es_hPa * Estar_hPa * gamma_hPa * rho_kgm3)
    gS = (-(2 * (phi_Wm2 * alpha * delta_hPa * Ea_hPa * gamma_hPa - phi_Wm2 * alpha * delta_hPa * Es_hPa * gamma_hPa)) / (gS_den + epsilon))
    gS = rt.clip(gS, 0.0001, 0.2)

    # difference between surface and air temperature (Celsius)
    dT_num = (2 * delta_hPa * Es_hPa - 2 * delta_hPa * Ea_hPa - 2 * Ea_hPa * gamma_hPa + Es_hPa * gamma_hPa + Estar_hPa * gamma_hPa - SM * Es_hPa * gamma_hPa + SM * Estar_hPa * gamma_hPa + 2 * alpha * delta_hPa * Ea_hPa - 2 * alpha * delta_hPa * Es_hPa)
    dT_den = (2 * alpha * delta_hPa * gamma_hPa)
    dT = dT_num / (dT_den + epsilon)
    dT = rt.clip(dT, -10, 50)

    # evaporative fraction
    EF_den = (2 * delta_hPa * Es_hPa - 2 * delta_hPa * Ea_hPa - 2 * Ea_hPa * gamma_hPa + Es_hPa * gamma_hPa + Estar_hPa * gamma_hPa - SM * Es_hPa * gamma_hPa + SM * Estar_hPa * gamma_hPa)
    EF = -(2 * alpha * delta_hPa * Ea_hPa - 2 * alpha * delta_hPa * Es_hPa) / (EF_den + epsilon)
    EF = rt.clip(EF, 0, 1)

    return gB, gS, dT, EF
