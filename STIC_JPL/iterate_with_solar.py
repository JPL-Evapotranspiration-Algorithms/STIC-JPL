from typing import Union, Tuple
import numpy as np

import rasters as rt

from rasters import Raster

from SEBAL_soil_heat_flux import calculate_SEBAL_soil_heat_flux

from .constants import *
from .canopy_air_stream import calculate_canopy_air_stream_vapor_pressure
from .soil_moisture_iteration import iterate_soil_moisture


def iterate_with_solar(
        seconds_of_day: Union[Raster, np.ndarray],  # Seconds of the day
        ST_C: Union[Raster, np.ndarray],  # Soil temperature (°C)
        NDVI: Union[Raster, np.ndarray],  # Normalized Difference Vegetation Index
        albedo: Union[Raster, np.ndarray],  # Albedo
        gB_ms: Union[Raster, np.ndarray],  # boundary layer conductance (m/s)
        gS_ms: Union[Raster, np.ndarray],  # stomatal conductance (m/s)
        LE_Wm2: Union[Raster, np.ndarray],  # latent heat flux (W/m^2)
        Rg_Wm2: Union[Raster, np.ndarray],  # Incoming solar radiation (W/m^2)
        Rn_Wm2: Union[Raster, np.ndarray],  # Net radiation (W/m^2)
        LWnet_Wm2: Union[Raster, np.ndarray],  # Net longwave radiation (W/m^2)
        Ta_C: Union[Raster, np.ndarray],  # Air temperature (°C)
        dTS_C: Union[Raster, np.ndarray],  # Change in soil temperature (°C)
        Td_C: Union[Raster, np.ndarray],  # Dew point temperature (°C)
        Tsd_C: Union[Raster, np.ndarray],  # Soil dew point temperature (°C)
        Ea_hPa: Union[Raster, np.ndarray],  # actual vapor pressure (hPa)
        Estar_hPa: Union[Raster, np.ndarray],  # saturation vapor pressure at surface temperature (hPa)
        VPD_hPa: Union[Raster, np.ndarray],  # Vapor pressure deficit (hPa)
        SVP_hPa: Union[Raster, np.ndarray],  # Saturation vapor pressure (hPa)
        delta_hPa: Union[Raster, np.ndarray],  # Slope of the saturation vapor pressure-temperature curve (hPa/°C)
        phi_Wm2: Union[Raster, np.ndarray],  # Net radiation minus soil heat flux (W/m^2)
        Es_hPa: Union[Raster, np.ndarray],  # Saturation vapor pressure (hPa)
        s1: Union[Raster, np.ndarray],  # Soil moisture parameter
        s3: Union[Raster, np.ndarray],  # Soil moisture parameter
        FVC: Union[Raster, np.ndarray],  # Fractional canopy cover
        T0_C: Union[Raster, np.ndarray],  # Reference temperature (°C)
        gB_by_gS: Union[Raster, np.ndarray],  # Ratio of boundary layer conductance to stomatal conductance
        rho_kgm3: Union[Raster, np.ndarray, float] = RHO_KGM3,  # Air density (kg/m^3)
        Cp_Jkg: Union[Raster, np.ndarray, float] = CP_JKG,  # Specific heat at constant pressure (J/kg/K)
        gamma_hPa: Union[Raster, np.ndarray, float] = GAMMA_HPA,  # Psychrometric constant (hPa/°C)
        G_method: str = "santanello"  # Method for calculating soil heat flux
    ) -> Tuple[Union[Raster, np.ndarray]]:
    """
    This function calculates the canopy air stream vapor pressures, vapor pressure deficit at source,
    vapor pressure at the reference height, soil moisture, soil heat flux, and recompute phi.

    Parameters:
    seconds_of_day (np.ndarray): Seconds of the day
    ST_C (np.ndarray): Soil temperature (°C)
    NDVI (np.ndarray): Normalized Difference Vegetation Index
    albedo (np.ndarray): Albedo
    gB (np.ndarray): Boundary layer conductance (m/s)
    gS (np.ndarray): Stomatal conductance (m/s)
    LE (np.ndarray): Latent heat flux (W/m^2)
    Rg (np.ndarray): Incoming solar radiation (W/m^2)
    Rn (np.ndarray): Net radiation (W/m^2)
    LWnet (np.ndarray): Net longwave radiation (W/m^2)
    Ta_C (np.ndarray): Air temperature (°C)
    dTS (np.ndarray): Change in soil temperature (°C)
    Td_C (np.ndarray): Dew point temperature (°C)
    Tsd_C (np.ndarray): Soil dew point temperature (°C)
    Ea_hPa (np.ndarray): Actual vapor pressure (hPa)
    Estar (np.ndarray): Saturation vapor pressure at surface temperature (hPa)
    VPD_hPa (np.ndarray): Vapor pressure deficit (hPa)
    SVP_hPa (np.ndarray): Saturation vapor pressure (hPa)
    delta (np.ndarray): Slope of the saturation vapor pressure-temperature curve (hPa/°C)
    phi (np.ndarray): Net radiation minus soil heat flux (W/m^2)
    Es (np.ndarray): Saturation vapor pressure (hPa)
    s1 (np.ndarray): Soil moisture parameter
    s3 (np.ndarray): Soil moisture parameter
    fc (np.ndarray): Fractional canopy cover
    T0 (np.ndarray): Reference temperature (°C)
    gB_by_gS (np.ndarray): Ratio of boundary layer conductance to stomatal conductance
    gamma (np.ndarray): Psychrometric constant (hPa/°C)
    rho (np.ndarray): Air density (kg/m^3)
    cp (np.ndarray): Specific heat at constant pressure (J/kg/K)
    G_method (str): Method for calculating soil heat flux (santanello or SEBAL)

    Returns:
    e0star (np.ndarray): Canopy air stream vapor pressures
    D0 (np.ndarray): Vapor pressure deficit at source
    e0 (np.ndarray): Vapor pressure at the reference height
    SM (np.ndarray): Soil moisture
    G (np.ndarray): Soil heat flux
    """
    # canopy air stream vapor pressures
    e0star = calculate_canopy_air_stream_vapor_pressure(
        LE=LE_Wm2,
        Ea_hPa=Ea_hPa,
        Estar=Estar_hPa,
        gB=gB_ms,
        gS=gS_ms
    )

    # vapor pressure deficit at source
    D0 = VPD_hPa + (delta_hPa * phi_Wm2 - (delta_hPa + gamma_hPa) * LE_Wm2) / (rho_kgm3 * Cp_Jkg * gB_ms)

    # Vapor pressure at the reference height (hPa)
    e0 = e0star - D0  
    e0 = rt.where(e0 < 0, Es_hPa, e0)
    e0 = rt.where(e0 > e0star, e0star, e0)

    # calculate soil moisture
    SM = iterate_soil_moisture(
        delta_hPa=delta_hPa,
        s1=s1,
        s3=s3,
        ST_C=ST_C,
        Ta_C=Ta_C,
        dTS_C=dTS_C,
        Td_C=Td_C,
        Tsd_C=Tsd_C,
        Rg_Wm2=Rg_Wm2,
        Rn_Wm2=Rn_Wm2,
        LWnet_Wm2=LWnet_Wm2,
        FVC=FVC,
        VPD_hPa=VPD_hPa,
        D0_hPa=D0,
        SVP_hPa=SVP_hPa,
        Ea_hPa=Ea_hPa,
        T0=T0_C,
        gamma_hPa=gamma_hPa
    )

    # calculate soil heat flux
    G = calculate_SEBAL_soil_heat_flux(
        ST_C=ST_C,
        NDVI=NDVI,
        albedo=albedo,
        Rn=Rn_Wm2
    )

    # recompute phi
    phi_Wm2 = Rn_Wm2 - G

    #
    alphaN = ((gS_ms * (e0star - Ea_hPa) * (2 * delta_hPa + 2 * gamma_hPa + gamma_hPa * gB_by_gS * (1 + SM))) / (
                2 * delta_hPa * (gamma_hPa * (T0_C - Ta_C) * (gB_ms + gS_ms) + gS_ms * (e0star - Ea_hPa))))

    return SM, G, e0, e0star, D0, alphaN
