from typing import Union, Tuple
import numpy as np

import rasters as rt
from rasters import Raster

from .FVC_from_NDVI import FVC_from_NDVI
from .LAI_from_NDVI import LAI_from_NDVI

from .constants import GAMMA_HPA
from .root_zone_initialization import calculate_root_zone_moisture

def initialize_soil_moisture(
        delta_hPa: Union[Raster, np.ndarray],  # Rate of change of saturation vapor pressure with temperature (hPa/°C)
        ST_C: Union[Raster, np.ndarray],  # Surface temperature (°C)
        Ta_C: Union[Raster, np.ndarray],  # Air temperature (°C)
        Td_C: Union[Raster, np.ndarray],  # Dewpoint temperature (°C)
        dTS: Union[Raster, np.ndarray],  # Difference between surface and air temperature (°C)
        Rn_Wm2: Union[Raster, np.ndarray],  # Net radiation (W/m²)
        LWnet_Wm2: Union[Raster, np.ndarray],  # Net longwave radiation (W/m²)
        FVC: Union[Raster, np.ndarray],  # Fractional vegetation cover (unitless)
        VPD_hPa: Union[Raster, np.ndarray],  # Vapor pressure deficit (hPa)
        SVP_hPa: Union[Raster, np.ndarray],  # Saturation vapor pressure (hPa)
        Ea_hPa: Union[Raster, np.ndarray],  # Actual vapor pressure (hPa)
        Estar_hPa: Union[Raster, np.ndarray],  # Saturation vapor pressure at surface temperature (hPa)
        gamma_hPa: Union[Raster, np.ndarray, float] = GAMMA_HPA  # Psychrometric constant (hPa/°C)
    ) -> Tuple[Union[Raster, np.ndarray]]:
    """
    This function estimates the soil moisture availability (SM) based on thermal IR and meteorological information. 
    The estimated SM is treated as initial SM, which will be later estimated through iteration in the actual ET estimation loop to establish feedback between SM and biophysical states.

    Parameters:
    delta (np.ndarray): Rate of change of saturation vapor pressure with temperature (kPa/°C)
    ST_C (np.ndarray): Surface temperature (°C)
    Ta_C (np.ndarray): Air temperature (°C)
    Td_C (np.ndarray): Dewpoint temperature (°C)
    dTS (np.ndarray): Difference between surface and air temperature (°C)
    Rn (np.ndarray): Net radiation (W/m²)
    LWnet (np.ndarray): Net longwave radiation (W/m²)
    NDVI (np.ndarray): Normalized Difference Vegetation Index (unitless)
    VPD_hPa (np.ndarray): Vapor pressure deficit (hPa)
    SVP_hPa (np.ndarray): Saturation vapor pressure (hPa)
    Ea_hPa (np.ndarray): Actual vapor pressure (hPa)
    Estar (np.ndarray): Saturation vapor pressure at surface temperature (hPa)

    Returns:
    SM (np.ndarray): soil moisture (m³/m³)
    Mrz (np.ndarray): The rootzone moisture (m³/m³)
    Ms (np.ndarray): The surface moisture (m³/m³)
    Ep_PT (np.ndarray): The potential evaporation (Priestley-Taylor eqn.) (mm/day)
    Ds (np.ndarray): The vapor pressure deficit at surface (hPa)
    s1 (np.ndarray): The slope of SVP at dewpoint temperature (hPa/K)
    s3 (np.ndarray): The slope of SVP at surface temperature (hPa/K)
    Tsd_C (np.ndarray): The surface dewpoint temperature (°C)
    """
    # Compute the surface dewpoint temperature
    s11 = (45.03 + 3.014 * Td_C + 0.05345 * Td_C ** 2 + 0.00224 * Td_C ** 3) * 1e-2  # slope of SVP at TD (hpa/K)
    s22 = (Estar_hPa - Ea_hPa) / (ST_C - Td_C)
    s33 = (45.03 + 3.014 * ST_C + 0.05345 * ST_C ** 2 + 0.00224 * ST_C ** 3) * 1e-2  # slope of SVP at TS (hpa/K)
    s44 = (SVP_hPa - Ea_hPa) / (Ta_C - Td_C)
    Tsd_C = (Estar_hPa - Ea_hPa - s33 * ST_C + s11 * Td_C) / (s11 - s33)  # Surface dewpoint temperature (degC)

    # Calculate the surface moisture or wetness
    Msurf = (s11 / s22) * ((Tsd_C - Td_C) / (ST_C - Td_C))  # Surface wetness
    Msurf = rt.clip(Msurf, 0.0001, 0.9999)

    # Calculate the surface vapor pressure and deficit
    esurf = Ea_hPa + Msurf * (Estar_hPa - Ea_hPa)
    Dsurf = esurf - Ea_hPa

    # Separate the soil and canopy wetness to form a composite surface moisture
    Ms = Msurf
    Mcan = FVC * Msurf
    Msoil = (1 - FVC) * Msurf

    TdewIndex = (ST_C - Tsd_C) / (Ta_C - Td_C)  # % TdewIndex > 1 signifies super dry condition
    Ep_PT = (1.26 * delta_hPa * Rn_Wm2) / (delta_hPa + gamma_hPa)  # Potential evaporation (Priestley-Taylor eqn.)

    # Adjust surface wetness based on certain conditions
    Ms = rt.where((FVC <= 0.25) & (TdewIndex < 1), Msoil, Ms)
    Mcan = rt.where((FVC <= 0.25) & (TdewIndex < 1), 0, Mcan)
    Ms = rt.where((FVC <= 0.25) & (Ta_C > 10) & (Td_C < 0) & (LWnet_Wm2 < -125), Msoil, Ms)
    Mcan = rt.where((FVC <= 0.25) & (Ta_C > 10) & (Td_C < 0) & (LWnet_Wm2 < -125), 0, Mcan)
    
    # Calculate the root-zone moisture
    SMrz = calculate_root_zone_moisture(
        delta_hPa=delta_hPa, 
        ST_C=ST_C, 
        Ta_C=Ta_C, 
        Td_C=Td_C, 
        s11=s11, 
        s33=s33, 
        s44=s44, 
        Tsd_C=Tsd_C,
        gamma_hPa=gamma_hPa
    )

    # Combine soil moisture to account for hysteresis and initial estimation of surface vapor pressure
    SM = Ms
    SM = rt.where((Ep_PT > Rn_Wm2) & (dTS > 0), SMrz, SM)
    SM = rt.where((Ep_PT > Rn_Wm2) & (FVC <= 0.25), SMrz, SM)
    SM = rt.where((Ep_PT > Rn_Wm2) & (Dsurf > VPD_hPa), SMrz, SM)
    SM = rt.where((FVC <= 0.25) & (dTS > 0) & (Ta_C > 10) & (Td_C < 0) & (LWnet_Wm2 < -125), SMrz, SM)
    SM = rt.where((FVC <= 0.25) & (dTS > 0) & (Ta_C > 10) & (Td_C < 0) & (Dsurf > VPD_hPa), SMrz, SM)
    SM = rt.where((Ep_PT < Rn_Wm2) & (FVC <= 0.25) & (Dsurf > VPD_hPa), SMrz, SM)
    
    es = Ea_hPa + SM * (Estar_hPa - Ea_hPa)
    
    # vapor pressure deficit at surface
    Ds = (Estar_hPa - es)  
    
    s1 = s11
    s3 = s33

    return SM, SMrz, Ms, Ep_PT, Ds, s1, s3, Tsd_C
