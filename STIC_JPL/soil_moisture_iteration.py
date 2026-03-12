from typing import Union
import numpy as np

import rasters as rt

from rasters import Raster

from .constants import GAMMA_HPA
from .root_zone_iteration import calculate_rootzone_moisture
from .priestley_taylor_potential_evaporation import priestley_taylor_potential_evaporation

def iterate_soil_moisture(
        delta_hPa: Union[Raster, np.ndarray], # Rate of change of saturation vapor pressure with temperature (hPa/°C)
        s1: Union[Raster, np.ndarray], # The slope of SVP at surface temperature (hPa/K)
        s3: Union[Raster, np.ndarray], # The slope of SVP at dewpoint temperature (hPa/K)
        ST_C: Union[Raster, np.ndarray], # surface temperature (°C)
        Ta_C: Union[Raster, np.ndarray], # air temperature (°C)
        dTS_C: Union[Raster, np.ndarray], # difference between surface and air temperature (°C)
        Td_C: Union[Raster, np.ndarray], # dewpoint temperature (°C)
        Tsd_C: Union[Raster, np.ndarray], # surface dewpoint temperature (°C)
        Rg_Wm2: Union[Raster, np.ndarray], # Incoming solar radiation (W/m^2)
        Rn_Wm2: Union[Raster, np.ndarray], # Net radiation (W/m^2)
        LWnet_Wm2: Union[Raster, np.ndarray], # Net longwave radiation (W/m^2)
        FVC: Union[Raster, np.ndarray], # Fractional vegetation cover (unitless)
        VPD_hPa: Union[Raster, np.ndarray], # Vapor pressure deficit (hPa)
        D0_hPa: Union[Raster, np.ndarray], # Vapor pressure deficit at source (hPa)
        SVP_hPa: Union[Raster, np.ndarray], # Saturation vapor pressure (hPa)
        Ea_hPa: Union[Raster, np.ndarray], # Actual vapor pressure (hPa)
        T0: Union[Raster, np.ndarray], # Temperature at source (°C)
        gamma_hPa: Union[Raster, np.ndarray, float] = GAMMA_HPA # Psychrometric constant (hPa/°C)
    ) -> Union[Raster, np.ndarray]:
    """
    Estimates the soil moisture availability (SM) (or wetness) (value 0 to 1) based on thermal IR and meteorological
    information. However, this M will be treated as initial M, which will be later on estimated through iteration in the
    actual ET estimation loop to establish feedback between M and biophysical states.

    Args:
        delta (np.ndarray): Rate of change of saturation vapor pressure with temperature (kPa/°C).
        s1 (np.ndarray): The slope of saturation vapor pressure at surface temperature (hPa/K).
        s3 (np.ndarray): The slope of saturation vapor pressure at dewpoint temperature (hPa/K).
        ST_C (np.ndarray): Surface temperature in degrees Celsius.
        Ta_C (np.ndarray): Air temperature in degrees Celsius.
        dTS (np.ndarray): Difference between surface and air temperature in degrees Celsius.
        Td_C (np.ndarray): Dewpoint temperature in degrees Celsius.
        Tsd_C (np.ndarray): Surface dewpoint temperature in degrees Celsius.
        Rg (np.ndarray): Incoming solar radiation in W/m^2.
        Rn (np.ndarray): Net radiation in W/m^2.
        LWnet (np.ndarray): Net longwave radiation in W/m^2.
        FVC (np.ndarray): Fractional vegetation cover (unitless).
        VPD_hPa (np.ndarray): Vapor pressure deficit in hPa.
        D0 (np.ndarray): Vapor pressure deficit at source in hPa.
        SVP_hPa (np.ndarray): Saturation vapor pressure in hPa.
        Ea_hPa (np.ndarray): Actual vapor pressure in hPa.
        T0 (np.ndarray): Temperature at source in degrees Celsius.

    Returns:
        np.ndarray: Soil moisture availability (SM) (or wetness) (value 0 to 1).
    """
    # calculate surface wetness (Msurf)
    kTSTD = (T0 - Td_C) / (ST_C - Td_C)
    Msurf = (s1 / s3) * ((Tsd_C - Td_C) / (kTSTD * (ST_C - Td_C)))  # Surface wetness
    Msurf = rt.where((Rn_Wm2 < 0) & (dTS_C < 0) & (Msurf < 0), np.abs(Msurf), Msurf)
    Msurf = rt.where((Rn_Wm2 < 0) & (dTS_C > 0) & (Msurf < 0), np.abs(Msurf), Msurf)
    Msurf = rt.where((Rn_Wm2 > 0) & (dTS_C < 0) & (Msurf < 0), np.abs(Msurf), Msurf)
    Msurf = rt.where((Rn_Wm2 > 0) & (dTS_C > 0) & (Msurf < 0), np.abs(Msurf), Msurf)
    
    # solar radiation is only used to correct negative values of surface wetness
    Msurf = rt.where((Rg_Wm2 > 0) & (Msurf < 0), np.abs(Msurf), Msurf)
    Msurf = rt.where((Rg_Wm2 < 0) & (Msurf < 0), np.abs(Msurf), Msurf)
    
    Msurf = rt.where((Td_C < 0) & (Msurf < 0), np.abs(Msurf), Msurf)
    Msurf = rt.clip(Msurf, 0.0001, 1)

    # Separating soil and canopy wetness to form a composite surface moisture
    Ms = Msurf
    Mcan = FVC * Msurf
    Msoil = (1 - FVC) * Msurf

    TdewIndex = (ST_C - Tsd_C) / (Ta_C - Td_C)

    # Potential evaporation (Priestley-Taylor eqn.)
    PET_PT_Wm2 = priestley_taylor_potential_evaporation(
        delta_hPa=delta_hPa,
        energy_Wm2=Rn_Wm2,
        alpha=1.26,
        gamma_hPa=gamma_hPa
    )

    # calculate surface wetness (Ms)
    Ms = rt.where((PET_PT_Wm2 > Rn_Wm2) & (FVC <= 0.25) & (TdewIndex < 1), Msoil, Ms)
    Ms = rt.where((FVC <= 0.25) & (Ta_C > 10) & (Td_C < 0) & (LWnet_Wm2 < -125), Msoil, Ms)
    Ms = rt.where((Rn_Wm2 > PET_PT_Wm2) & (FVC <= 0.25) & (TdewIndex < 1) & (Td_C <= 0), Msoil, Ms)
    Ms = rt.where((Rn_Wm2 > PET_PT_Wm2) & (FVC <= 0.25) & (TdewIndex < 1), Msoil, Ms)
    Ms = rt.where((D0_hPa > VPD_hPa) & (FVC <= 0.25) & (TdewIndex < 1), Msoil, Ms)

    # calculate canopy wetness (Mcan)
    Mcan = rt.where((PET_PT_Wm2 > Rn_Wm2) & (FVC <= 0.25) & (TdewIndex < 1), 0, Mcan)
    Mcan = rt.where((FVC <= 0.25) & (Ta_C > 10) & (Td_C < 0) & (LWnet_Wm2 < -125), 0, Mcan)
    Mcan = rt.where((Rn_Wm2 > PET_PT_Wm2) & (FVC <= 0.25) & (TdewIndex < 1) & (Td_C <= 0), 0, Mcan)
    Mcan = rt.where((Rn_Wm2 > PET_PT_Wm2) & (FVC <= 0.25) & (TdewIndex < 1), 0, Mcan)
    Mcan = rt.where((D0_hPa > VPD_hPa) & (FVC <= 0.25) & (TdewIndex < 1), 0, Mcan)

    # calculate rootzone moisture (Mrz)
    Mrz = calculate_rootzone_moisture(
        delta_hPa=delta_hPa, 
        s1_hPa=s1, 
        s3_hPa=s3, 
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
        D0_hPa=D0_hPa, 
        SVP_hPa=SVP_hPa, 
        Ea_hPa=Ea_hPa,
        gamma_hPa=gamma_hPa
    )

    TdewIndex = (ST_C - Tsd_C) / (Ta_C - Td_C)

    # Potential evaporation (Priestley-Taylor eqn.)
    PET_PT_Wm2 = priestley_taylor_potential_evaporation(
        delta_hPa=delta_hPa,
        energy_Wm2=Rn_Wm2,
        alpha=1.26,
        gamma_hPa=gamma_hPa
    )

    # COMBINE M to account for Hysteresis and initial estimation of surface vapor pressure
    SM = Msurf
    SM = rt.where((PET_PT_Wm2 > Rn_Wm2) & (dTS_C > 0) & (FVC <= 0.25) & (D0_hPa > VPD_hPa) & (TdewIndex < 1), Mrz, SM)
    SM = rt.where((FVC <= 0.25) & (dTS_C > 0) & (Ta_C > 10) & (Td_C < 0) & (LWnet_Wm2 < -125) & (D0_hPa > VPD_hPa), Mrz, SM)

    return SM



