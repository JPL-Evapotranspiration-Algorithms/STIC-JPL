from typing import Union, Tuple
import numpy as np

import rasters as rt

from rasters import Raster

from SEBAL_soil_heat_flux import calculate_SEBAL_soil_heat_flux

from .constants import *
from .soil_moisture_initialization import initialize_soil_moisture
from .net_radiation import calculate_net_longwave_radiation


def initialize_with_solar(
        seconds_of_day: Union[Raster, np.ndarray],  # time of day in seconds since midnight
        Rg_Wm2: Union[Raster, np.ndarray],  # solar radiation (W/m^2)
        Rn_Wm2: Union[Raster, np.ndarray],  # net radiation (W/m^2)
        ST_C: Union[Raster, np.ndarray],  # surface temperature (Celsius)
        emissivity: Union[Raster, np.ndarray],  # emissivity of the surface
        Ta_C: Union[Raster, np.ndarray],  # air temperature (Celsius)
        dTS_C: Union[Raster, np.ndarray],  # surface air temperature difference (Celsius)
        Td_C: Union[Raster, np.ndarray],  # dew point temperature (Celsius)
        VPD_hPa: Union[Raster, np.ndarray],  # vapor pressure deficit (hPa)
        SVP_hPa: Union[Raster, np.ndarray],  # saturation vapor pressure at given air temperature (hPa)
        Ea_hPa: Union[Raster, np.ndarray],  # actual vapor pressure at air temperature (hPa)
        Estar_hPa: Union[Raster, np.ndarray],  # saturation vapor pressure at surface temperature (hPa)
        delta_hPa: Union[Raster, np.ndarray],  # slope of saturation vapor pressure to air temperature (hpa/K)
        NDVI: Union[Raster, np.ndarray],  # normalized difference vegetation index
        FVC: Union[Raster, np.ndarray],  # fractional vegetation cover
        LAI: Union[Raster, np.ndarray],  # leaf area index
        albedo: Union[Raster, np.ndarray],  # albedo of the surface
        gamma_hPa: Union[Raster, np.ndarray, float] = GAMMA_HPA,  # psychrometric constant (hPa/°C)
        G_method: str = DEFAULT_G_METHOD  # method for calculating soil heat flux
    ) -> Tuple[Union[Raster, np.ndarray]]:
    # Rn SOIL
    kRN = 0.6
    Rn_soil = Rn_Wm2 * np.exp(-kRN * LAI)

    LWnet = calculate_net_longwave_radiation(
        Ta_C=Ta_C, 
        Ea_hPa=Ea_hPa, 
        ST_C=ST_C, 
        emissivity=emissivity, 
        albedo=albedo
    )
    
    # initialize soil moisture
    SM, SMrz, Ms, Ep_PT, Ds, s1, s3, Tsd_C = initialize_soil_moisture(
        delta_hPa=delta_hPa,  # slope of saturation vapor pressure to air temperature (hpa/K)
        ST_C=ST_C,  # surface temperature (Celsius)
        Ta_C=Ta_C,  # air temperature (Celsius)
        Td_C=Td_C,  # dew point temperature (Celsius)
        dTS=dTS_C,  # surface air temperature difference (Celsius)
        Rn_Wm2=Rn_Wm2,  # net radiation (W/m^2)
        LWnet_Wm2=LWnet,  # net longwave radiation (W/m^2)
        FVC=FVC,  # fractional vegetation cover
        VPD_hPa=VPD_hPa,  # vapor pressure deficit (hPa)
        SVP_hPa=SVP_hPa,  # saturation vapor pressure at given air temperature (hPa)
        Ea_hPa=Ea_hPa,  # actual vapor pressure at air temperature (hPa)
        Estar_hPa=Estar_hPa,  # saturation vapor pressure at surface temperature (hPa)
        gamma_hPa=gamma_hPa # psychrometric constant (hPa/°C)
    )

    # calculate soil heat flux
    G = calculate_SEBAL_soil_heat_flux(
        ST_C=ST_C,  # Surface temperature in Celsius
        NDVI=NDVI,  # Normalized Difference Vegetation Index
        albedo=albedo,  # Albedo of the surface
        Rn=Rn_Wm2  # Net radiation (W/m^2)
    )

    # get phi with new comp
    phi = Rn_Wm2 - G

    Es = rt.where((Ep_PT > phi) & (dTS_C > 0) & (Td_C <= 0), Ea_hPa + SMrz * (Estar_hPa - Ea_hPa),
                    Ea_hPa + Ms * (Estar_hPa - Ea_hPa))
    
    # Return all the created variables
    return SM, SMrz, Ms, s1, s3, Ep_PT, Rn_soil, LWnet, G, Tsd_C, Ds, Es, phi
