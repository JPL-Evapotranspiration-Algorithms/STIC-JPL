import logging

import numpy as np
from dateutil import parser
from pandas import DataFrame
from rasters import Point
from sentinel_tiles import sentinel_tiles
from solar_apparent_time import UTC_to_solar
from SEBAL_soil_heat_flux import calculate_SEBAL_soil_heat_flux

from .model import STIC_JPL, MAX_ITERATIONS, USE_VARIABLE_ALPHA

logger = logging.getLogger(__name__)

def process_STIC_table(
        input_df: DataFrame, 
        max_iterations = MAX_ITERATIONS, 
        use_variable_alpha = USE_VARIABLE_ALPHA) -> DataFrame:
    hour_of_day = np.float64(np.array(input_df.hour_of_day))
    lon = np.float64(np.array(input_df.lon))
    ST_C = np.float64(np.array(input_df.ST_C))
    emissivity = np.float64(np.array(input_df.EmisWB))
    NDVI = np.float64(np.array(input_df.NDVI))
    albedo = np.float64(np.array(input_df.albedo))
    Ta_C = np.float64(np.array(input_df.Ta_C))
    RH = np.float64(np.array(input_df.RH))
    Rn = np.float64(np.array(input_df.Rn))
    Rg = np.float64(np.array(input_df.Rg))

    if "G" in input_df:
        G = np.array(input_df.G)
    else:
        G = calculate_SEBAL_soil_heat_flux(
            Rn=Rn,
            ST_C=ST_C,
            NDVI=NDVI,
            albedo=albedo
        )
    
    results = STIC_JPL(
        hour_of_day=hour_of_day,
        # longitude=lon,
        ST_C = ST_C,
        emissivity=emissivity,
        NDVI=NDVI,
        albedo=albedo,
        Ta_C=Ta_C,
        RH=RH,
        Rn_Wm2=Rn,
        G_Wm2=G,
        # Rg_Wm2=Rg,
        max_iterations=max_iterations,
        use_variable_alpha=use_variable_alpha
    )

    output_df = input_df.copy()

    for key, value in results.items():
        output_df[key] = value

    return output_df
