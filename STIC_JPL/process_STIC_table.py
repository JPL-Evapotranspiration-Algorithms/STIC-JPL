import logging

import numpy as np
from dateutil import parser
from pandas import DataFrame
from rasters import Point
from sentinel_tiles import sentinel_tiles
from solar_apparent_time import UTC_to_solar
from SEBAL_soil_heat_flux import calculate_SEBAL_soil_heat_flux

from .constants import *
from .model import STIC_JPL, MAX_ITERATIONS, USE_VARIABLE_ALPHA

logger = logging.getLogger(__name__)

def process_STIC_table(
        input_df: DataFrame, 
        max_iterations = MAX_ITERATIONS, 
        use_variable_alpha = USE_VARIABLE_ALPHA,
        constrain_negative_LE = CONSTRAIN_NEGATIVE_LE,
        supply_SWin = SUPPLY_SWIN
    ) -> DataFrame:
    ST_C = np.float64(np.array(input_df.ST_C))
    emissivity = np.float64(np.array(input_df.EmisWB))
    NDVI = np.float64(np.array(input_df.NDVI))
    albedo = np.float64(np.array(input_df.albedo))
    Ta_C = np.float64(np.array(input_df.Ta_C))
    RH = np.float64(np.array(input_df.RH))
    Rn_Wm2 = np.float64(np.array(input_df.Rn_Wm2))
    Rg = np.float64(np.array(input_df.Rg))

    if "G" in input_df:
        G_Wm2 = np.array(input_df.G)
    else:
        G_Wm2 = calculate_SEBAL_soil_heat_flux(
            Rn=Rn_Wm2,
            ST_C=ST_C,
            NDVI=NDVI,
            albedo=albedo
        )
    
    if "SWin_Wm2" in input_df and supply_SWin:
        SWin_Wm2 = np.float64(np.array(input_df.SWin_Wm2))
    else:
        SWin_Wm2 = None

    # --- Handle geometry and time columns ---
    import pandas as pd
    from rasters import MultiPoint, WGS84
    from shapely.geometry import Point

    def ensure_geometry(df):
        if "geometry" in df:
            if isinstance(df.geometry.iloc[0], str):
                def parse_geom(s):
                    s = s.strip()
                    if s.startswith("POINT"):
                        coords = s.replace("POINT", "").replace("(", "").replace(")", "").strip().split()
                        return Point(float(coords[0]), float(coords[1]))
                    elif "," in s:
                        coords = [float(c) for c in s.split(",")]
                        return Point(coords[0], coords[1])
                    else:
                        coords = [float(c) for c in s.split()]
                        return Point(coords[0], coords[1])
                df = df.copy()
                df['geometry'] = df['geometry'].apply(parse_geom)
        return df

    input_df = ensure_geometry(input_df)

    logger.info("started extracting geometry from PT-JPL-SM input table")

    if "geometry" in input_df:
        # Convert Point objects to coordinate tuples for MultiPoint
        if hasattr(input_df.geometry.iloc[0], "x") and hasattr(input_df.geometry.iloc[0], "y"):
            coords = [(pt.x, pt.y) for pt in input_df.geometry]
            geometry = MultiPoint(coords, crs=WGS84)
        else:
            geometry = MultiPoint(input_df.geometry, crs=WGS84)
    elif "lat" in input_df and "lon" in input_df:
        lat = np.array(input_df.lat).astype(np.float64)
        lon = np.array(input_df.lon).astype(np.float64)
        geometry = MultiPoint(x=lon, y=lat, crs=WGS84)
    else:
        raise KeyError("Input DataFrame must contain either 'geometry' or both 'lat' and 'lon' columns.")

    logger.info("completed extracting geometry from PT-JPL-SM input table")

    logger.info("started extracting time from PT-JPL-SM input table")
    time_UTC = pd.to_datetime(input_df.time_UTC).tolist()
    logger.info("completed extracting time from PT-JPL-SM input table")
    
    results = STIC_JPL(
        geometry=geometry,
        ST_C = ST_C,
        emissivity=emissivity,
        NDVI=NDVI,
        albedo=albedo,
        Ta_C=Ta_C,
        RH=RH,
        Rn_Wm2=Rn_Wm2,
        G_Wm2=G_Wm2,
        SWin_Wm2=SWin_Wm2,
        time_UTC=time_UTC,
        max_iterations=max_iterations,
        use_variable_alpha=use_variable_alpha,
        constrain_negative_LE=constrain_negative_LE
    )

    output_df = input_df.copy()

    for key, value in results.items():
        output_df[key] = value

    return output_df
