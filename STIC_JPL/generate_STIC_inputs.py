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

def generate_STIC_inputs(STIC_inputs_from_calval_df: DataFrame) -> DataFrame:
    """
    STIC_inputs_from_claval_df:
        Pandas DataFrame containing the columns: tower, lat, lon, time_UTC, albedo, elevation_km
    return:
        Pandas DataFrame containing the columns: tower, lat, lon, time_UTC, doy, albedo, elevation_km, AOT, COT, vapor_gccm, ozone_cm, SZA, KG
    """
    # output_rows = []
    STIC_inputs_df = STIC_inputs_from_calval_df.copy()

    hour_of_day = []
    doy = []
    Topt = []
    fAPARmax = []

    for i, input_row in STIC_inputs_from_calval_df.iterrows():
        tower = input_row.tower
        lat = input_row.lat
        lon = input_row.lon
        time_UTC = input_row.time_UTC
        albedo = input_row.albedo
        elevation_km = input_row.elevation_km
        logger.info(f"collecting STIC inputs for tower {tower} lat {lat} lon {lon} time {time_UTC} UTC")
        time_UTC = parser.parse(str(time_UTC))
        time_solar = UTC_to_solar(time_UTC, lon)
        hour_of_day.append(time_solar.hour)
        doy.append(time_UTC.timetuple().tm_yday)
        date_UTC = time_UTC.date()
        tile = sentinel_tiles.toMGRS(lat, lon)[:5]
        tile_grid = sentinel_tiles.grid(tile=tile, cell_size=70)
        rows, cols = tile_grid.shape
        row, col = tile_grid.index_point(Point(lon, lat))
        geometry = tile_grid[max(0, row - 1):min(row + 2, rows - 1),
                             max(0, col - 1):min(col + 2, cols - 1)]
    
    if not "hour_of_day" in STIC_inputs_df.columns:
        STIC_inputs_df["hour_of_day"] = hour_of_day

    if not "doy" in STIC_inputs_df.columns:
        STIC_inputs_df["doy"] = doy
    
    if not "Topt" in STIC_inputs_df.columns:
        STIC_inputs_df["Topt"] = Topt
    
    if not "fAPARmax" in STIC_inputs_df.columns:
        STIC_inputs_df["fAPARmax"] = fAPARmax
    
    if "Ta" in STIC_inputs_df and "Ta_C" not in STIC_inputs_df:
        STIC_inputs_df.rename({"Ta": "Ta_C"}, inplace=True)
    
    return STIC_inputs_df
