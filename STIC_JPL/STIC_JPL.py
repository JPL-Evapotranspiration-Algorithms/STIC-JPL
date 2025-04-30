from typing import Union, Callable
import logging
from datetime import datetime, timedelta
from os.path import join, abspath, expanduser
from typing import Dict, List
import numpy as np
import warnings
from .diagnostic import diagnostic
import colored_logging as cl
from .meteorology_conversion import calculate_air_density, calculate_specific_heat, calculate_specific_humidity, calculate_surface_pressure, celcius_to_kelvin
import rasters as rt
from GEOS5FP import GEOS5FP
from solar_apparent_time import solar_day_of_year_for_area, solar_hour_of_day_for_area

from .timer import Timer

from rasters import Raster, RasterGeometry

from .vegetation_conversion.vegetation_conversion import FVC_from_NDVI, LAI_from_NDVI

from .constants import *
from .closure import STIC_closure
from .soil_moisture_initialization import initialize_soil_moisture
from .soil_moisture_iteration import iterate_soil_moisture
from .net_radiation import calculate_net_longwave_radiation
from .initialize_with_solar import initialize_with_solar
from .canopy_air_stream import calculate_canopy_air_stream_vapor_pressure
from .initialize_without_solar import initialize_without_solar
from .iterate_with_solar import iterate_with_solar
from .iterate_without_solar import iterate_without_solar
from .root_zone_initialization import calculate_root_zone_moisture

from .soil_heat_flux import calculate_SEBAL_soil_heat_flux

from .model import STIC_JPL
