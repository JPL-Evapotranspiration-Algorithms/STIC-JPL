from typing import Union

import numpy as np

import rasters as rt

from rasters import Raster


def calculate_santanello_soil_heat_flux(
        Rn_Wm2: Union[Raster, np.ndarray],
        seconds_of_day: Union[Raster, np.ndarray],
        M: Union[Raster, np.ndarray]
    ) -> Union[Raster, np.ndarray]:
    cg_min = 0.05
    cg_max = 0.35
    tg_min = 74000.0
    tg_max = 100000.0
    solar_noon_seconds = 12.0 * 60.0 * 60.0

    time_offset_seconds = solar_noon_seconds - seconds_of_day

    cg = (1 - M) * cg_max + M * cg_min
    tg = (1 - M) * tg_max + M * tg_min

    G_Wm2 = Rn_Wm2 * cg * np.cos(2 * np.pi * (time_offset_seconds + 10800.0) / tg)
    G_Wm2 = rt.where(Rn_Wm2 < 0, -G_Wm2, G_Wm2)

    return G_Wm2