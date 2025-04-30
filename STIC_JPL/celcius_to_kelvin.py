from typing import Union
import numpy as np
from rasters import Raster

def celcius_to_kelvin(T_C: Union[Raster, np.ndarray]) -> Union[Raster, np.ndarray]:
    """
    convert temperature in celsius to kelvin.
    :param T_C: temperature in celsius
    :return: temperature in kelvin
    """
    return T_C + 273.15
