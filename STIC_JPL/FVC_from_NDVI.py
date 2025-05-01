from typing import Union
import numpy as np
import rasters as rt
from rasters import Raster

KPAR = 0.5
MIN_FIPAR = 0.0
MAX_FIPAR = 1.0
MIN_LAI = 0.0
MAX_LAI = 10.0

def FVC_from_NDVI(NDVI: Union[Raster, np.ndarray]) -> Union[Raster, np.ndarray]:
    """
    Estimate Fractional Vegetation Cover (FVC) from Normalized Difference Vegetation Index (NDVI)
    using a scaled NDVI approach.

    This method linearly scales NDVI values between two endmembers:
        - NDVIs: NDVI value for bare soil (typically ~0.04 ± 0.03)
        - NDVIv: NDVI value for full vegetation (typically ~0.52 ± 0.03)

    The resulting Fractional Vegetation Cover (FVC) is calculated as:

        FVC = clip((NDVI - NDVIs) / (NDVIv - NDVIs), 0.0, 1.0)

    This approach is based on the assumption that NDVI increases linearly with vegetation cover
    between these two extremes, and is well-supported in the literature.

    References:
        Carlson, T. N., & Ripley, D. A. (1997). On the relation between NDVI, fractional vegetation cover,
        and leaf area index. Remote Sensing of Environment, 62(3), 241–252.
        https://doi.org/10.1016/S0034-4257(97)00104-1

        Gutman, G., & Ignatov, A. (1998). The derivation of the green vegetation fraction from NOAA/AVHRR
        data for use in numerical weather prediction models. International Journal of Remote Sensing,
        19(8), 1533–1543. https://doi.org/10.1080/014311698215333

    Parameters:
        NDVI (Union[Raster, np.ndarray]): Input NDVI data.

    Returns:
        Union[Raster, np.ndarray]: Estimated Fractional Vegetation Cover (FVC).
    """
    NDVIv = 0.52  # NDVI for fully vegetated pixel
    NDVIs = 0.04  # NDVI for bare soil pixel

    # Scale NDVI to FVC using a linear model and clip to [0, 1]
    FVC = rt.clip((NDVI - NDVIs) / (NDVIv - NDVIs), 0.0, 1.0)

    return FVC