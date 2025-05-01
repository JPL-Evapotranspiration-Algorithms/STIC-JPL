from typing import Union
import numpy as np
import rasters as rt
from rasters import Raster

# Constants
KPAR = 0.5  # Extinction coefficient for PAR, assumed average for broadleaf canopies (Weiss & Baret, 2010)
MIN_FIPAR = 0.0
MAX_FIPAR = 1.0
MIN_LAI = 0.0
MAX_LAI = 10.0

def LAI_from_NDVI(
        NDVI: Union[Raster, np.ndarray],
        min_fIPAR: float = MIN_FIPAR,
        max_fIPAR: float = MAX_FIPAR,
        min_LAI: float = MIN_LAI,
        max_LAI: float = MAX_LAI) -> Union[Raster, np.ndarray]:
    """
    Estimate Leaf Area Index (LAI) from NDVI using a simplified two-step empirical model.

    This method first approximates the fraction of absorbed photosynthetically active radiation (fIPAR)
    from NDVI, and then estimates LAI using the Beer–Lambert Law. The extinction coefficient for PAR (KPAR)
    is assumed to be 0.5, which is typical for broadleaf canopies under diffuse light conditions.

    Steps:
    1. fIPAR ≈ NDVI - 0.05 (empirical offset to account for soil background and sensor noise)
       - Based on observed relationships in Myneni & Williams (1994)
    2. LAI = -ln(1 - fIPAR) / KPAR (Beer–Lambert Law)
       - From Sellers (1985)

    All outputs are clipped to user-defined minimum and maximum thresholds to ensure physical realism.

    Parameters:
        NDVI (Union[Raster, np.ndarray]): Input NDVI data.
        min_fIPAR (float): Minimum fIPAR value for clipping (default 0.0).
        max_fIPAR (float): Maximum fIPAR value for clipping (default 1.0).
        min_LAI (float): Minimum LAI value for clipping (default 0.0).
        max_LAI (float): Maximum LAI value for clipping (default 10.0).

    Returns:
        Union[Raster, np.ndarray]: Estimated LAI values.

    References:
        - Sellers, P. J. (1985). Canopy reflectance, photosynthesis and transpiration. 
          *International Journal of Remote Sensing*, 6(8), 1335–1372.
        - Myneni, R. B., & Williams, D. L. (1994). On the relationship between FAPAR and NDVI. 
          *Remote Sensing of Environment*, 49(3), 200–211.
        - Weiss, M., & Baret, F. (2010). CAN-EYE V6.1 User Manual. INRA-CSE.

    """
    # Empirical conversion from NDVI to fIPAR (adjusted for background noise)
    fIPAR = rt.clip(NDVI - 0.05, min_fIPAR, max_fIPAR)

    # Avoid division by zero or log of 0 by masking zero fIPAR values
    fIPAR = np.where(fIPAR == 0, np.nan, fIPAR)

    # Apply Beer–Lambert law to estimate LAI
    LAI = rt.clip(-np.log(1 - fIPAR) * (1 / KPAR), min_LAI, max_LAI)

    return LAI