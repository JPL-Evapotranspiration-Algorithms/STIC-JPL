from typing import Union
import numpy as np

from rasters import Raster

from .constants import SB_SIGMA

def calculate_net_longwave_radiation(
        Ta_C: Union[Raster, np.ndarray], 
        Ea_hPa: Union[Raster, np.ndarray], 
        ST_C: Union[Raster, np.ndarray], 
        emissivity: Union[Raster, np.ndarray], 
        sigma: float = SB_SIGMA,
        apply_surface_emissivity_to_LWin: bool = False) -> Union[Raster, np.ndarray]:
    """
    Calculates the net longwave radiation at the Earth's surface under clear-sky conditions.

    Net longwave radiation is the energy balance between incoming downward longwave radiation 
    from the atmosphere and the outgoing upward longwave radiation from the surface 
    (which includes both actively emitted radiation and reflected incoming radiation).

    Parameters:
    -----------
    Ta_C : Union[Raster, np.ndarray]
        Air temperature at reference height (typically 2 meters) in degrees Celsius.
    Ea_hPa : Union[Raster, np.ndarray]
        Actual vapor pressure of the air in hectopascals (hPa) or millibars (mb).
    ST_C : Union[Raster, np.ndarray]
        Land surface temperature (LST) or radiometric surface temperature in degrees Celsius.
    emissivity : Union[Raster, np.ndarray]
        Broadband surface emissivity (dimensionless, typically between 0.9 and 1.0).
    sigma : float, optional
        Stefan-Boltzmann constant (W * m^-2 * K^-4). Defaults to SB_SIGMA.
    apply_surface_emissivity_to_LWin : bool, optional
        If True, applies Kirchhoff's Law of thermal radiation to account for the fraction 
        of incoming longwave radiation that is reflected by the surface. If False, assumes 
        the surface acts as a perfect blackbody for incoming radiation (absorbing 100%). 
        Defaults to False.

    Returns:
    --------
    LWnet : Union[Raster, np.ndarray]
        Net longwave radiation at the surface in Watts per square meter (W/m^2). 
        Positive values indicate a net energy gain by the surface, negative indicates a loss.

    References:
    -----------
    - Brutsaert, W. (1975). On a derivable formula for long-wave radiation from clear skies. 
      Water Resources Research, 11(5), 742-744. (Used for clear-sky air emissivity equation).
    - Kirchhoff's Law of Thermal Radiation: States that for an arbitrary body in thermal 
      equilibrium, its emissivity is equal to its absorptivity. Thus, a surface reflects 
      (1 - emissivity) of incoming atmospheric longwave radiation.
    - Stefan-Boltzmann Law: Describes the thermal power radiated by a black body as 
      proportional to the fourth power of its absolute temperature (E = sigma * T^4).
    """
    # 1. Atmospheric Emissivity (etaa)
    # Calculated using the empirical model proposed by Brutsaert (1975) for clear skies.
    # Note: Temperatures must be converted from Celsius to Kelvin (+ 273.15).
    etaa = 1.24 * (Ea_hPa / (Ta_C + 273.15)) ** (1.0 / 7.0)  
    
    # 2. Incoming Longwave Radiation (LWin)
    # The thermal radiation emitted downward by the atmosphere. 
    # Derived via the Stefan-Boltzmann Law using atmospheric properties (Ta and etaa).
    LWin = sigma * etaa * (Ta_C + 273.15) ** 4
    
    # 3. Outgoing Emitted Longwave Radiation (LWout_emitted)
    # The thermal radiation actively emitted by the Earth's surface.
    # Derived via the Stefan-Boltzmann Law using surface properties (ST_C and surface emissivity).
    LWout_emitted = sigma * emissivity * (ST_C + 273.15) ** 4
    
    # 4. Net Longwave Radiation (LWnet)
    if apply_surface_emissivity_to_LWin:
        # Applies Kirchhoff's Law: 
        # The surface only absorbs a fraction of LWin equal to its 'emissivity'.
        # The rest (1 - emissivity) is physically reflected back into the atmosphere.
        # LWnet = (Absorbed LWin) - (Emitted LWout)
        LWnet = (emissivity * LWin) - LWout_emitted
    else:
        # Assumes the surface is a blackbody regarding incoming radiation (absorptivity = 1.0).
        # It absorbs 100% of LWin and reflects nothing.
        # LWnet = (Total LWin) - (Emitted LWout)
        LWnet = LWin - LWout_emitted

    return LWnet
