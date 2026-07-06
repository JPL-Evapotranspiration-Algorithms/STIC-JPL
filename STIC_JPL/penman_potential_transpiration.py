from typing import Union

import numpy as np

from rasters import Raster

from .constants import GAMMA_HPA


def penman_potential_transpiration(
        delta_hPa: Union[Raster, np.ndarray],
        phi_Wm2: Union[Raster, np.ndarray],
        rho_kgm3: Union[Raster, np.ndarray],
        Cp_Jkg: Union[Raster, np.ndarray],
        gB_ms: Union[Raster, np.ndarray],
        VPD_hPa: Union[Raster, np.ndarray],
        gamma_hPa: Union[Raster, np.ndarray, float] = GAMMA_HPA,
        SM: Union[Raster, np.ndarray, float] = 1.0,
        gB_by_gS: Union[Raster, np.ndarray, float] = 0.0
    ) -> Union[Raster, np.ndarray]:
    """
    Compute potential transpiration as latent heat flux using a Penman-Monteith-type expression.

    Equation
    --------
    PT = (delta * phi + rho * Cp * gB * VPD) / (delta + gamma * (1 + SM * gB_by_gS))

    where:
    - delta [hPa K^-1] is the slope of saturation vapor pressure curve,
    - phi [W m^-2] is available energy (Rn - G),
    - rho [kg m^-3] is air density,
    - Cp [J kg^-1 K^-1] is specific heat of air,
    - gB [m s^-1] is boundary-layer conductance,
    - VPD [hPa] is vapor pressure deficit,
    - gamma [hPa K^-1] is psychrometric constant,
    - SM [-] is moisture availability state used by STIC,
    - gB_by_gS [-] is conductance ratio gB / gS.

    Notes
    -----
    This form is the STIC implementation's potential-transpiration diagnostic,
    derived from Penman-Monteith conductance structure with moisture-state weighting.

    References
    ----------
    - Penman, H. L. (1948). Natural evaporation from open water, bare soil and grass.
      Proceedings of the Royal Society A, 193, 120-145. https://doi.org/10.1098/rspa.1948.0037
    - Monteith, J. L. (1965). Evaporation and environment.
      Symposia of the Society for Experimental Biology, 19, 205-234.
    - Mallick, K. et al. (2015). Reintroducing radiometric surface temperature into the Penman-Monteith
      formulation. Water Resources Research, 51, 6214-6243. https://doi.org/10.1002/2014WR016106

    Args:
        delta_hPa: Slope of saturation vapor pressure curve [hPa K^-1].
        phi_Wm2: Available energy [W m^-2].
        rho_kgm3: Air density [kg m^-3].
        Cp_Jkg: Specific heat at constant pressure [J kg^-1 K^-1].
        gB_ms: Boundary-layer conductance [m s^-1].
        VPD_hPa: Vapor pressure deficit [hPa].
        gamma_hPa: Psychrometric constant [hPa K^-1].
        SM: Moisture availability state [-].
        gB_by_gS: Conductance ratio gB/gS [-].

    Returns:
        Potential transpiration latent heat flux [W m^-2].
    """
    return (delta_hPa * phi_Wm2 + rho_kgm3 * Cp_Jkg * gB_ms * VPD_hPa) / (delta_hPa + gamma_hPa * (1 + SM * gB_by_gS))