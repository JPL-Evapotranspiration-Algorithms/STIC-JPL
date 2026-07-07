from typing import Union

import numpy as np

from rasters import Raster

from .constants import GAMMA_HPA, PT_ALPHA


def priestley_taylor_potential_evaporation(
        delta_hPa: Union[Raster, np.ndarray],
        energy_Wm2: Union[Raster, np.ndarray],
        alpha: Union[Raster, np.ndarray, float] = PT_ALPHA,
        gamma_hPa: Union[Raster, np.ndarray, float] = GAMMA_HPA
    ) -> Union[Raster, np.ndarray]:
        """
        Compute Priestley-Taylor potential evaporation as latent heat flux.

        Equation
        --------
        Ep_PT = alpha * delta * A / (delta + gamma)

        where:
        - alpha [-] is the Priestley-Taylor coefficient,
        - delta [hPa K^-1] is the slope of saturation vapor pressure curve,
        - A [W m^-2] is the available energy term supplied by the caller,
        - gamma [hPa K^-1] is psychrometric constant.

        References
        ----------
        - Priestley, C. H. B., & Taylor, R. J. (1972). On the assessment of surface heat flux and
            evaporation using large-scale parameters. Monthly Weather Review, 100, 81-92.
            https://doi.org/10.1175/1520-0493(1972)100<0081:OTAOSH>2.3.CO;2
        - Brutsaert, W., & Stricker, H. (1979). An advection-aridity approach to estimate actual regional
            evapotranspiration. Water Resources Research, 15(2), 443-450. https://doi.org/10.1029/WR015i002p00443
        - Mallick, K. et al. (2014). A Surface Temperature Initiated Closure (STIC) for
            surface energy balance fluxes. Remote Sensing of Environment, 141, 243-261.
            Priestley-Taylor potential evaporation form shown in Eq. (12).
        - Mallick, K. et al. (2015). Reintroducing radiometric surface temperature into the Penman-Monteith
            formulation. Water Resources Research, 51, 6214-6243. https://doi.org/10.1002/2014WR016106

        Args:
                delta_hPa: Slope of saturation vapor pressure curve [hPa K^-1].
                energy_Wm2: Available energy proxy [W m^-2] (for this codebase typically Rn or phi).
                alpha: Priestley-Taylor coefficient [-].
                gamma_hPa: Psychrometric constant [hPa K^-1].

        Returns:
                Priestley-Taylor potential evaporation latent heat flux [W m^-2], preserving input shape.
        """
        return (alpha * delta_hPa * energy_Wm2) / (delta_hPa + gamma_hPa)
