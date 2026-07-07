from typing import Union, Tuple
import numpy as np

import rasters as rt
from rasters import Raster

from .constants import *
from .priestley_taylor_potential_evaporation import priestley_taylor_potential_evaporation

def initialize_without_solar(
            ST_C: Union[Raster, np.ndarray],  # Surface temperature in Celsius
            Ta_C: Union[Raster, np.ndarray],  # Air temperature in Celsius
            dTS: Union[Raster, np.ndarray],  # Temperature difference between surface and air in Celsius
            Td_C: Union[Raster, np.ndarray],  # Dewpoint temperature in Celsius
            Ea_hPa: Union[Raster, np.ndarray],  # Actual vapor pressure in hPa
            Estar_hPa: Union[Raster, np.ndarray],  # Saturation vapor pressure at surface temperature in hPa
            SVP_hPa: Union[Raster, np.ndarray],  # Saturation vapor pressure at the surface in hPa
            delta_hPa: Union[Raster, np.ndarray],  # Slope of the saturation vapor pressure-temperature curve in hPa
            phi_Wm2: Union[Raster, np.ndarray],  # Available energy in W/m2
            gamma_hPa: Union[Raster, np.ndarray, float] = GAMMA_HPA,  # Psychrometric constant in hPa/°C
            alpha: Union[Raster, np.ndarray, float] = PT_ALPHA  # Priestley-Taylor alpha
      ) -> Tuple[Union[Raster, np.ndarray]]:
      """
      Initializes the variables related to moisture and vapor pressure without considering solar radiation.

      Source equations:
      - Mallick et al. (2014), Eq. (2): e_s = e_a + M (e_s* - e_a)
      - Mallick et al. (2014), Eq. (5): moisture availability M from slope/temperature ratios
      - Mallick et al. (2014), Eq. (8): surface dewpoint temperature T_sd
      - Mallick et al. (2014), Eq. (12): Priestley-Taylor potential evaporation

      Parameters:
          ST_C (np.ndarray): Surface temperature in Celsius.
          Ta_C (np.ndarray): Air temperature in Celsius.
          dTS (np.ndarray): Temperature difference between surface and air in Celsius.
          Td_C (np.ndarray): Dewpoint temperature in Celsius.
          Ea_hPa (np.ndarray): Actual vapor pressure in hPa.
          Estar (np.ndarray): saturation vapor pressure at surface temperature (hPa/K)
          SVP_hPa (np.ndarray): Saturation vapor pressure at the surface in hPa.
          delta (np.ndarray): Slope of the saturation vapor pressure-temperature curve in hPa/K.
          phi (np.ndarray): available energy in W/m2.

      Returns:
          Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
              - SM (np.ndarray): Soil moisture.
              - SMrz (np.ndarray): Rootzone moisture.
              - s1 (np.ndarray): Slope of saturation vapor pressure and dewpoint temperature.
              - s3 (np.ndarray): Slope of saturation vapor pressure and temperature.
              - s33 (np.ndarray): Slope of saturation vapor pressure and surface temperature.
              - s44 (np.ndarray): s44    
              - Ms (np.ndarray): Surface moisture.
              - Tsd_C (np.ndarray): Surface dewpoint temperature in Celsius.
              - Es (np.ndarray): Surface vapor pressure in hPa.
              - Ds (np.ndarray): Vapor pressure deficit at the surface.
      """
      s33 = (45.03 + 3.014 * ST_C + 0.05345 * ST_C ** 2 + 0.00224 * ST_C ** 3) * 1e-2  # hpa/K
      s1 = (45.03 + 3.014 * Td_C + 0.05345 * Td_C ** 2 + 0.00224 * Td_C ** 3) * 1e-2  # hpa/K
      
      # Surface dewpoint temperature from STIC derivation (Mallick et al., 2014, Eq. 8).
      Tsd_C = ((Estar_hPa - Ea_hPa) - (s33 * ST_C) + (s1 * Td_C)) / (s1 - s33)
      
      # slope of saturation vapor pressure and temperature
      s3 = rt.where((dTS > -20) & (dTS < 5), (Estar_hPa - Ea_hPa) / (ST_C - Td_C),
            (45.03 + 3.014 * ST_C + 0.05345 * ST_C ** 2 + 0.00224 * ST_C ** 3) * 1e-2)  # hpa/K
      
      # Surface moisture availability M (Mallick et al., 2014, Eq. 5).
      Ms = (s1 / s3) * ((Tsd_C - Td_C) / (ST_C - Td_C))
      Ms = rt.clip(rt.where((dTS < 0) & (Ms < 0) & (phi_Wm2 < 0), np.abs(Ms), Ms), 0, 1)

      # Rootzone Moisture (Mrz)
      s44 = (SVP_hPa - Ea_hPa) / (Ta_C - Td_C)
      SMrz = (gamma_hPa * s1 * (Tsd_C - Td_C)) / (
            delta_hPa * s3 * (ST_C - Td_C) + gamma_hPa * s44 * (Ta_C - Td_C) - delta_hPa * s1 * (
            Tsd_C - Td_C))  # rootzone wetness
      SMrz = rt.clip(rt.where((dTS < 0) & (SMrz < 0) & (phi_Wm2 < 0), np.abs(SMrz), SMrz), 0, 1)

      # now the limits of both Ms and Mrz are consistent
      # combine M to account for Hysteresis and initial estimation of surface vapor pressure
      # Potential evaporation (Priestley-Taylor eqn.)
      Ep_PT = priestley_taylor_potential_evaporation(
            delta_hPa=delta_hPa,
            energy_Wm2=phi_Wm2,
            alpha=alpha,
            gamma_hPa=gamma_hPa
      )
      # Vapor pressure at the evaporating front using STIC mixing relation (Mallick et al., 2014, Eq. 2).
      Es = rt.where((Ep_PT > phi_Wm2) & (dTS > 0) & (Td_C <= 0), Ea_hPa + SMrz * (Estar_hPa - Ea_hPa),
            Ea_hPa + Ms * (Estar_hPa - Ea_hPa))

      # soil moisture
      SM = rt.where((Ep_PT > phi_Wm2) & (dTS > 0) & (Td_C <= 0), SMrz, Ms)

      # hysteresis logic
      # vapor pressure deficit at surface (Ds is later replaced by D0)
      Ds = (Estar_hPa - Es)

      return SM, SMrz, s1, s3, s33, s44, Ms, Tsd_C, Es, Ds
