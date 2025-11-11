from typing import Union, Callable
import logging
from datetime import datetime, timedelta
from os.path import join, abspath, expanduser
from typing import Dict, List
import numpy as np
import warnings
from verma_net_radiation import daylight_Rn_integration_verma
from sun_angles import calculate_daylight

from pytictoc import TicToc

import colored_logging as cl
from check_distribution import check_distribution
import rasters as rt
from GEOS5FP import GEOS5FP
from solar_apparent_time import calculate_solar_day_of_year, calculate_solar_hour_of_day
from SEBAL_soil_heat_flux import calculate_SEBAL_soil_heat_flux

from rasters import Raster, RasterGeometry

from daylight_evapotranspiration import daylight_ET_from_instantaneous_LE

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
from .FVC_from_NDVI import FVC_from_NDVI
from .LAI_from_NDVI import LAI_from_NDVI
from .celcius_to_kelvin import celcius_to_kelvin

__author__ = 'Kaniska Mallick, Madeleine Pascolini-Campbell, Gregory Halverson'

logger = logging.getLogger(__name__)

def STIC_JPL(
        ST_C: Union[Raster, np.ndarray],
        emissivity: Union[Raster, np.ndarray],
        NDVI: Union[Raster, np.ndarray],
        albedo: Union[Raster, np.ndarray],
        Rn_Wm2: Union[Raster, np.ndarray],
        geometry: RasterGeometry = None,
        time_UTC: datetime = None,
        hour_of_day: np.ndarray = None,
        day_of_year: np.ndarray = None,
        GEOS5FP_connection: GEOS5FP = None,
        Ta_C: Union[Raster, np.ndarray] = None,
        RH: Union[Raster, np.ndarray] = None,
        G_Wm2: Union[Raster, np.ndarray] = None,
        G_method: str = DEFAULT_G_METHOD,
        SM: Union[Raster, np.ndarray] = None,
        SWin_Wm2: Union[Raster, np.ndarray] = None,
        FVC: Union[Raster, np.ndarray] = None,
        LAI: Union[Raster, np.ndarray] = None,
        elevation_m: Union[Raster, np.ndarray] = None,
        delta_hPa: Union[Raster, np.ndarray] = None,
        gamma_hPa: Union[Raster, np.ndarray, float] = GAMMA_HPA,
        rho_kgm3: Union[Raster, np.ndarray] = RHO_KGM3,
        Cp_Jkg: Union[Raster, np.ndarray] = CP_JKG,
        alpha: float = PT_ALPHA,
        LE_convergence_target: float = LE_CONVERGENCE_TARGET_WM2,
        max_iterations: int = MAX_ITERATIONS,
        show_distributions: bool = SHOW_DISTRIBUTIONS,
        use_variable_alpha: bool = USE_VARIABLE_ALPHA,
        upscale_to_daylight: bool = UPSCALE_TO_DAYLIGHT,
        constrain_negative_LE: bool = CONSTRAIN_NEGATIVE_LE,
        resampling: str = RESAMPLING) -> Dict[str, Union[Raster, np.ndarray]]:
    results = {}
    # For daily upscaling
    Rn_daily_Wm2 = None
    EF = None
    LE_daylight_Wm2 = None
    ET_daily_kg = None

    if geometry is None and isinstance(ST_C, Raster):
        geometry = ST_C.geometry

    if GEOS5FP_connection is None:
        GEOS5FP_connection = GEOS5FP()

    if (day_of_year is None or hour_of_day is None) and time_UTC is not None and geometry is not None:
        day_of_year = calculate_solar_day_of_year(time_UTC=time_UTC, geometry=geometry)
        hour_of_day = calculate_solar_hour_of_day(time_UTC=time_UTC, geometry=geometry)

    if time_UTC is None and day_of_year is None and hour_of_day is None:
        raise ValueError("no time given between time_UTC, day_of_year, and hour_of_day")

    seconds_of_day = hour_of_day * 3600.0

    # load air temperature in Celsius if not provided
    if Ta_C is None:
        Ta_C = GEOS5FP_connection.Ta_C(time_UTC=time_UTC, geometry=geometry, resampling=resampling)

    # load relative humidity if not provided
    if RH is None:
        RH = GEOS5FP_connection.RH(time_UTC=time_UTC, geometry=geometry, resampling=resampling)

    # calculate fraction of vegetation cover if it's not given
    if FVC is None:
        FVC = FVC_from_NDVI(NDVI)
    
    # calculate leaf area index if it's not given
    if LAI is None:
        LAI = LAI_from_NDVI(NDVI)

    # saturation air pressure in hPa
    SVP_hPa = 6.13753 * (np.exp((17.27 * Ta_C) / (Ta_C + 237.3)))

    # calculate delta term if it's not given
    if delta_hPa is None:
        # slope of saturation vapor pressure to air temperature (hpa/K)
        delta_hPa = 4098 * SVP_hPa / (Ta_C + 237.3) ** 2

    Ta_K = celcius_to_kelvin(Ta_C)

    # actual vapor pressure at TA (hpa/K)
    Ea_hPa = SVP_hPa * (RH)
    Ea_Pa = Ea_hPa * 100.0

    # vapor pressure deficit (hPa)
    VPD_hPa = SVP_hPa - Ea_hPa

    # swapping in the dew-point calculation from PT-JPL
    Td_C = Ta_C - ((100 - RH * 100) / 5.0)

    # difference between surface and air temperature (Celsius)
    dTS_C = ST_C - Ta_C

    # saturation vapor pressure at surface temperature (hPa/K)
    Estar_hPa = 6.13753 * np.exp((17.27 * ST_C) / (ST_C + 237.3))

    if SWin_Wm2 is None:
        # if G is None and SM is None:
        #     raise ValueError("soil heat flux or soil moisture prior required if solar radiation is not given")

        if G_Wm2 is None:
            G_Wm2 = calculate_SEBAL_soil_heat_flux(
                ST_C=ST_C,
                NDVI=NDVI,
                albedo=albedo,
                Rn=Rn_Wm2,
            )
        
        phi_Wm2 = Rn_Wm2 - G_Wm2

        # initialize without solar radiation
        SM, SMrz, s1, s3, s33, s44, Ms, Tsd_C, Es_hPa, Ds = initialize_without_solar(
            ST_C = ST_C,  # Surface temperature in Celsius
            Ta_C = Ta_C,  # Air temperature in Celsius
            dTS = dTS_C,  # Temperature difference between surface and air in Celsius
            Td_C = Td_C,  # Dewpoint temperature in Celsius
            Ea_hPa = Ea_hPa,  # Actual vapor pressure in hPa
            Estar_hPa = Estar_hPa,  # Saturation vapor pressure at surface temperature (hPa/K)
            SVP_hPa = SVP_hPa,  # Saturation vapor pressure at the surface in hPa
            delta_hPa = delta_hPa,  # Slope of the saturation vapor pressure-temperature curve in hPa/K
            phi_Wm2 = phi_Wm2,  # Available energy in W/m2
            gamma_hPa = gamma_hPa,  # Psychrometric constant in hPa/°C
            alpha = alpha  # Priestley-Taylor alpha
        )
    else:
        SM, SMrz, Ms, s1, s3, Ep_PT, Rnsoil, LWnet_Wm2, G_Wm2, Tsd_C, Ds, Es_hPa, phi_Wm2 = initialize_with_solar(
            seconds_of_day = seconds_of_day,  # time of day in seconds since midnight
            Rg_Wm2 = SWin_Wm2,  # solar radiation (W/m^2)
            Rn_Wm2 = Rn_Wm2,  # net radiation (W/m^2)
            ST_C = ST_C,  # surface temperature (Celsius)
            emissivity = emissivity,  # emissivity of the surface
            Ta_C = Ta_C,  # air temperature (Celsius)
            dTS_C = dTS_C,  # surface air temperature difference (Celsius)
            Td_C = Td_C,  # dew point temperature (Celsius)
            VPD_hPa = VPD_hPa,  # vapor pressure deficit (hPa)
            SVP_hPa = SVP_hPa,  # saturation vapor pressure at given air temperature (hPa)
            Ea_hPa = Ea_hPa,  # actual vapor pressure at air temperature (hPa)
            Estar_hPa = Estar_hPa,  # saturation vapor pressure at surface temperature (hPa)
            delta_hPa = delta_hPa,  # slope of saturation vapor pressure to air temperature (hpa/K)
            NDVI=NDVI,  # normalized difference vegetation index
            FVC = FVC,  # fractional vegetation cover
            LAI = LAI,  # leaf area index
            albedo = albedo,  # albedo of the surface
            gamma_hPa=gamma_hPa,  # psychrometric constant (hPa/°C)
            G_method = DEFAULT_G_METHOD,  # method for calculating soil heat flux
        )
    
    check_distribution(Ms, "Ms")

    # STIC analytical equations (convergence on LE)
    gB_ms, gS_ms, dT_C, EF = STIC_closure(
        delta_hPa=delta_hPa, 
        phi_Wm2=phi_Wm2, 
        Es_hPa=Es_hPa, 
        Ea_hPa=Ea_hPa, 
        Estar_hPa=Estar_hPa, 
        SM=SM,
        gamma_hPa=gamma_hPa,
        rho_kgm3=rho_kgm3,
        Cp_Jkg=Cp_Jkg,
        alpha=alpha
    )
    
    gBB = gB_ms
    gSS = gS_ms
    gBB_by_gSS = rt.where(gSS == 0, 0, gBB / gSS)
    gB_by_gS = rt.where(gS_ms == 0, 0, gB_ms / gS_ms)
    dT_C = dT_C
    T0_C = dT_C + Ta_C
    
    PET_Wm2 = ((delta_hPa * phi_Wm2 + rho_kgm3 * Cp_Jkg * gB_ms * VPD_hPa) / (delta_hPa + gamma_hPa))  # Penman potential evaporation
    
    gR = (4 * SB_SIGMA * (Ta_C + 273) ** 3 * emissivity) / (rho_kgm3 * Cp_Jkg)
    omega = ((delta_hPa / gamma_hPa) + 1) / ((delta_hPa / gamma_hPa) + 1 + gB_by_gS)
    LE_eq = (phi_Wm2 * (delta_hPa / gamma_hPa)) / ((delta_hPa / gamma_hPa) + 1)
    LE_imp = (Cp_Jkg * 0.0289644 / gamma_hPa) * gS_ms * 40 * VPD_hPa
    LE_init = omega * LE_eq + (1 - omega) * LE_imp
    dry = (Ds > VPD_hPa) & (PET_Wm2 > phi_Wm2) & (dTS_C > 0) & (Td_C <= 0)
    omega = rt.where(dry,
                        ((delta_hPa / gamma_hPa) + 1 + gR / gB_ms) / ((delta_hPa / gamma_hPa) + 1 + gB_ms / gS_ms + gR / gS_ms + gR / gB_ms),
                        omega)
    LE_eq = rt.where(dry, (phi_Wm2 * (delta_hPa / gamma_hPa)) / ((delta_hPa / gamma_hPa) + 1 + gR / gB_ms), LE_eq)
    LE_init = rt.where(dry, omega * LE_eq + (1 - omega * LE_imp), LE_init)
    
    # sensible heat flux
    H_Wm2 = ((gamma_hPa * phi_Wm2 * (1 + gB_by_gS) - rho_kgm3 * Cp_Jkg * gB_ms * VPD_hPa) / (delta_hPa + gamma_hPa * (1 + (gB_by_gS))))
    
    LE_Wm2_new = LE_init
    LE_Wm2_change = LE_convergence_target
    LE_Wm2_old = LE_Wm2_new
    LE_transpiration_Wm2 = None
    PT_Wm2 = None
    iteration = 1
    LE_Wm2_max_change = 0

    t = TicToc()
    t.tic()

    while (np.nanmax(LE_Wm2_change) >= LE_convergence_target and iteration <= max_iterations):
        logger.info(f"running STIC iteration {cl.val(iteration)} / {cl.val(max_iterations)}")

        if SWin_Wm2 is None:
            SM, SMrz, Ms, s1, e0, e0star, Tsd_C, D0, alphaN = iterate_without_solar(
                LE = LE_Wm2_new,  # Latent heat flux (W/m^2)
                PET = PET_Wm2,  # Potential evapotranspiration (W/m^2)
                SM = SM,
                ST_C = ST_C,  # Surface temperature (°C)
                Ta_C = Ta_C,  # Air temperature (°C)
                dTS = dTS_C,  # Surface-air temperature difference (°C)
                T0 = T0_C,  # Reference temperature (°C)
                gB = gB_ms,  # Boundary layer conductance (m/s)
                gS = gS_ms,  # Stomatal conductance (m/s)
                Ea_hPa = Ea_hPa,  # Actual vapor pressure (hPa)
                Td_C = Td_C,  # Dew point temperature (°C)
                VPD_hPa = VPD_hPa,  # Vapor pressure deficit (hPa)
                Estar = Estar_hPa,  # Saturation vapor pressure at surface temperature (hPa)
                delta = delta_hPa,  # Slope of the saturation vapor pressure-temperature curve (hPa/°C)
                phi = phi_Wm2,  # available energy (W/m^2)
                Ds = Ds,  # Vapor pressure deficit at source (hPa)
                Es = Es_hPa,  # Saturation vapor pressure (hPa)
                s3 = s3,  # Slope of the saturation vapor pressure and temperature 
                s4 = s44,  # Slope of the saturation vapor pressure and temperature 
                gB_by_gS = gB_by_gS,  # Ratio of boundary layer conductance to stomatal conductance
                gamma_hPa = gamma_hPa,  # Psychrometric constant (hPa/°C)
                rho_kgm3 = rho_kgm3,  # Air density (kg/m^3)
                Cp_Jkg = Cp_Jkg  # Specific heat at constant pressure (J/kg/K)
            )
        else:
            SM, G_Wm2, e0, e0star, D0, alphaN = iterate_with_solar(
                seconds_of_day = seconds_of_day,  # Seconds of the day
                ST_C = ST_C,  # Soil temperature (°C)
                NDVI = NDVI,  # Normalized Difference Vegetation Index
                albedo = albedo,  # Albedo
                gB_ms = gB_ms,  # boundary layer conductance (m/s)
                gS_ms = gS_ms,  # stomatal conductance (m/s)
                LE_Wm2 = LE_Wm2_new,  # latent heat flux (W/m^2)
                Rg_Wm2 = SWin_Wm2,  # Incoming solar radiation (W/m^2)
                Rn_Wm2 = Rn_Wm2,  # Net radiation (W/m^2)
                LWnet_Wm2 = LWnet_Wm2,  # Net longwave radiation (W/m^2)
                Ta_C = Ta_C,  # Air temperature (°C)
                dTS_C = dTS_C,  # Change in soil temperature (°C)
                Td_C = Td_C,  # Dew point temperature (°C)
                Tsd_C = Tsd_C,  # Soil dew point temperature (°C)
                Ea_hPa = Ea_hPa,  # actual vapor pressure (hPa)
                Estar_hPa = Estar_hPa,  # saturation vapor pressure at surface temperature (hPa)
                VPD_hPa = VPD_hPa,  # Vapor pressure deficit (hPa)
                SVP_hPa = SVP_hPa,  # Saturation vapor pressure (hPa)
                delta_hPa = delta_hPa,  # Slope of the saturation vapor pressure-temperature curve (hPa/°C)
                phi_Wm2 = phi_Wm2,  # Net radiation minus soil heat flux (W/m^2)
                Es_hPa = Es_hPa,  # Saturation vapor pressure (hPa)
                s1 = s1,  # Soil moisture parameter
                s3 = s3,  # Soil moisture parameter
                FVC = FVC,  # Fractional canopy cover
                T0_C = T0_C,  # Reference temperature (°C)
                gB_by_gS = gB_by_gS,  # Ratio of boundary layer conductance to stomatal conductance
                gamma_hPa = gamma_hPa,  # Psychrometric constant (hPa/°C)
                rho_kgm3 = rho_kgm3,  # Air density (kg/m^3)
                Cp_Jkg = Cp_Jkg,  # Specific heat at constant pressure (J/kg/K)
                G_method = "santanello"  # Method for calculating soil heat flux
            )

        if use_variable_alpha:
            alpha = alphaN
            logger.info(f"using variable Priestley-Taylor alpha with mean: {cl.val(np.round(np.nanmean(alpha), 3))}")

        # re-estimated conductances and states
        gB_ms, gS_ms, dT_C, EF = STIC_closure(
            delta_hPa=delta_hPa,  # Slope of the saturation vapor pressure-temperature curve (hPa/°C)
            phi_Wm2=phi_Wm2,  # available energy (W/m^2)
            Es_hPa=Es_hPa,  # Vapor pressure at the reference height (hPa)
            Ea_hPa=Ea_hPa,  # Actual vapor pressure (hPa)
            Estar_hPa=Estar_hPa,  # Saturation vapor pressure at the reference height (hPa)
            SM=SM,  # Soil moisture (m³/m³)
            gamma_hPa=gamma_hPa,  # Psychrometric constant (hPa/°C)
            rho_kgm3=rho_kgm3,  # Air density (kg/m³)
            Cp_Jkg=Cp_Jkg,  # Specific heat capacity of air (J/kg/°C)
            alpha=alpha  # Stability correction factor for conductance 
        )

        gB_by_gS = rt.where(gS_ms == 0, 0, gB_ms / gS_ms)
        T0_C = dT_C + Ta_C
        # latent heat flux
        LE_Wm2_new = ((delta_hPa * phi_Wm2 + rho_kgm3 * Cp_Jkg * gB_ms * VPD_hPa) / (delta_hPa + gamma_hPa * (1 + gB_by_gS)))
        LE_Wm2_new = rt.where(LE_Wm2_new > phi_Wm2, phi_Wm2, LE_Wm2_new)

        if constrain_negative_LE:
            LE_Wm2_new = rt.where(LE_Wm2_new < 0, 0, LE_Wm2_new)

        # Sensible Heat Flux
        H_Wm2 = ((gamma_hPa * phi_Wm2 * (1 + gB_by_gS) - rho_kgm3 * Cp_Jkg * gB_ms * VPD_hPa) / (delta_hPa + gamma_hPa * (1 + (gB_by_gS))))
        # potential evaporation (Penman)
        PET_Wm2 = ((delta_hPa * phi_Wm2 + rho_kgm3 * Cp_Jkg * gB_ms * VPD_hPa) / (delta_hPa + gamma_hPa))
        # Potential Transpiration
        PT_Wm2 = (delta_hPa * phi_Wm2 + rho_kgm3 * Cp_Jkg * gB_ms * VPD_hPa) / (delta_hPa + gamma_hPa * (1 + SM * gB_by_gS))  # potential transpiration
        # ET PARTITIONING
        LE_soil_Wm2 = rt.clip(SM * PET_Wm2, 0, None)
        LE_transpiration_Wm2 = rt.clip(LE_Wm2_new - LE_soil_Wm2, 0, None)
        # change in latent heat flux estimate
        LE_Wm2_change = np.abs(LE_Wm2_old - LE_Wm2_new)
        LE_Wm2_new = rt.where(np.isnan(LE_Wm2_new), LE_Wm2_old, LE_Wm2_new)
        LE_Wm2_old = LE_Wm2_new
        LE_Wm2_max_change = np.nanmax(LE_Wm2_change)
        logger.info(
            f"completed STIC iteration {cl.val(iteration)} / {cl.val(max_iterations)} with max LE change: {cl.val(np.round(LE_Wm2_max_change, 3))} ({t.tocvalue()} seconds)")
        
        check_distribution(SM, f"SM_{iteration}")
        check_distribution(G_Wm2, f"G_{iteration}")
        check_distribution(LE_Wm2_new, f"LE_{iteration}")

        if LE_Wm2_max_change <= LE_convergence_target:
            logger.info(f"max LE change {cl.val(np.round(LE_Wm2_max_change, 3))} within convergence target {cl.val(np.round(LE_convergence_target, 3))} with {cl.val(iteration)} iteration{'s' if iteration > 1 else ''}")

        iteration += 1

    iteration -= 1
    results["LE_max_change"] = LE_Wm2_max_change
    results["iteration"] = iteration

    LE_Wm2 = LE_Wm2_new

    results["LE_Wm2"] = LE_Wm2
    results["LE_change"] = LE_Wm2_change
    results["LE_transpiration_Wm2"] = LE_transpiration_Wm2
    results["PT"] = PT_Wm2
    results["PET_Wm2"] = PET_Wm2
    results["G_Wm2"] = G_Wm2

    if isinstance(geometry, RasterGeometry):
        for name, array in results.items():
            try:
                results[name] = Raster(array.reshape(geometry.shape), geometry=geometry)
            except Exception as e:
                pass
        
        results["LE_Wm2"].cmap = ET_COLORMAP
        results["PET_Wm2"].cmap = ET_COLORMAP

    warnings.resetwarnings()

    # --- Daily Upscaling (if requested) ---
    if upscale_to_daylight and time_UTC is not None:
        logger.info("started daily ET upscaling (STIC-JPL)")
        t_et = TicToc()
        t_et.tic()

                # Use new upscaling function from daylight_evapotranspiration
        daylight_results = daylight_ET_from_instantaneous_LE(
            LE_instantaneous_Wm2=LE_Wm2,
            Rn_instantaneous_Wm2=Rn_Wm2,
            G_instantaneous_Wm2=G_Wm2,
            day_of_year=day_of_year,
            time_UTC=time_UTC,
            geometry=geometry
        )
        # Add all returned daylight results to output
        results.update(daylight_results)

        elapsed_et = t_et.tocvalue()
        logger.info(f"completed daily ET upscaling (elapsed: {elapsed_et:.2f} seconds)")

    return results
