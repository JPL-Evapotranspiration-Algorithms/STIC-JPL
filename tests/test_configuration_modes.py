import numpy as np
import pandas as pd
import importlib

from STIC_JPL.model import _resolve_mode_defaults
from STIC_JPL.process_STIC_table import process_STIC_table
from santanello_soil_heat_flux import santanello_soil_heat_flux as calculate_santanello_soil_heat_flux


def test_resolve_mode_defaults_ecov002():
    defaults = _resolve_mode_defaults("ECOv002")

    assert defaults["configuration"] == "ECOv002"
    assert defaults["LE_convergence_target"] == 1.0
    assert defaults["max_iterations"] == 3
    assert defaults["g_method"] == "santanello"
    assert defaults["constrain_negative_LE"] is False
    assert defaults["constrain_PET"] is True
    assert defaults["apply_surface_emissivity_to_LWin"] is True
    assert defaults["run_iterative_convergence"] is True
    assert defaults["use_buck_dewpoint"] is False
    assert defaults["default_use_variable_alpha"] is True
    assert defaults["constrain_LE_to_available_energy"] is True


def test_resolve_mode_defaults_mallick2014():
    defaults = _resolve_mode_defaults("MALLICK2014")

    assert defaults["configuration"] == "MALLICK2014"
    assert defaults["closure_version"] == "2014"
    assert defaults["g_method"] == "sebal"
    assert defaults["constrain_negative_LE"] is False
    assert defaults["constrain_PET"] is False
    assert defaults["apply_surface_emissivity_to_LWin"] is False
    assert defaults["run_iterative_convergence"] is False
    assert defaults["use_buck_dewpoint"] is True
    assert defaults["default_use_variable_alpha"] is False
    assert defaults["constrain_LE_to_available_energy"] is False


def test_resolve_mode_defaults_mallick2015():
    defaults = _resolve_mode_defaults("MALLICK2015")

    assert defaults["configuration"] == "MALLICK2015"
    assert defaults["closure_version"] == "2015"
    assert defaults["g_method"] == "sebal"
    assert defaults["max_iterations"] == 25
    assert defaults["constrain_negative_LE"] is False
    assert defaults["constrain_PET"] is False
    assert defaults["apply_surface_emissivity_to_LWin"] is False
    assert defaults["run_iterative_convergence"] is True
    assert defaults["use_buck_dewpoint"] is True
    assert defaults["default_use_variable_alpha"] is True
    assert defaults["constrain_LE_to_available_energy"] is False


def test_resolve_mode_defaults_mallick2016():
    defaults = _resolve_mode_defaults("MALLICK2016")

    assert defaults["configuration"] == "MALLICK2016"
    assert defaults["closure_version"] == "2015"  # 2016 uses same STIC1.2 closure as 2015
    assert defaults["g_method"] == "sebal"
    assert defaults["max_iterations"] == 25
    assert defaults["constrain_negative_LE"] is False
    assert defaults["constrain_PET"] is False
    assert defaults["apply_surface_emissivity_to_LWin"] is False
    assert defaults["run_iterative_convergence"] is True
    assert defaults["use_buck_dewpoint"] is True
    assert defaults["default_use_variable_alpha"] is True
    assert defaults["constrain_LE_to_available_energy"] is False


def test_resolve_mode_defaults_invalid_falls_back_to_ecov003(caplog):
    caplog.set_level("WARNING")
    defaults = _resolve_mode_defaults("bad_mode")

    assert defaults["configuration"] == "ECOv003"
    assert "falling back to ECOv003" in caplog.text


def test_process_stic_table_forwards_configuration(monkeypatch):
    call_kwargs = {}
    process_module = importlib.import_module("STIC_JPL.process_STIC_table")

    def fake_stic_jpl(**kwargs):
        call_kwargs.update(kwargs)
        return {"LE_Wm2": np.array([1.0])}

    monkeypatch.setattr(process_module, "STIC_JPL", fake_stic_jpl)

    input_df = pd.DataFrame(
        {
            "ST_C": [30.0],
            "EmisWB": [0.97],
            "NDVI": [0.7],
            "albedo": [0.2],
            "Ta_C": [25.0],
            "RH": [0.5],
            "Rn_Wm2": [450.0],
            "Rg": [900.0],
            "G": [80.0],
            "lat": [34.0],
            "lon": [-118.0],
            "time_UTC": [pd.Timestamp("2023-06-15T19:00:00")],
        }
    )

    output_df = process_STIC_table(
        input_df=input_df,
        configuration="ECOv002",
        offline_mode=True,
    )

    assert call_kwargs["configuration"] == "ECOv002"
    assert "LE_Wm2" in output_df.columns


def test_process_stic_table_ecov002_maps_rg_to_swin(monkeypatch):
    call_kwargs = {}
    process_module = importlib.import_module("STIC_JPL.process_STIC_table")

    def fake_stic_jpl(**kwargs):
        call_kwargs.update(kwargs)
        return {"LE_Wm2": np.array([1.0])}

    monkeypatch.setattr(process_module, "STIC_JPL", fake_stic_jpl)

    input_df = pd.DataFrame(
        {
            "ST_C": [30.0],
            "EmisWB": [0.97],
            "NDVI": [0.7],
            "albedo": [0.2],
            "Ta_C": [25.0],
            "RH": [0.5],
            "Rn_Wm2": [450.0],
            "Rg": [900.0],
            "G": [80.0],
            "lat": [34.0],
            "lon": [-118.0],
            "time_UTC": [pd.Timestamp("2023-06-15T19:00:00")],
        }
    )

    process_STIC_table(
        input_df=input_df,
        configuration="ECOv002",
        supply_SWin=False,
        offline_mode=True,
    )

    assert np.allclose(call_kwargs["SWin_Wm2"], input_df["Rg"].to_numpy())


def test_santanello_soil_heat_flux_helper():
    rn_wm2 = np.array([100.0])
    seconds_of_day = np.array([12.0 * 60.0 * 60.0])
    moisture = np.array([0.5])

    g_wm2 = calculate_santanello_soil_heat_flux(
        seconds_of_day=seconds_of_day,
        Rn=rn_wm2,
        SM=moisture,
    )

    cg = (1 - moisture) * 0.35 + moisture * 0.05
    tg = (1 - moisture) * 100000.0 + moisture * 74000.0
    expected = rn_wm2 * cg * np.cos(2 * np.pi * (0 + 10800.0) / tg)

    assert np.allclose(g_wm2, expected)
