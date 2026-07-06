import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname("__file__"), '..')))

from typing import Callable
from os import makedirs
from os.path import join
import numpy as np
import pandas as pd
from verma_net_radiation import verma_net_radiation_table
from BESS_JPL import process_BESS_table, load_ECOv002_calval_BESS_inputs
from STIC_JPL import process_STIC_table, load_ECOv002_calval_STIC_inputs
from monte_carlo_sensitivity import perturbed_run, sensitivity_analysis, divide_absolute_by_unperturbed
import matplotlib.pyplot as plt
from scipy.stats import mstats
import seaborn as sns
from matplotlib.ticker import FuncFormatter

normalization_function = divide_absolute_by_unperturbed

input_df = load_ECOv002_calval_STIC_inputs()
BESS_input_df = load_ECOv002_calval_BESS_inputs()

for column in BESS_input_df.columns:
    if column not in input_df.columns:
        input_df[column] = BESS_input_df[column]

input_df = input_df[input_df.ST_C <= 50]
input_df

def with_SWin(input_df: pd.DataFrame) -> pd.DataFrame:
    return process_STIC_table(
        # input_df=verma_net_radiation_table(input_df), 
        input_df=process_BESS_table(
            input_df=input_df,
            offline_mode=True
        ),
        constrain_negative_LE=False,
        supply_SWin=True,
        upscale_to_daylight=True,
        offline_mode=True
    )

def without_SWin(input_df: pd.DataFrame) -> pd.DataFrame:
    return process_STIC_table(
        # input_df=verma_net_radiation_table(input_df), 
        input_df=process_BESS_table(
            input_df=input_df,
            offline_mode=True
        ),
        constrain_negative_LE=False,
        supply_SWin=False,
        upscale_to_daylight=True,
        offline_mode=True
    )

processed_without_SWin = without_SWin(input_df)
print(processed_without_SWin)
