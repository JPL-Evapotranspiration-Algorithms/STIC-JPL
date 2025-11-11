from matplotlib.colors import LinearSegmentedColormap

RESAMPLING = "cubic"
CONSTRAIN_NEGATIVE_LE = False
SUPPLY_SWIN = False
UPSCALE_TO_DAYLIGHT = False

RHO_KGM3 = 1.2  # Air density (kg m-3)
CP_JKG = 1013  # Specific heat of air at constant pressure (J/kg/K)
GAMMA_HPA = 0.67  # Psychrometric constant (hpa/K)
PT_ALPHA = 1.26  # Priestley-Taylor coefficient
SB_SIGMA = 5.67e-8  # Stefann Boltzmann constant
DEFAULT_G_METHOD = "santanello"
LE_CONVERGENCE_TARGET_WM2 = 2.0
MAX_ITERATIONS = 30
USE_VARIABLE_ALPHA = True
SHOW_DISTRIBUTIONS = True

ET_COLORMAP = LinearSegmentedColormap.from_list("ET", [
    "#f6e8c3",
    "#d8b365",
    "#99974a",
    "#53792d",
    "#6bdfd2",
    "#1839c5"
])
