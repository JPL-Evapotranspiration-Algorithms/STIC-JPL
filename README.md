# `STIC-JPL` 
## Surface Temperature Initiated Closure (STIC) Evapotranspiration Model Python Implementation

[![CI](https://github.com/JPL-Evapotranspiration-Algorithms/STIC/actions/workflows/ci.yml/badge.svg)](https://github.com/JPL-Evapotranspiration-Algorithms/STIC/actions/workflows/ci.yml)

This repository contains the python implementation for the Surface Temperature Initiated Closure (STIC) evapotranspiration model, used by the ECOsystem Spaceborne Thermal Radiometer Experiment on Space Station (ECOSTRESS) and Surface Biology and Geology (SBG) missions.

This software package is a Python implementation of the Surface Temperature Initiated Closure (STIC) version 1.3 model designed to implement LST to solve the aerodynamic temperature, which is critical for ET modeling. The original STIC model produced by Kaniska Mallick (Mallick et al. 2015; 2018; 2022) was re-implemented from MATLAB code to Python by Gregory Halverson and Madeleine Pascolini-Campbell. The software was developed under a research grant by the NASA Research Opportunities in Space and Earth Sciences (ROSES) program. It is intended for use by the Hyperspectral Thermal Emission Spectrometer (HyTES), MODIS/ASTER (MASTER) Airborne Simulator, Ecosystem Spaceborne Thermal Radiometer Experiment on Space Station (ECOSTRESS) mission and as a precursor for the Surface Biology and Geology (SBG) mission. 

The software was developed as part of a research grant by the NASA Research Opportunities in Space and Earth Sciences (ROSES) program. It was designed for use by the Hyperspectral Thermal Emission Spectrometer (HyTES), MODIS/ASTER (MASTER) Airborne Simulator, Ecosystem Spaceborne Thermal Radiometer Experiment on Space Station (ECOSTRESS) mission as a precursor for the Surface Biology and Geology (SBG) mission. It may also be useful for general remote sensing and GIS projects in Python. This package can be utilized for remote sensing research in Jupyter notebooks and deployed for operations in data processing pipelines.
 
The software is being released according to the SPD-41 open-science requirements of NASA-funded ROSES projects.

Gregory H. Halverson (they/them)<br>
[gregory.h.halverson@jpl.nasa.gov](mailto:gregory.h.halverson@jpl.nasa.gov)<br>
Lead developer<br>
NASA Jet Propulsion Laboratory 329G

Kaniska Mallick (he/him)<br>
[kaniska.mallick@list.lu](mailto:kaniska.mallick@list.lu)<br>
Algorithm inventor<br>
Luxembourg Institute of Science and Technology

Tian Hu (he/him)<br>
[tian.hu@list.lu](mailto:tian.hu@list.lu)<br>
Algorithm developer<br>
Luxembourg Institute of Science and Technology

Madeleine Pascolini-Campbell (she/her)<br>
[madeleine.a.pascolini-campbell@jpl.nasa.gov](mailto:madeleine.a.pascolini-campbell@jpl.nasa.gov)<br>
Algorithm developer<br>
NASA Jet Propulsion Laboratory 329F

Claire Villanueva-Weeks (she/her)<br>
[claire.s.villanueva-weeks@jpl.nasa.gov](mailto:claire.s.villanueva-weeks@jpl.nasa.gov)<br>
Code maintenance<br>
NASA Jet Propulsion Laboratory 329G

## Installation

Use the pip package manager to install the `STIC-JPL` PyPi package with dashes in the name.

```
pip install STIC-JPL
```

## Usage

Import the `STIC_JPL` function from the `STIC_JPL` package with underscores in the name.

```
from STIC_JPL import STIC_JPL
```

See the [ECOSTRESS example](ECOSTRESS%20Example.ipynb) notebook for usage.

See the [STIC sensitivity](STIC%20Sensitivity.ipynb) notebook for sensitivity analysis.

## References
 
Mallick, K., Boegh, E., Trebs, I., Alfieri, J. G., Kustas, W. P., Prueger, J. H., ... & Jarvis, A. J. (2015). Reintroducing radiometric surface temperature into the P enman‐M onteith formulation. Water Resources Research, 51(8), 6214-6243. https://doi.org/10.1002/2014WR016106 

Mallick, K., Toivonen, E., Trebs, I., Boegh, E., Cleverly, J., Eamus, D., ... & Garcia, M. (2018). Bridging Thermal Infrared Sensing and Physically‐Based Evapotranspiration Modeling: From Theoretical Implementation to Validation Across an Aridity Gradient in Australian Ecosystems. Water Resources Research, 54(5), 3409-3435. https://doi.org/10.1029/2017WR021357

Mallick, K., Baldocchi, D., Jarvis, A., Hu, T., Trebs, I., Sulis, M., ... & Kustas, W. P. (2022). Insights Into the Aerodynamic Versus Radiometric Surface Temperature Debate in Thermal‐Based Evaporation Modeling. Geophysical Research Letters, 49(15), e2021GL097568. https://doi.org/10.1029/2021GL097568


