# A Surface Temperature Initiated Closure (STIC) for surface energy balance fluxes

**Remote Sensing of Environment 141 (2014) 243-261** *Contents lists available at ScienceDirect* *Elsevier* *Journal homepage: [www.elsevier.com/locate/rse](http://www.elsevier.com/locate/rse)*

**Authors:** Kaniska Mallick $^{a,*}$, Andrew J. Jarvis $^b$, Eva Boegh $^c$, Joshua B. Fisher $^a$, Darren T. Drewry $^a$, Kevin P. Tu $^d$, Simon J. Hook $^a$, Glynn Hulley $^a$, Jonas Ardö $^e$, Jason Beringer $^f$, Altaf Arain $^g$, Dev Niyogi $^h$

$^a$ *Jet Propulsion Laboratory, California Institute of Technology, 4800 Oak Grove Drive, Pasadena, CA 91109, United States* $^b$ *Lancaster Environment Centre, Lancaster University, LA1 4YQ, United Kingdom* $^c$ *Department of Environmental Social and Spatial Change, Roskilde University, Denmark* $^d$ *Pioneer Hi-Bred International, Woodland, CA, United States* $^e$ *Department of Physical Geography and Ecosystems Science, Lund University, Sweden* $^f$ *School of Geography and Environmental Science, Monash University, Clayton, Victoria 3800, Australia* $^g$ *School of Geography and Earth Sciences, McMaster University, Ontario, Canada* $^h$ *Department of Agronomy and Earth Atmospheric Planetary Sciences, Purdue University, West Lafayette, IN, United States*

**Article History:** * Received 19 November 2012  
* Received in revised form 11 September 2013  
* Accepted 25 October 2013  

**Keywords:** Surface energy balance, Penman-Monteith equation, Advection-aridity hypothesis, Boundary layer conductance, Surface conductance, MODIS, Land surface temperature, FLUXNET, Evapotranspiration

---

## Abstract
The use of Penman-Monteith (PM) equation in thermal remote sensing based surface energy balance modeling is not prevalent due to the unavailability of any direct method to integrate thermal data into the PM equation and due to the lack of physical models expressing the surface (or stomatal) and boundary layer conductances ($g_S$ and $g_B$) as a function of surface temperature. Here we demonstrate a new method that physically integrates the radiometric surface temperature ($T_S$) into the PM equation for estimating the terrestrial surface energy balance fluxes (sensible heat, $H$ and latent heat, $\lambda E$). The method combines satellite $T_S$ data with standard energy balance closure models in order to derive a hybrid closure that does not require the specification of surface to atmosphere conductance terms. We call this the Surface Temperature Initiated Closure (STIC), which is formed by the simultaneous solution of four state equations. Taking advantage of the psychrometric relationship between temperature and vapor pressure, the present method also estimates the near surface moisture availability ($M$) from $T_S$, air temperature ($T_A$) and relative humidity ($R_H$), thereby being capable of decomposing $\lambda E$ into evaporation ($\lambda E_E$) and transpiration ($\lambda E_T$). 

STIC is driven with $T_S$, $T_A$, $R_H$, net radiation ($R_N$), and ground heat flux ($G$). $T_S$ measurements from both MODIS Terra (MOD11A2) and Aqua (MYD11A2) were used in conjunction with FLUXNET $R_N$, $G$, $T_A$, $R_H$, $\lambda E$ and $H$ measurements corresponding to the MODIS equatorial crossing time. The performance of STIC has been evaluated in comparison to the eddy covariance measurements of $\lambda E$ and $H$ at 30 sites that cover a broad range of biomes and climates. We found a RMSE of 37.79 (11%) (with MODIS Terra $T_S$) and 44.27 $\text{W}\ \text{m}^{-2}$ (15%) (with MODIS Aqua $T_S$) in $\lambda E$ estimates, while the RMSE was 37.74 (9%) (with Terra) and 44.72 $\text{W}\ \text{m}^{-2}$ (8%) (with Aqua) in $H$. STIC could efficiently capture the $\lambda E$ dynamics during the dry down period in the semi-arid landscapes where $\lambda E$ is strongly governed by the subsurface soil moisture and where the majority of other $\lambda E$ models generally show poor results. Sensitivity analysis revealed a high sensitivity of both the fluxes to the uncertainties in $T_S$. A realistic response and modest relationship was also found when partitioned $\lambda E$ components ($\lambda E_E$ and $\lambda E_T$) were compared to the observed soil moisture and rainfall. This is the first study to report the physical integration of $T_S$ into the PM equation and finding analytical solution of the physical ($g_B$) and physiological conductances ($g_S$). The performance of STIC over diverse biomes and climates points to its potential to benefit future NASA and NOAA missions having thermal sensors, such as HyspIRI, GeoSTAR and GOES-R for mapping multi-scale $\lambda E$ and drought.  
© 2013 Elsevier Inc. All rights reserved.

---

## 1. Introduction
Quantitative understanding of the behavior of the latent and sensible heat fluxes is central to a broad range of applications including water resource management (Anderson, Norman, Mecikalski, Otkin & Kustas, 2007; French, Schmugge, Kustas, Brubaker & Prueger, 2003; Norman, Kustas & Humes, 1995), weather forecasting and climate change prediction (Entekhabi, Asrar, Betts, Beven, Bras, Duffy, et al., 1999).

The Penman-Monteith (PM) equation (Monteith, 1965) has become the pre-eminent method for specifying the surface to air latent heat flux, $\lambda E$, (or evapotranspiration, $E$), in a broad range of applications (Choudhury, 1997; Droogers & Allen, 2002; Eslamian, Gohari, Zareian & Firoozfar, 2012; Moran, Rahman, Washburne, Goodrich, Weltz & Kustas, 1996). The central concept in the derivation of this combination equation is the elimination of the need to specify surface temperature ($T_S$) when estimating $\lambda E$. This was originally motivated by the fact that observations of $T_S$ were invariably lacking for the scales at which $\lambda E$ is needed to be specified (Monteith, 1965). However, in solving this problem the PM approach creates further problems, particularly in relation to the specification of the boundary layer (or aerodynamic) and surface (or stomatal) conductance terms ($g_B$ and $g_S$) which are generally not measureable at scales in which the PM equation is applied. The principal solution to this has been the adoption of somewhat speculative models for these conductances, degrading the predictive quality of the PM framework. The speculative nature of the models for $g_B$ stems from several sources (Cleugh, Leuning, Mu & Running, 2007; Mu, Heinsch, Zhao & Running, 2007; Su, 2002) and although credible models for mass, momentum and energy transfer exist (Liu, Lu, Mao & Jia, 2007) these require parameterization in order to make them applicable to the intended scale of use (Raupach & Finnigan, 1995). Furthermore, these parameterizations are not stationary due to the dynamics of the near surface and its boundary layer. For $g_S$ the situation is more problematic in that, in addition to the scale dependence and spatiotemporal heterogeneity, there remains no universally agreed quantitative model that describes the biological controls on water and carbon fluxes (Landsberg, Kaufman, Binkeley, Isebrands & Jarvis, 1991) despite mechanistic (Dewar, 1995; Katul, Manzoni, Palmroth & Oren, 2010; Leuning, 1995) and empirical (Ball, Woodrow & Berry, 1987; Jarvis, 1976) models being widely used.

The nature of the elimination of $T_S$ from the PM equation arises from the fact that the physical feedbacks regulating $\lambda E$ are heavily temperature dependent, making $T_S$ a primary state variable of surface energy balance closures (Kustas & Anderson, 2009). As a result, in situations where observations of $T_S$ are available, such as from remote sensing platforms, this provides a rich source of information that can be exploited to characterize the components of surface energy balance. Recognizing this, and motivated by the advent of thermal remote sensing, an alternative modeling strategy for $\lambda E$ has been to use $T_S$ to solve for the sensible heat flux ($H$) and then estimate $\lambda E$ as a residual of the surface energy balance (Anderson et al., 1997, 2007; Bastiaansssen, Menenti, Feddes & Holtslag, 1998; Kustas & Norman, 1999; Norman et al., 1995; Su, 2002). So, although Hall et al. (1992) expressed serious doubts over using remotely sensed $T_S$ measurements in aerodynamic transfer equations because of the large differences that can exist between the radiometric, $T_{SR}$, and aerodynamic, $T_{SA}$, surface temperatures (Troufleau et7 al., 1997), using a two-source soil-canopy framework (Anderson et al., 1997, 2007; Norman et al., 1995) and inclusion of an 'extra conductance' concept in the single-source framework (Boegh, Soegaard & Thomsen, 2002; Su, 2002) appeared to accommodate the effects of differences between $T_{SR}$ and $T_{SA}$.

However, although useful, the above approaches still rely on the specification of $g_B$ as an exogenous input when clearly it is an internal variable whose state is closely related to that of $T_S$. Therefore, an alternative approach may be to revisit the PM framework to see if $T_S$ can be used to eliminate the need to specify $g_B$ and $g_S$ altogether. This paper focuses on attempting to solve the PM closure problem in such a way that treats $g_B$ and $g_S$ as internal unobserved components and instead exploits remotely sensed observations of $T_S$. For this we combine the PM framework with the advection-aridity hypothesis of Brutsaert & Stricker (1979), which is a modification of Bouchet's complementary hypothesis (Bouchet, 1963). This approach is based on the PMBL (Penman-Monteith-Bouchet-Lhomme) method proposed by Mallick, Jarvis, Fisher, Tu, Boegh & Niyogi (2013) where relative humidity ($R_H$) and atmospheric vapor pressure deficit ($D_A$) information were used to describe the 'wetness' of the land surface. The $T_S$ variant of PMBL considered here uses $T_S$ as an additional data source to obtain the system closure through retrieving the 'near surface' moisture availability ($M$) and 'effective' vapor pressure at the evaporating front ($e_S$) (vapor pressure at the level of $T_S$ measurements). We refer to the scheme as the Surface Temperature Initiated Closure (STIC).

The study objectives of the paper are:
1. To derive a surface energy balance closure expression for $\lambda E$ building on the PM equation framework but re-incorporating $T_S$ in order to eliminate the exogenous expressions for $g_B$ and $g_S$.
2. To evaluate the estimates of $\lambda E$ and $H$ produced by this closure using tower flux data from a broad range of biome and climate types, including a sensitivity analysis to input errors.
3. To evaluate some of the components generated by STIC, in particular evaporation, transpiration, and land surface wetness in comparison with observed hydrological variables (e.g., soil moisture and rainfall).

Section 2 sets out the derivation of STIC. Section 3 describes the data sources used in the study. Section 4 describes a sensitivity analysis of STIC, the evaluation of STIC outputs against the tower measurements and a putative evaluation of the ancillary STIC outputs. Study results are discussed in Section 5, while the strengths and limitations of the proposed approach are described in Section 6. A list of symbols used in the present study is given in Table 1.

---

### Table 1
**Symbols and their description used in the study.**

| Symbol | Description |
| :--- | :--- |
| $\lambda$ | Latent heat of vaporization of water ($\text{J}\ \text{kg}^{-1}\ \text{K}^{-1}$) |
| $\lambda E$ | Evaporation (evaporation + transpiration) as latent heat flux ($\text{W}\ \text{m}^{-2}$) |
| $H$ | Sensible heat flux ($\text{W}\ \text{m}^{-2}$) |
| $R_N$ | Net radiation ($\text{W}\ \text{m}^{-2}$) |
| $R_G$ | Shortwave radiation ($\text{W}\ \text{m}^{-2}$) |
| $G$ | Ground heat flux ($\text{W}\ \text{m}^{-2}$) |
| $\Phi$ | Net available energy ($\text{W}\ \text{m}^{-2}$) |
| $E$ | Evapotranspiration (evaporation + transpiration) as depth of water (mm) |
| $\lambda E_E$ | Evaporation as flux ($\text{W}\ \text{m}^{-2}$) |
| $\lambda E_T$ | Transpiration as flux ($\text{W}\ \text{m}^{-2}$) |
| $\lambda E^*$ | Potential evaporation as flux ($\text{W}\ \text{m}^{-2}$) |
| $\lambda E_T^*$ | Potential transpiration as flux ($\text{W}\ \text{m}^{-2}$) |
| $\lambda E_W$ | Wet environment evaporation as flux ($\text{W}\ \text{m}^{-2}$) |
| $\lambda E_P$ | Potential evaporation as flux ($\text{W}\ \text{m}^{-2}$) according to Penman |
| $\lambda E_{PT}$ | Potential evaporation as flux ($\text{W}\ \text{m}^{-2}$) according to Priestley-Taylor |
| $E_E$ | Evaporation as depth of water (mm) |
| $E_T$ | Transpiration as depth of water (mm) |
| $E^*$ | Potential evaporation as depth of water (mm) |
| $E_T^*$ | Potential transpiration as depth of water (mm) |
| $E_P$ | Potential evaporation as depth of water (mm) according to Penman |
| $E_{PT}$ | Potential evaporation as depth of water (mm) according to Priestley-Taylor |
| $E_W$ | Wet environment evaporation as depth of water (mm) |
| $g_B$ | Boundary layer conductance ($\text{m}\ \text{s}^{-1}$) |
| $g_S$ | Stomatal/surface conductance ($\text{m}\ \text{s}^{-1}$) |
| $M$ | Surface moisture availability (0-1) |
| $T_A$ | Air temperature (°C) |
| $T_D$ | Dewpoint temperature (°C) |
| $T_S$ | Radiometric surface temperature (°C) |
| $T_{SD}$ | Surface dew-point temperature (°C) |
| $T_{SA}$ | Aerodynamic surface temperature (°C) |
| $T_{SR}$ | Radiometric surface temperature (°C) |
| $R_H$ | Relative humidity (%) |
| $e_A$ | Atmospheric vapor pressure (hPa) at the level of $T_A$ measurement |
| $D_A$ | Atmospheric vapor pressure deficit (hPa) at the level of $T_A$ measurement |
| $e_S$ | 'effective' vapor pressure of evaporating front near the surface (hPa) |
| $e_S^*$ | Saturation vapor pressure of surface (hPa) at $T_S$ |
| $s$ | Slope of saturation vapor pressure versus temperature curve ($\text{hPa}\ \text{K}^{-1}$) |
| $\gamma$ | Psychrometric constant ($\text{Pa}\ \text{K}^{-1}$) |
| $\rho$ | Density of air ($\text{kg}\ \text{m}^{-3}$) |
| $c_P$ | Specific heat of dry air ($\text{MJ}\ \text{kg}^{-1}\ \text{K}^{-1}$) |
| $\Lambda$ | Evaporative fraction |
| $\beta$ | Bowen ratio |
| $\theta$ | Surface (0-5 cm) soil moisture ($\text{m}^3\ \text{m}^{-3}$) |
| $\theta_d$ | Deep layer (15-30 cm) soil moisture ($\text{m}^3\ \text{m}^{-3}$) |
| $L$ | Leaf area index |

---

## 2. Derivation of STIC

![Fig. 1. Schematic representation of one-dimensional description of STIC.](https://placeholder-for-figure1)  
*Fig. 1. Schematic representation of one-dimensional description of STIC. Here $g_B$ and $g_S$ are the aerodynamic and stomatal conductances, $e_S^*$ is the saturation vapor pressure of the surface, $T_{SA}$ is the aerodynamic surface temperature that is responsible for transferring the sensible heat ($H$), $e_S$ is the near surface vapor pressure of the evaporating front, $T_S$ is the radiometric surface temperature, $T_{SD}$ is the surface dewpoint temperature, $M$ is the surface moisture availability, $T_A$ and $e_A$ are temperature and vapor pressure of reference height, and $\lambda E$ is the latent heat flux, respectively.*

The conceptual framework for STIC is given in Fig. 1. The aim here is to define a surface energy balance equation set that treats $g_B$ and $g_S$ as internal states and $T_S$ as external to the solution. This set will consist of four equations and four unknowns which will be solved analytically. The four unknowns are: $g_B$, $g_S$, the vertical difference between the source/sink height temperature ($T_{SA}$) (also called the aerodynamic temperature) and air temperature ($T_A$) $\Delta T$, and the evaporative fraction, $\Lambda$.

The PM equation states (Monteith, 1965),
$$\lambda E = \frac{s\Phi + \rho c_P g_B D_A}{s + \gamma \left(1 + \frac{g_B}{g_S}\right)} \qquad (1)$$

where $\rho$ is the density of dry air ($\text{kg}\ \text{m}^{-3}$), $c_P$ is the specific heat of dry air ($\text{MJ}\ \text{kg}^{-1}\ \text{K}^{-1}$), $\gamma$ is the psychrometric constant ($\text{hPa}\ \text{K}^{-1}$), $s$ is the slope of the saturation vapor pressure versus air temperature ($\text{hPa}\ \text{K}^{-1}$), $D_A$ is the saturation deficit of the air (hPa) at the reference level and $\Phi$ is the net available energy ($\text{W}\ \text{m}^{-2}$) (the difference between net radiation, $R_N$, and ground heat flux, $G$). The two unknowns in Eq. (1) are $g_B$ and $g_S$ and derivations of the physical expressions for these two unknowns are given in Appendix A (see also Mallick et al., 2013). However, these expressions require specification of the near surface moisture availability, $M$, and evaporative fraction, $\Lambda$.

$M$ is a unitless quantity which signifies the relative dryness or wetness of the surface and hence controls the potential evaporation rate. Considering the general case of evaporation from any nonsaturated surface at a rate less than the potential, $M$ is the ratio of the evaporation rate, $E$, to the potential evaporation rate, $E^*$. $M$ is assumed to be homogeneous between the surface and the evaporation front and its contribution to the effective vapor pressure of the evaporating front near the surface ($e_S$) is given as follows (Lee & Pielke, 1992; Segal, Jia, Ye & Pielke, 1990):
$$e_S = e_A(1 - M) + M e_S^* \qquad (2)$$

where $e_S^*$ is the saturation vapor pressure at the surface (hPa), normally expressed as a function of $T_S$, and $e_A$ is the atmospheric vapor pressure (hPa) at the level at which $T_A$ is measured. $M = 0$ indicates complete dryness whereas $M = 1$ reflects a saturated evaporating front ($e_S \rightarrow e_S^*$) (i.e., after a heavy rainfall event or irrigation). When the evaporating front is absolutely dry $e_S = e_A$. Since $e_A > 0$, the initial moisture gradient towards the evaporating front eventually results in $e_S = e_A$ (no dew is formed during the day). When $e_S = e_A$, no latent heat flux is transported from the surface to the atmosphere, as anticipated in the case of a dry surface.

Fisher, Tu & Baldocchi (2008) proposed an empirical expression for $M$ in terms of $D_A$ and $R_H$ ($M = R_H^{D_A}$), assuming equilibrium between the surface moisture and atmospheric moisture. An atmosphere with low $D_A$ and high $R_H$ indicates a moist humid atmosphere, but the underlying surface may be intermediately dry to very dry. According to the Fisher et al. (2008) expression for $M$, such cases will portray a wet surface condition when actually it may be dry leading to $\lambda E$ being overestimated. Similarly an atmosphere with high $D_A$ and low $R_H$ will portray a dry surface condition when actually it is wet, leading to $\lambda E$ being underestimated. The surface and atmospheric moisture may only correspond when soil moisture is tightly coupled to the atmospheric moisture or radiation (Small & Kurc, 2003). As a result, although useful, the Fisher et al. (2008) expression generally results in the overestimation of $\lambda E$ for $D_A$ ranging from 10 to 20 hPa and above (Mallick et al., 2013). Because $T_S$ is far more sensitive to the land surface moisture conditions (Kustas & Anderson, 2009), here we use $T_S$ in conjunction with $T_A$ and $R_H$ to attempt to retrieve $M$ within a physical framework. $T_S$ serves as a direct metric for the soil moisture status and vegetation conditions, which in turn influences the surface energy fluxes and their partitioning (Kustas & Anderson, 2009). The retrieval of $M$ and the expression for $\Lambda$ is a key novelty of the STIC framework and hence these are described in detail in the main body of the paper. The equations for rest of the variables and their derivations are described in Appendix A.

### 2.1. Derivation of M
Like Fisher et al. (2008) we hypothesize that the moisture availability at the surface and at the evaporating front are similar and, therefore, $M$ is derived from the surface-atmosphere information. According to Noilhan & Planton (1989), Ye & Pielke (1993) and Boegh et al. (2002), the transfer of $\lambda E$ from the surface can be written as follows:
$$\lambda E = \frac{\rho c_P}{\gamma}g_B(e_S - e_A) = M E_P = \frac{\rho c_P}{\gamma} M g_B(e_S^* - e_A) \qquad (3)$$

From Eq. (3), a physical expression for $M$ is given in terms of the vapor pressure gradients:
$$M = \frac{(e_S - e_A)}{(e_S^* - e_A)} \qquad (4)$$

The main difficulty in applying Eq. (4) is that $e_S$ is unknown and there is no straightforward way to relate $e_S$ to $T_S$. Although the Clausius-Clapeyron relationship between vapor pressure and surface temperature is not linear, it is generally linearized for small temperature differences (Monteith, 1965). Fig. 2(a and b) shows the relationships between $e_S$, $e_S^*$ and $e_A$ and their corresponding temperatures.

![Fig. 2. Conceptual diagrams of the saturation vapor pressure curve.](https://placeholder-for-figure2)  
*Fig. 2. (a) Conceptual diagram of the saturation vapor pressure curve and the relationship among $T_{SD}$, $T_S$, $e_S^*$ and $e_S$ in the context of surface at temperature $T_S$ according to Venturini et al. (2008). The air and dewpoint temperature of the overlying air is characterized by $T_A$ and $T_D$ with the vapor pressure $e_A$ and $e_A^*$. Here $T_{SD}$ is the surface dew point temperature, $T_S$ is the radiometric surface temperature, $e_S$ is the near surface vapor pressure and $e_S^*$ is the surface saturation vapor pressure. (b) Conceptual diagram of the linearized saturation vapor pressure curve to demonstrate the relationship between $(e_S - e_A)$ with $s_1(T_{SD} - T_D)$, $(e_S^* - e_S)$ with $s_3(T_S - T_{SD})$ and $(e_S^* - e_A)$ with $s_2(T_S - T_D)$. Here $s_1$, $s_2$ and $s_3$ are the slope of the saturation vapor pressure and temperature curve linearized according to Monteith (1965).*

By analogy to the dew point temperature, $T_D$, if the surface air is brought to saturation without affecting $e_S$ then $e_S = f(T_{SD})$. Thus, $T_{SD} < T_S$ for unsaturated surface and $T_{SD} \rightarrow T_S$ as the surface tends to saturation. Therefore, $e_S$ is the surface saturation vapor pressure at some notional surface dew point $T_{SD}$. Although $T_{SD}$ cannot be directly measured, it can be derived from the slope of the relevant temperature-saturation vapor pressure curve. Therefore, from Eq. (4), Fig. 2, and following Monteith (1965), $M$ can be written as follows,
$$M = \frac{(e_S - e_A)}{(e_S^* - e_A)} = \frac{s_1(T_{SD} - T_D)}{s_2(T_S - T_D)} \qquad (5)$$

where $s_1$ and $s_2$ are the slopes of the saturation vapor pressure and temperature between $(T_{SD} - T_D)$ versus $(e_S - e_A)$ and $(T_S - T_D)$ versus $(e_S^* - e_A)$ relationship. For a dry surface $T_S \gg T_{SD}$, $e_S^* \gg e_S$, $(T_S - T_D) \gg (T_{SD} - T_D)$ and $M \rightarrow 0$. For any wet surface, $e_S^* \cong e_S$, $T_S \cong T_{SD}$, $(T_S - T_D) \cong (T_{SD} - T_D)$ and $M \rightarrow 1$. According to Eq. (5), when condensation occurs, the surface temperature goes below the ambient air temperature and $e_S < e_A$, otherwise $e_S \ge e_A$ and $e_S \le e_S^*$. Hence, $T_{SD}$ falls somewhere between $T_S$ and $T_D$ (Fig. 2a) and can be computed by linearizing the saturation vapor pressure curve between the two levels $(T_D, T_{SD})$ and $(T_{SD}, T_S)$ (Monteith, 1965).

Earlier, Venturini, Islam & Rodriguez (2008) assumed $s_1 = s_2$ and hence they cancel each other. Here, both $s_1$ and $s_3$ are approximated from the corresponding Clausius-Clapeyron relationships. Since $T_S$ and $e_A$ are available, $s_2$ can be calculated directly. According to Fig. 2,
$$s_1 = \frac{(e_S - e_A)}{(T_{SD} - T_D)} \qquad (6)$$
$$s_3 = \frac{(e_S^* - e_S)}{(T_S - T_{SD})} \qquad (7)$$

Combining Eqs. (6) and (7) we can find an expression of $T_{SD}$:
$$T_{SD} = \frac{(e_S^* - e_A) - s_3 T_S + s_1 T_D}{(s_1 - s_3)} \qquad (8)$$

Among the various exponential functions that relate the saturation vapor pressure and temperature, Buck's (Buck, 1981) equation was used for its simplicity. After estimating the saturation vapor pressure at different temperatures, the initial values of $s_1$ and $s_3$ are specified at $T_D$ and $T_S$. Then using Eq. (8) an initial value of $T_{SD}$ is obtained. The process is then iterated by updating $s_3$ with the first estimate of $T_{SD}$ and computing a new $T_{SD}$. A large proportion of the saturation vapor pressure curve is assumed to be linearized in the first calculation and therefore $s_1$ is not updated. Repeating this process produces stable values of $T_{SD}$ in 3-5 iterations. The final $T_{SD}$ values are used in Eq. (5) for obtaining the estimates of $M$.

### 2.2. Expression of $\Lambda$
To close the system of equations we need an expression for the evaporative fraction, $\Lambda$, and this expression must include the dependence of $\Lambda$ on the conductances. For this we exploited two different representations of evaporation: the Penman (P) equation (Penman, 1948), and the Priestley-Taylor (PT) equation (Priestley & Taylor, 1972). Here, these two expressions are related to each other through the complementary relationship advection-aridity hypothesis developed by Brutsaert & Stricker (1979).

According to the original complementary hypothesis (Bouchet, 1963), for a large homogeneous area of 1 to 10 km and away from sharp environmental discontinuities there exists a complementary feedback mechanism between potential evaporation ($\lambda E^*$), evaporation ($\lambda E$) and sensible heat flux ($H$). $\lambda E^*$ is defined as the evaporation from a wet surface under the prevailing atmospheric condition, limited only by the amount of available energy. If moisture at the surface is unlimited (i.e., when $M = 1$), $\lambda E = \lambda E^*$ and this condition is referred to as the wet-environment evaporation ($\lambda E_W$). Bouchet (1963) hypothesized that under the conditions of constant energy supply to a given surface-atmosphere system, when the surface moisture becomes limited ($M \rightarrow 0$), $\lambda E$ falls below $\lambda E_W$ and some amount of energy becomes 'available'. This extra energy increases the temperature and humidity gradient of the overlying air (in the form of sensible heat or longwave back radiation) and leads to an increase in $\lambda E^*$ whose magnitude is equal to the decrease in $\lambda E$. If moisture availability is increased, $\lambda E$ again starts increasing and $\lambda E^*$ decreases accordingly. As a result $\lambda E^*$ is predicated upon the prevailing moisture availability conditions. If the energy budget remains unchanged and all the excess energy is converted into sensible heat, a complementary relationship exists of the following form:
$$\lambda E + \lambda E^* = 2\lambda E_W \qquad (9)$$

This feedback mechanism helps bypass the need for detailed knowledge of the complex processes and interaction between soil, vegetation, and near surface boundary layer. Based on Bouchet's earlier work, Brutsaert & Stricker (1979) proposed an advection-aridity hypothesis that allows the formulation of $\lambda E$ under non-potential conditions. Their approach is based on a conceptual framework where the aridity (i.e., surface dryness) was deduced from the drying power of large scale advection as implied by the atmospheric conditions in the turbulent surface layer. According to Brutsaert & Stricker (1979), $\lambda E_W$ was approximated as the potential evaporation according to Priestley & Taylor (1972), $\lambda E_{PT}^*$, which represents the potential evaporation under the conditions of minimal advection and $\lambda E^*$ was approximated as the potential evaporation according to Penman (1948), $\lambda E_P^*$, in order to capture the effects of large scale advection. Therefore, any large scale advection effects causing $\lambda E_P^* > \lambda E_{PT}^*$ were assumed to be due to the regional surface moisture conditions and also due to forced convection that increases the atmospheric vapor pressure deficit (Brutsaert & Stricker, 1979). Thus actual evapotranspiration could be computed by means of Eq. (9) assuming $\lambda E^* = \lambda E_P^*$ and $\lambda E_W = \lambda E_{PT}^*$:
$$\lambda E + \lambda E_P^* = 2\lambda E_{PT}^* \qquad (10)$$

This approach was found to yield reasonable estimates of $\lambda E$ across both well watered and drought conditions (Brutsaert & Stricker, 1979; Huntington, Szilagyi, Tyler & Pohll, 2011; Parlange & Katul, 1992). The main advantage of this approach is that it does not require a sub-model for representing the stomatal conductance, soil moisture or any other land surface measures of aridity. Taking advantage of this advection-aridity hypothesis we are able to express the evaporative fraction in terms of the two conductances ($g_B$ and $g_S$) and therefore we are able to close the system of equations in the present scheme (see Section 2.3 and Appendix A).

One obvious question arises over the interdependence of the P and PT equations because the PT equation was derived from Penman's equation. According to Penman (1948), $\lambda E_P^*$ is written as follows:
$$\lambda E_P^* = \frac{s\Phi + \rho c_P g_B D_A}{s + \gamma} \qquad (11)$$

The first term of Eq. (11) may be considered to represent the lower limit of evaporation from a moist surface, which was referred to as 'equilibrium evaporation' by Slatyer & McIlroy (1961) or diabatic evaporation by Penman (1948). The second term of Eq. (11) is a measure of the departure from this equilibrium state in the atmosphere and is often referred to as the adiabatic evaporation (Penman, 1948) or 'drying power' of air (Granger & Gray, 1989). In the absence of radiative divergence, this adiabatic departure would stem from large-scale advection involving horizontal variations of surface variables (e.g., $T_S$, soil moisture etc.) and atmospheric conditions. The atmospheric boundary layer is not uniform and tends to maintain a humidity deficit, even over the ocean. Therefore, true equilibrium conditions are rarely encountered over a wet surface. Analyzing data collected over various surfaces Priestley & Taylor (1972) proposed an empirical adjustment of the equilibrium evaporation expression,
$$\lambda E_{PT}^* = \frac{\alpha s \Phi}{s + \gamma} \qquad (12)$$

where $\alpha$ is a compensation factor introduced due to neglecting the adiabatic term in the Penman equation. They decided that for a large saturated landscape with minimum advection the best estimate of $\alpha$ is 1.26 (Priestley & Taylor, 1972). Therefore, the interdependence of $\lambda E_{PT}^*$ and $\lambda E_P^*$ is not by chance. When the underlying surface is wet, the air is close to saturation and evaporation is dominated by the radiation term of the Penman equation. But, $\lambda E_{PT}^*$ never equals $\lambda E_P^*$ because of the involvement of the 'drying power' term in $\lambda E_P^*$.

It should also be emphasized that although the Priestley-Taylor $\alpha$ is an input into our scheme, the sensitivity analysis over a large range of $\alpha$ values revealed that STIC is only mildly sensitive to $\alpha$ (see Section 4.1 and Fig. 5b).

The main restriction of Bouchet's hypothesis can be relaxed through the introduction of this advection-aridity hypothesis (Brutsaert & Stricker, 1979) that accounts for advection through the adiabatic term of $\lambda E_P^*$. Under some circumstances $\lambda E_P > \Phi$, which is not only due to the availability of excess energy that increases $\lambda E_P^*$, but also because of horizontally advected dry air which increases the 'drying power' and hence $\lambda E_P^*$. The horizontal advection is pronounced when a dry front passes over a surface and increases both $D_A$ and $g_B$ and consequently the drying power of air, all of which are accounted for in the estimation of $\lambda E_P^*$ (Brutsaert & Stricker, 1979; Huntington et al., 2011; Parlange & Katul, 1992). Moving away from any idealized condition, the stronger the advection and the greater the surface dryness, the more $\lambda E$ differs from its equilibrium rate. Taken as a whole, these findings would suggest that the advection-aridity hypothesis is relatively robust.

According to Budyko (Budyko, Efimova, Zubenok & Strokina, 1962; see also Roderick & Farquhar, 2004), when complementarity holds, the regional evaporation rate is limited by moisture availability in arid climates and energy availability in humid climates. However, the complementary relationship allows regional $\lambda E^*$ to depend on regional $\lambda E$ in a complementary manner throughout any range of moisture and energy availability (Ramirez, Hobbins & Brown, 2005). Some theoretical arguments suggest that the hypothesis of 1:1 compensation between $\lambda E$ and $\lambda E^*$ is only partially fulfilled (Lhomme, 1997; Sugita, Usui, Tamagawa & Kaihotsu, 2001). However, more recently Ramirez et al. (2005) found observational evidence for 1:1 compensation between $\lambda E$ and $\lambda E^*$ and, therefore, we also assume perfect compensation here.

This preamble is important because Eqs. (1), (10), (11) and (12) can be combined and algebraically reorganized to give the following expression for $\Lambda$ (please see Appendix A for details):
$$\Lambda = \frac{2\alpha s}{2s + \gamma \left(2 + \frac{g_B}{g_S}\right)} \qquad (13)$$

### 2.3. The STIC closure equations
The closure equations for STIC are as follows and their derivations are detailed in Appendix A:
$$g_B = \frac{\Phi}{\rho c_P \left(\Delta T + \frac{e_S - e_A}{\gamma}\right)} \qquad (14)$$
$$g_S = g_B \frac{(e_S - e_A)}{(e_S^* - e_S)} \qquad (15)$$
$$\Delta T = \left(\frac{e_S - e_A}{\gamma}\right)\left(\frac{1 - \Lambda}{\Lambda}\right) \qquad (16)$$
$$\Lambda = \frac{2\alpha s}{2s + \gamma \left(2 + \frac{g_B}{g_S}\right)} \qquad (17)$$

These four equations provide constraints for the four unknowns $g_B$, $g_S$, $\Delta T$, and $\Lambda$. The computational diagram is given in Fig. 3. Eqs. (14)-(17) were solved to retrieve the analytical expressions of the four unobserved components.

![Fig. 3. Diagrammatic representation of the core equations.](https://placeholder-for-figure3)  
*Fig. 3. Diagrammatic representation of the core equations used to recover the internal state variables in STIC. All the variables are explained in the main body of the manuscript.*

### 2.4. Partitioning $\lambda E$ estimates
The terrestrial latent heat flux is an aggregate of both evaporation ($\lambda E_E$) and transpiration ($\lambda E_T$). When it rains the land surface becomes wetter and $\lambda E$ tends to the potential evaporation ($\lambda E^*$), while surface drying after rainfall leads to $\lambda E$ tending to the potential transpiration rate ($\lambda E_T^*$) in the presence of vegetation, or zero without any vegetation. $\lambda E$ at any time is, therefore, some blend of these two end member conditions depending on the degree of surface moisture availability (Bosveld & Bouten, 2003; Loescher, Gholz, Jacobs & Oberbauer, 2005). Given moisture availability ($M$) is a STIC output, the separation of lumped $\lambda E$ into $\lambda E_E$ and $\lambda E_T$ was also tested:
$$\lambda E = \lambda E_E + \lambda E_T = M \lambda E^* + (1 - M)\lambda E_T \qquad (18)$$

Having recovered $g_B$, $\lambda E^*$ can be estimated according to the Penman equation (Eq. 11) and $\lambda E_T$ can be estimated as the residual in Eq. (18). $\lambda E_T^*$ can be estimated as $\lambda E_T/(1 - M)$.

---

## 3. Datasets
Estimation of $\lambda E$ and $H$ through STIC requires information on $R_N$, $G$, $T_A$, $R_H$ or $e_A$ and $T_S$. The first four variables were available from FLUXNET sites and $T_S$ was obtained from MODIS. Given the absence of dewpoint temperature ($T_D$), it was calculated from $T_A$ and $R_H$ according to Buck (1981) as follows:
$$T_D = \frac{c \cdot \gamma_M(T_A, R_H)}{b - \gamma_M(T_A, R_H)} \qquad (19)$$
$$\gamma_M(T_A, R_H) = \log\left(R_H \cdot \exp\left(\left(b - \frac{T_A}{d}\right)\left(\frac{T_A}{c + T_A}\right)\right)\right) \qquad (20)$$

Where $b$, $c$, and $d$, are empirical coefficients having values of 18.678, 257.14 and 234.5, respectively. Detailed descriptions of the different datasets are given below.

### 3.1. Eddy Covariance (EC) tower data (FLUXNET data)
These data include 30 EC sites covering different sub-networks of FLUXNET (Ameriflux, Euroflux, Ozflux, FLUXNET-Canada and CarboAfrica) (Baldocchi et al., 2001) (Table 2). The data covers a broad spectrum of biomes and climate types such as the semi-arid tropics, subtropics, Mediterranean savannas, temperate grasslands, temperate forests, boreal forests and arctic wetlands. Distribution of sites (Fig. 4 and Table 2) shows that more than 50% of the total sites represent semi-arid landscapes. The main reason for selecting majority of the semi-arid sites is because these regions are characterized by strong land-atmosphere coupling (Diremeyer, Cash, Kinter, Stan, Jung, Marx, et al., 2012; Kurc & Small, 2004), which is often confounded by extreme heterogeneity in soil moisture and vegetation (Strasser & Mauser, 2001). Such features make the specifications of the surface energy balance components more sensitive to the boundary conditions. In addition, limited water resources and growing needs of urban as well as agricultural water requirements in semi-arid sites are continuously increasing the pressure on accessible groundwater and could introduce dynamical changes in the hydrological cycle (Niyogi, Kishtawal, Tripathi & Govindaraju, 2010; Pielke, Pitman, Niyogi, Mahmood, McAlpine, Hossain, et al., 2011).

EC data were obtained from the FLUXNET LaThuile database (www.fluxdata.org) and CarboAfrica network website (www.gaia.agraria.unitus.it/database/carboafrica). Independent measurements of $\lambda E$ and $H$ were available along with the measurements of $R_N$, $G$, $T_A$ and $R_H$. Soil moisture and precipitation data were also available through time domain reflectometer and rain gauge measurements. The energy balance closure fraction was estimated as $C_f = (\lambda E + H)/(R_N - G)$ by linear regression between the turbulent fluxes against the net available energy, forced through the origin (Barr, Morgenstern, Black, McCaughey, & Nesic, 2006). The closure fraction was found to vary from 75 to 80%. The inherent errors associated with measurements of $H$, $\lambda E$, $R_N$ and $G$ using eddy covariance are reported to be 15-20% (Weaver, 1990), 15-20% (Field, Fritschen, Kanemasu, Smith, Stewart, Verma, et al., 1992), 5-10% (Twine et al., 2000) and 20-30% (Twine et al., 2000), respectively. This leaves a considerable proportion of energy unaccounted for in the partitioning of $R_N$ to $H$ and $\lambda E$. This could lead to significant discrepancies when estimated fluxes are compared with the observations. Therefore, the measured Bowen ratio (Bowen, 1926) of the EC towers was used to adjust the measured fluxes and surface energy balance was closed as suggested by Twine et al. (2000) and later adopted by Chavez, Neale, Hipps, Prueger & Kustas (2005) and Anderson et al. (2007). For all the sites, the gap-filled half-hourly level 3 data were used. These include meteorological ($T_A$ and $R_H$), micrometeorological ($\lambda E$, $H$, $R_N$, and $G$), and associated hydrological variables (soil moisture and rainfall). The meteorological and micrometeorological data co-incident to MODIS Terra and Aqua overpass times were averaged for every eight-day (corresponding to MODIS $T_S$ data which are eight-day products) and the eight-day data were used for running STIC. Similarly, eight-day averages of the $\lambda E$ and $H$ measurements corresponding to MODIS overpass were used to validate the estimated fluxes.

### 3.2. MODIS data
MODIS terra has the local equatorial crossing time at approximately 10:30 AM in a descending node while Aqua has the local equatorial crossing time at approximately 1:30 PM in an ascending node with a sun-synchronous, near-polar, circular orbit. Tiles of the MODIS Terra and Aqua eight-day level-3 $T_S$ data (MOD11A2 and MYD11A2 version 005) at 1 km spatial resolution was acquired through the Oak Ridge National Laboratory (ORNL) data gateway (http://daac.ornl.gov/MODIS/). This is the harmonized $T_S$ data maintained by ORNL over all the flux measurement sub-networks. The data is a $7 \times 7$ pixel cutout centred on the EC flux measurement sites. From these cutouts, $T_S$ of $2 \times 2$ pixels from the central location were again extracted and averaged. Only quality controlled data (with a QA flag of 0) were used in the pixel averaging. Any bad quality data within $2 \times 2$ pixels were ignored while averaging. Despite MODIS $T_S$ having 1 km spatial resolution and the tower meteorological variables (e.g., $T_A$ and $R_H$) having a fetch <1 km, still MODIS $T_S$ is the best available time-series data that covers multiple seasonal cycles with minimum data gaps. Similarly, FLUXNET is the best available database to estimate and validate the surface fluxes covering multiple biomes and climate regimes.

---

### Table 2
**Description of the sites including latitude, longitude, biome types, climate, annual precipitation (mm) and years of data used in the study.**

| Site name (Code) | Country | Latitude | Longitude | Biome | Climate type | EC height (m) | Annual rainfall (mm) | Years | Reference |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Audubon Ranch (AUD) | US | 31.59 | -110.51 | GRA | Dry arid and Semi-arid | 5 | 440 | 2002-2006 | Krishnan et al. (2012) |
| Freeman Ranch (FR2) | US | 29.95 | -97.99 | WSA | Subtropical Mediterranean | 7 | 865 | 2004-2006 | Heinsch et al. (2004) |
| Santarita Mesquite (SRM) | US | 31.82 | -110.87 | WSA | Dry arid and Semi-arid | 7 | 480 | 2004-2006 | Scott et al. (2009) |
| Sky Oaks (SO2) | US | 33.37 | -116.62 | CSH | Subtropical Mediterranean | n/a | 555 | 2004-2006 | Lipson et al. (2005) |
| Tonzi Ranch (TON) | US | 38.43 | -120.97 | WSA | Subtropical Mediterranean | 23 | 525 | 2002-2006 | Baldocchi and Xu (2007) |
| Walnut Gulch (WKG) | US | 31.74 | -109.94 | GRA | Dry arid and Semi-arid | 7 | 410 | 2004-2006 | Scott et al. (2010) |
| Flagstaff-Managed Forest (FMF) | US | 35.14 | -111.73 | WSA | Dry arid and Semi-arid | 23 | 625 | 2005-2006 | Dore et al. (2010) |
| Flagstaff-Unmanaged Forest (FUF) | US | 35.09 | -111.76 | WSA | Dry arid and Semi-arid | 23 | 650 | 2005-2006 | Dore et al. (2010) |
| Flagstaff-Wildfire (FWF) | US | 35.45 | -111.77 | WSA | Dry arid and Semi-arid | 23 | 620 | 2005-2006 | Dore et al. (2010) |
| Sao Paulo Cerrado (SP1) | Brazil | -21.62 | -47.65 | WSA | Tropical | 50 | 1420 | 2001-2002 | Santos et al. (2004) |
| Howard Springs (HOW) | Australia | -12.49 | 131.15 | WSA | Tropical | n/a | 1450 | 2002-2006 | Beringer et al. (2003) |
| Foggdam (FOG) | Australia | -12.54 | 131.31 | WSA | Tropical | n/a | 1325 | 2006-2007 | - |
| Las Majadas del Tieter (LMA) | Spain | 39.94 | -5.77 | SAV | Subtropical Mediterranean | n/a | 370 | 2004-2006 | - |
| Demokeya (DEM) | Sudan | 13.28 | 30.48 | GRA | Semi-arid tropics | n/a | 320 | 2007-2009 | Ardö et al. (2008) |
| Maun Mopane (MA1) | Botswana | -19.92 | 23.56 | WSA | Dry arid and Semi-arid | 13.5 | 495 | 2000-2001 | Veenendaal et al. (2004) |
| Kruger National Park (KRU) | South Africa | -25.02 | 31.49 | SAV | Subtropical Mediterranean | n/a | 525 | 2001-2003 | Williams et al. (2009) |
| Willow Creek (WCR) | US | 45.81 | -90.08 | DBF | Temperate continental | 30 | 790 | 2001-2006 | Cook et al. (2004) |
| Hainich (HAI) | Germany | 51.08 | 10.45 | DBF | Temperate | 42 | 780 | 2002-2006 | Grünwald and Bernhofer (2007) |
| Univ. of Michigan Station (UMB) | US | 45.56 | -84.71 | DBF | Temperate continental | 46 | 805 | 2001-2003 | Curtis et al. (2002) |
| Hesse forest (HES) | France | 48.67 | 7.06 | DBF | Temperate | 22 | 795 | 2002-2006 | Granier et al. (2000) |
| Cabauw (CAI) | Netherlands | 51.97 | 4.93 | GRA | Temperate | n/a | 775 | 2003-2006 | Sulkava et al. (2011) |
| Wetztein (WET) | Germany | 50.45 | 11.45 | ENF | Temperate | 25 | 870 | 2002-2006 | Rebmann et al. (2010) |
| Mehrstadt (MEH) | Germany | 51.27 | 10.66 | GRA | Temperate | n/a | 695 | 2004-2006 | Sulkava et al. (2011) |
| Grillenburg (GRI) | Germany | 50.95 | 13.51 | GRA | Temperate | 5 | 670 | 2005-2005 | Sulkava et al. (2011) |
| SSA Old Aspen (OAS) | Canada | 53.63 | -106.20 | DBF | Boreal | n/a | 430 | 2001-2005 | Hilker et al. (2011) |
| Lethbridge (LET) | Canada | 49.709 | -112.940 | GRA | Temperate | 5 | 400 | 2001-2005 | Flanagan et al. (2002) |
| Loobos (LOO) | Netherlands | 52.17 | 5.74 | ENF | Temperate | 52 | 785 | 2001-2006 | Sulkava et al. (2011) |
| Saskatchewan Fire 1977 (SF1) | Canada | 53.628 | -106.198 | ENF | Boreal | n/a | 425 | 2004-2005 | Amiro et al. (2006) |
| Ontario Turkey Point (TP1) | Canada | 42.661 | -80.559 | ENF | Boreal | 10 | 910 | 2004-2005 | Peichl and Arain (2006) |
| Ivotuk (IVO) | US | 68.486 | -155.75 | WET | Arctic | n/a | 305 | 2004-2006 | Epstein et al. (2004) |

*GRA = grassland; WSA = woody savannah; SAV = savannah; CSH = closed shrubland; DBF = deciduous broadleaf forest; ENF = evergreen needleleaf forest; WET = arctic wetland.*

---

![Fig. 4. Distribution of eddy covariance sites.](https://placeholder-for-figure4)  
*Fig. 4. Distribution of eddy covariance sites (small yellow dots) used in the present analysis spanned over different biome and climatic settings of FLUXNET covering five continents.*

---

## 4. Results
The estimates of $\lambda E$ and $H$ made through STIC were first subjected to a sensitivity analysis, validation and a residual error analysis of the estimated $\lambda E$. This was followed by an evaluation of $M$ and $\lambda E$ components based on independently measured hydrological variables (soil moisture and rainfall) in the EC tower sites. Detailed evaluation of the individual conductance and $\Delta T$ is beyond the scope of this paper because no direct measurements of these variables were available.

### 4.1. Sensitivity analysis
The accuracy of STIC heavily depends on the quality of MODIS land surface temperature data due to its role in retrieving $T_{SD}$ and $M$. Therefore, a sensitivity analysis of the scheme was first carried out to quantify the impacts of uncertainty in $T_S$ on $\lambda E$ and $H$. The nature of the sensitivity analysis used here is similar to that of Anderson et al. (1997): the absolute sensitivity ($S_V$) of any of the output variables ($V$) to $X$ uncertainty in $T_S$ was assigned as $S_V = |(V_{X+} - V_{X-})/V_{Xt}|$. Here $V_{X+}$ and $V_{X-}$ are the estimated variables when the value of $T_S$ is increased or decreased by $X$ and $V_{Xt}$ is the value of the estimated variable at actual $T_S$. Sensitivity was tested over four different biome types (cropland, grassland, savanna and forest) using actual MODIS $T_S$ uncertainty over the North American FLUXNET sites as reported by Hulley, Hughes & Hook (2012). One representative site from each biome class was selected where the actual MODIS $T_S$ uncertainties were 1.78 K (cropland, e.g., Bondville), 1.34 K (grassland, e.g., Walnut Gulch), 2.66 K (savanna, e.g., Freeman ranch) and 2.24 K (forest, e.g., Duke forest) (Hulley et al., 2012), respectively. Both $\lambda E$ and $H$ were found to be sensitive to the reported $T_S$ uncertainties with a magnitude of $S_V$ to the order of 24-33% (for $\lambda E$) and 17-44% (for $H$) (Fig. 5a).

Since the Priestley-Taylor parameter ($\alpha$) was involved in the entire scheme and there are reports of wide variability of $\alpha$ (Komatsu, 2005), the sensitivity of STIC to $\alpha$ was also tested for large variations in $\alpha$. The reported values of $\alpha$ in the literature are found to vary between 1 and 1.5 (Baldocchi & Xu, 2007). Taking this range we found that for $\lambda E$ the magnitude of $S_V$ varied from 10 to 13% and for $H$ it varied from 9 to 16% (Fig. 5b). This suggests that STIC is not particularly sensitive to $\alpha$ and assigning $\alpha = 1.25$ does not introduce significant error in the $\lambda E$ and $H$ outputs. Maximum RMSE difference in $\lambda E$ and $H$ was 10 and 15 $\text{W}\ \text{m}^{-2}$, respectively.

Sensitivity of STIC was also tested for the atmospheric vapor pressure ($e_A$) and both the fluxes were not very sensitive to $\pm 10\%$ uncertainty in the $e_A$ values (Fig. 5c). The magnitude of $S_V$ varied between 5-6% for $\lambda E$ and 3-9% for $H$. The typical errors of water vapor retrieval from the current MODIS atmospheric data and Atmospheric Infrared Sounder (AIRS) are 5-10% (Gao & Kaufman, 2003) and 10% (Tobin, Revercomb, Knuteson, Lesht, Strow, Hannon, et al., 2006), respectively. With the advent of new sensor technology and the scheduled launch of a future geostationary interferometer sounder GEOSTAR, the water vapor retrieval error will be within 5%, which will facilitate the spatial implementation of STIC for large area mapping of surface fluxes.

![Fig. 5. Sensitivity of the STIC derived fluxes.](https://placeholder-for-figure5)  
*Fig. 5. Sensitivity of the STIC derived fluxes ($\lambda E$ and $H$) to (a) actual uncertainties in MODIS Terra $T_S$, (b) the Priestley-Taylor parameter and (c) uncertainties in $e_A$ over four different biomes. Actual $T_S$ uncertainties over these biomes were obtained from Hulley et al. (2012). One representative site from every biome is chosen. The numbers of data points were 230 for cropland, 138 for grassland, 138 for forest, and 46 for savanna. Average sensitivity of all data points is reported.*

### 4.2. Evaluation of $\lambda E$ and $H$ estimates at MODIS Terra and Aqua overpass time
The predicted versus measured $\lambda E$ and $H$ show good overall agreement in all the 30 tower sites (Table 3) with a relative cumulative error $(\text{observed} - \text{estimated})/\text{observed}$ in $\lambda E$ and $H$ of the order of 4% and 3% respectively. The root mean square error (RMSE) of individual sites varied between 10.70 and 59.67 $\text{W}\ \text{m}^{-2}$ with measured versus predicted correlation coefficients ($r$) between 0.55 and 0.97 (Table 3). The RMSE and $r$ of $\lambda E$ from all the sites were 37.79 $\text{W}\ \text{m}^{-2}$ (11% error) and 0.89 (Table 3). The performance of STIC was consistent over all sites except FOG and DEM where the RMSE was relatively higher (see Table 3). Pooled RMSE and $r$ for $H$ for the 30 sites were 37.74 $\text{W}\ \text{m}^{-2}$ (9% error) and 0.91. RMSE's of individual sites varied between 12.11 and 58.48 $\text{W}\ \text{m}^{-2}$ with $r$ varying between 0.55 and 0.96 (Table 3).

The pooled evaluation of $\lambda E$ and $H$ using all the tower sites revealed the credible estimates of STIC fluxes over a broad range of biomes and climatic settings. An even distribution of points around 1:1 validation line (Fig. 6a and b) indicates low bias, having a measured to predicted regression slopes of 0.99 ($\pm0.005$) and 0.88 ($\pm0.003$) respectively. Given the footprint size of one MODIS pixel is 1 km it was not possible to obtain the $T_S$ heterogeneity within the tower footprint. However the $2 \times 2$ cut-out of vegetation index (i.e., NDVI) (250 m spatial resolution) around the tower sites revealed the standard deviation of NDVI to be varying from 0.01 to 0.10. Therefore, we assumed these comparisons to be reasonable although errors may be propagated due to $T_S$ heterogeneity which is governed by the surface dryness-wetness patterns of the landscape around the EC towers, a feature not captured in NDVI data.

Consistent results are also obtained with MODIS Aqua data where the overall RMSE and $r$ of the predicted $\lambda E$ was 42.27 $\text{W}\ \text{m}^{-2}$ (15% error) and 0.88 (Fig. 6c). The RMSE of individual sites varied between 26.32 and 63.77 $\text{W}\ \text{m}^{-2}$ with $r$ values between 0.48 and 0.90. The pooled RMSE and $r$ for $H$ was 44.31 $\text{W}\ \text{m}^{-2}$ (8% error) and 0.90 (Fig. 6d). RMSE of individual sites varied between 23.67 and 66.65 $\text{W}\ \text{m}^{-2}$ with $r$ varying between 0.47 and 0.97. Fig. 6c and d showed the pooled evaluations which again indicate reasonable fits between predictions and measurements for $\lambda E$ and $H$, having regression slopes of 0.93 ($\pm0.006$) and 0.93 ($\pm0.005$) respectively. Given MODIS Terra $T_S$ is frequently used in the terrestrial studies (Wan, Zhang, Zhang & Li, 2004), further analyses of STIC were performed with the MODIS Terra $T_S$ only.

Fig. 7 shows the time series comparisons between STIC fluxes and their measured FLUXNET counterparts for six different climatic groups. These reveal that the temporal dynamics of $\lambda E$ and $H$ are consistently captured by STIC throughout the entire study period. In SP1 (Fig. 7b), relatively less seasonality was found in both measured and predicted $\lambda E$, and almost no seasonality was found in $H$. This is because SP1 is a tropical site with an annual rainfall of 850-1100 mm evenly distributed throughout the year and $T_S$ also showed less seasonality (with a coefficient of variation of 0.13). This combination might be responsible for lesser seasonal component of the surface fluxes in SP1.

The residual $\lambda E$ error (predicted-observed) was negatively correlated ($r = -0.47$) with the observed $\lambda E$ (Fig. 8a). Fig. 8a also highlights slight overestimation for the lower range of $\lambda E$ ($<50\ \text{W}\ \text{m}^{-2}$), and improvements with increasing $\lambda E$. This was further confirmed when the residual error was found to be negatively correlated ($r = -0.30$) with the measured surface soil moisture ($\theta$) (Fig. 8b) and $R_N/D_A$ ratio ($r = -0.31$) (Fig. 8c) for their lower range values (for $\theta = 0\text{-}10\ \text{m}^3\ \text{m}^{-3}$ and $R_N/D_A = 0\text{-}25$). The residual error was also positively correlated to $T_S$ ($r = 0.13$) and $D_A$ ($r = 0.31$) for a $D_A$ range of 20-50 hPa (Fig. 8d, e and inset of 8e). The residual error was evenly distributed across the entire range of $T_A$ (Fig. 8f). Like $\lambda E$, the residual errors in $H$ were also found to be influenced by the same variables, but in the opposite direction.

---

### Table 3
**Error statistics of STIC derived eight-day average of MODIS Terra equatorial crossing time $\lambda E$ and $H$ over the eddy covariance sites of FLUXNET.**

| Site | N | $\lambda E$ RMSE ($\text{W}\ \text{m}^{-2}$) | $\lambda E$ MB ($\text{W}\ \text{m}^{-2}$) | $\lambda E$ r | Mean Obs $\lambda E$ | Mean STIC $\lambda E$ | $\lambda E$ Err (%) | $H$ RMSE ($\text{W}\ \text{m}^{-2}$) | $H$ MB ($\text{W}\ \text{m}^{-2}$) | $H$ r | Mean Obs $H$ | Mean STIC $H$ | $H$ Err (%) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| AUD | 230 | 33.05 | 13.66 | 0.77 | 60.44 | 70.73 | 17 | 35.32 | -13.02 | 0.87 | 183.83 | 169.15 | 8 |
| FR2 | 138 | 33.02 | 27.72 | 0.90 | 122.26 | 153.68 | 26 | 39.07 | -25.72 | 0.85 | 163.93 | 140.31 | 14 |
| SRM | 138 | 21.23 | 3.44 | 0.91 | 56.16 | 59.60 | 6 | 24.16 | -13.63 | 0.95 | 184.73 | 171.10 | 7 |
| SO2 | 138 | 33.86 | 9.26 | 0.55 | 81.91 | 91.19 | 11 | 40.70 | -21.40 | 0.89 | 225.56 | 204.15 | 9 |
| TON | 230 | 40.63 | 13.32 | 0.71 | 72.11 | 86.32 | 20 | 46.84 | -31.73 | 0.95 | 186.60 | 158.36 | 15 |
| WKG | 138 | 23.34 | 2.83 | 0.81 | 44.28 | 45.30 | 2 | 24.01 | -2.52 | 0.93 | 161.03 | 158.00 | 2 |
| FMF | 92 | 34.88 | 2.27 | 0.85 | 97.87 | 100.14 | 2 | 35.26 | -14.13 | 0.92 | 219.94 | 205.81 | 6 |
| FUF | 92 | 33.54 | 5.73 | 0.83 | 82.03 | 87.76 | 7 | 36.76 | -16.42 | 0.92 | 231.91 | 215.48 | 7 |
| FWF | 92 | 20.18 | 3.88 | 0.90 | 67.21 | 71.10 | 6 | 24.43 | -13.19 | 0.96 | 168.40 | 155.21 | 8 |
| SP1 | 92 | 39.39 | 6.29 | 0.92 | 198.05 | 203.68 | 3 | 39.28 | -8.09 | 0.72 | 148.60 | 138.41 | 7 |
| HOW | 230 | 45.76 | 0.34 | 0.79 | 199.01 | 197.58 | 1 | 41.87 | -15.29 | 0.75 | 139.06 | 126.36 | 9 |
| FOG | 92 | 59.67 | 4.96 | 0.79 | 242.03 | 250.67 | 4 | 55.95 | 3.92 | 0.65 | 91.65 | 89.87 | 2 |
| LMA | 138 | 31.32 | 0.09 | 0.89 | 95.71 | 96.37 | 1 | 31.67 | 0.10 | 0.93 | 132.46 | 134.78 | 2 |
| DEM | 138 | 52.45 | 19.11 | 0.85 | 99.21 | 117.68 | 19 | 58.48 | -35.00 | 0.55 | 204.89 | 169.10 | 17 |
| MA1 | 92 | 44.58 | 35.34 | 0.83 | 86.93 | 126.37 | 45 | 52.47 | -42.41 | 0.84 | 191.25 | 144.45 | 24 |
| KRU | 138 | 32.96 | 13.77 | 0.75 | 68.98 | 78.61 | 14 | 34.23 | -15.56 | 0.82 | 216.84 | 202.75 | 7 |
| WCR | 276 | 44.24 | 25.23 | 0.90 | 74.22 | 89.15 | 20 | 44.17 | -25.15 | 0.79 | 122.80 | 105.90 | 14 |
| HAI | 230 | 34.88 | 17.11 | 0.89 | 64.29 | 79.12 | 23 | 35.06 | -17.44 | 0.89 | 97.49 | 86.39 | 11 |
| UMB | 138 | 31.09 | 11.88 | 0.97 | 99.61 | 114.63 | 15 | 31.19 | -11.89 | 0.90 | 107.50 | 109.61 | 2 |
| HES | 230 | 38.81 | 20.69 | 0.89 | 62.93 | 85.96 | 37 | 40.76 | -23.64 | 0.86 | 97.34 | 80.02 | 18 |
| CAI | 184 | 31.20 | -24.27 | 0.96 | 107.30 | 94.56 | 12 | 32.36 | 24.92 | 0.86 | 55.20 | 74.73 | 35 |
| WET | 230 | 42.85 | 16.93 | 0.89 | 65.22 | 85.77 | 31 | 42.85 | -16.93 | 0.91 | 107.69 | 96.04 | 11 |
| MEH | 138 | 25.49 | 14.64 | 0.96 | 61.13 | 77.70 | 27 | 26.96 | -18.14 | 0.95 | 76.91 | 60.64 | 21 |
| GRI | 92 | 19.21 | -11.87 | 0.97 | 70.47 | 58.24 | 17 | 18.99 | 11.68 | 0.85 | 39.12 | 41.01 | 5 |
| OAS | 230 | 32.85 | 12.57 | 0.91 | 54.21 | 67.90 | 25 | 40.71 | -19.62 | 0.91 | 105.32 | 87.97 | 16 |
| LET | 230 | 32.34 | 4.39 | 0.90 | 61.37 | 67.87 | 10 | 33.56 | -6.29 | 0.91 | 107.33 | 109.29 | 2 |
| LOO | 276 | 39.04 | 8.41 | 0.86 | 87.69 | 105.03 | 20 | 40.29 | -12.77 | 0.87 | 108.02 | 90.20 | 17 |
| SF1 | 92 | 31.90 | 11.62 | 0.94 | 74.37 | 84.52 | 14 | 32.40 | -12.31 | 0.96 | 120.79 | 109.79 | 9 |
| TP1 | 92 | 45.96 | -35.58 | 0.96 | 152.90 | 117.43 | 23 | 45.13 | 34.64 | 0.85 | 74.09 | 113.65 | 53 |
| IVO | 138 | 10.70 | -2.82 | 0.94 | 16.53 | 13.49 | 18 | 12.11 | 0.82 | 0.95 | 18.59 | 18.99 | 2 |
| Pooled | 4646 | 37.79 | 10.10 | 0.89 | 93.20 | 103.31 | 11 | 37.74 | -12.18 | 0.91 | 138.52 | 126.34 | 9 |

---

![Fig. 6. Comparison of measured versus STIC estimates.](https://placeholder-for-figure6)  
*Fig. 6. Comparison of measured versus STIC estimates of (a) $\lambda E$ using MODIS Terra $T_S$ and tower meteorology over the 30 EC sites of FLUXNET with $\lambda E_{\text{predicted}} = 0.99 (\pm0.005) E_{\text{measured}}$, (b) $H$ using MODIS Terra $T_S$ and tower meteorology over the 30 EC sites of FLUXNET with $H_{\text{predicted}} = 0.88 (\pm0.003) H_{\text{measured}}$, (c) $\lambda E$ using MODIS Aqua $T_S$ and tower meteorology over the 28 EC sites of FLUXNET with $\lambda E_{\text{predicted}} = 0.93 (\pm0.006) E_{\text{measured}}$, (d) $H$ using MODIS Aqua $T_S$ and tower meteorology over the 28 EC sites of FLUXNET with $H_{\text{predicted}} = 0.93 (\pm0.005) H_{\text{measured}}$.*

---

### 4.3. Evaluation of daily $\lambda E$ and $H$ estimates using MODIS Terra
Daily $\lambda E$ was computed following the principle of diurnal conservation of evaporative fraction following Anderson et al. (2007). The estimated evaporative fraction $(\Lambda = \lambda E / (\lambda E + H))$ from STIC was multiplied by the daily net available energy $(\Phi = R_N - G)$ which was available at each EC site. Thus daily $\lambda E = \Lambda \Phi$ and daily $H = (1 - \Lambda)\Phi$. Analysis of the daily predicted versus measured $\lambda E$ using MODIS Terra revealed RMSE to be 5.59-34.89 $\text{W}\ \text{m}^{-2}$ with $r$ varying between 0.61 and 0.98 (Table 4). Overall RMSE and $r$ was 15.70 $\text{W}\ \text{m}^{-2}$ and 0.89, with percent error ranging from 1 to 31% (Table 4). Pooled RMSE and $r$ for $H$ over all the sites were 14.52 $\text{W}\ \text{m}^{-2}$ (6% error) and 0.88, with a range varying from 5.64-25.72 $\text{W}\ \text{m}^{-2}$ and 0.48-0.96, respectively.

### 4.4. Evaluation of annual $\lambda E$
Given the immediate connection between ecosystem water use and water resource management (Anderson et al., 2012), annual evapotranspiration ($E$) over the flux towers was also evaluated. Annual $E$ was computed by summing the eight-day $\lambda E$ values. If in any year, observed or estimated $E$ value in an eight-day block was missing, that particular year was not included in the computation. If such gaps exist for all the years in any of the study sites, that particular site was not accounted for in this analysis. The analysis revealed a good agreement for all sites (Fig. 9a) with $r = 0.95$ and RMSE of 41 mm $\text{yr}^{-1}$. The magnitude of $E_T$ had a wide range (Fig. 9b) between 114 and 720 mm, with a minimum $E_T$ over the WKG and maximum $E_T$ over FOG. STIC successfully captured this range.

### 4.5. Surface moisture availability (M) and partitioned $\lambda E$
Fig. 10(a and b) shows the surface moisture availability ($M$) along with the partitioned $\lambda E$ ($\lambda E_E$ and $\lambda E_T$) for representative grassland (WKG) and savanna (HOW) sites having contrasting annual rainfall ($R_F$). This reveals the higher proportion of $\lambda E_E$ (52% of total $\lambda E$) in HOW as compared to WKG (34% of total $\lambda E$). This is because the annual rainfall in HOW is nearly three times higher than the annual rainfall in WKG. This might have led to high surface wetness in HOW and as a result $\lambda E_E$ might have dominated. Seasonality was evident in $M$, $\lambda E_E$ and $\lambda E_T$ dynamics. $\lambda E_T$ was dominant particularly during the growing season while $\lambda E_E$ dominated immediately after the rainfall. However, there tend to be some $\lambda E_E$ signals even when there was no rainfall (Fig. 10). Given the MODIS Terra overpass time, the morning dew might have contributed to some degree of $\lambda E_E$ in the absence of rainfall. Dew is an important source of moisture in the semi-arid regions and it surrounds the surface nearly every morning for approximately two to three hours past sunrise (Baier, 1966; Malek, McCurdy & Giles, 1999). The weather data of most of the sites confirmed the differences between midnight to early morning $T_A$ and $T_D$ to be quite narrow, which might have developed a favorable environment for the dew formation.

Over the entire time period, $\Sigma \lambda E_E/\Sigma \lambda E$ varied from 33 to 63% across all the sites (Table 5). Highest proportion (63%) of $\lambda E_E$ was found in the SP1 which also had the highest $R_F$ (1420 mm) during the entire study period. The proportion of $\Sigma \lambda E_T/\Sigma \lambda E$ varied between 40 and 67% for all the sites (Table 5).

Significant correlation ($r = 0.31\text{-}0.88$) was found between $M$ and the observed surface soil moisture ($\theta$) across all the grassland and savanna sites and very low or no correlation ($r = 0.0\text{-}0.18$) was found in the forest sites (Table 5). Examples of the scatterplots between $M$ and $\theta$ (Fig. 10c and d) over the representative grassland and savanna sites showed the reliability of retrieved $M$ to reasonably capture the dynamics of $\theta$ over the landscapes having low or sparse vegetation cover. Table 5 also shows significant cross correlation ($r = 0.33\text{-}0.62$) between $M$ and eight-day cumulative $R_F$ for the grasslands and savannas and no or poor correlation ($r = 0.0\text{-}0.17$) over the forests.

The retrieved $\lambda E_E$ was significantly related to $\theta$ over almost every grassland and savanna site with $r$ varying from 0.24 to 0.79 (Table 5). No relationship between $\lambda E_E$ and $\theta$ was found in the temperate and boreal forest sites. Examples of the scatterplots (Fig. 10e and f) between $\lambda E_E$ and $\theta$ over the representative grassland and savanna sites confirmed the dominance of $\lambda E_E$ under higher surface wetness conditions. The effects of atmospheric water demand and annual $R_F$ were reflected in the relationships between $M$, $\theta$, and $\lambda E_E$. Regions having high $D_A$ and low annual $R_F$ (<500 mm) generally showed higher correlation between $M$, $\theta$, and $\lambda E_E$ (Table 5). An exception was HOW, where the correlation was substantially high despite receiving very high annual $R_F$ (1450 mm). This may be due to the existence of high atmospheric 'drying power' and strong land atmosphere interaction hot-spot over Australia that helps maintaining such close relationship between $M$, $\theta$, and $\lambda E_E$ (Koster et al., 2006), despite receiving high annual $R_F$.

Given the close link between transpiration and root zone soil moisture ($\theta_d$) (Anderson et al., 2007; Guswa, 2005; Noilhan & Planton, 1989; Small & McConnell, 2008), the retrieved $\lambda E_T$ was also linked with the tower measured $\theta_d$ (wherever available). A significant relationship ($r = 0.30\text{-}0.67$) (Table 5) was found in the semi-arid grasslands and savannas under the conditions of moderate to high atmospheric water demand ($D_A \ge 10\ \text{hPa}$) (Fig. 11a and inset). But no relationship between $\lambda E_T$ and $\theta_d$ was apparent under low $D_A$ ($D_A < 10\ \text{hPa}$) despite the availability of unlimited radiation ($R_G$) and high leaf area index ($L$) (Fig. 11a and c). A direct response of $\lambda E_T$ to $R_G$ was also evident ($r = 0.70\text{-}0.90$) when $\theta_d$ exceeded certain threshold limit ($>15\ \text{m}^3\ \text{m}^{-3}$) (Fig. 11b). No relationship between $\lambda E_T$ and $\theta_d$ was found in the high latitude forests and temperate grasslands due to the low humidity deficit of air.

---

![Fig. 7. Illustrative examples of the time series comparison.](https://placeholder-for-figure7)  
*Fig. 7. Illustrative examples of the time series comparison of $\lambda E$ and $H$ observations and STIC predictions using MODIS Terra $T_S$ and tower meteorology at selected climate types representing (a) dry & semi-arid (SRM), (b) tropical (SP1), (c) subtropical Mediterranean (LMA), (d) temperate (MEH), (e) boreal (SF1), and (f) arctic (IVO) climate group according to Köppen's climate classification. Any gap in the time series is caused either due to the absence of flux measurements or missing input data.*

---

## 5. Discussion
This study has shown that the surface flux estimates from STIC are able to capture the observed dynamics of the fluxes over diverse biomes and climate types and provide reasonable estimates of $\lambda E$ (and $H$) and its components. A systematic overestimation of $\lambda E$ is noted when very low surface soil moisture is combined with a very high radiative load and atmospheric water demand (see Fig. 8b and c). This is typically found in the semi-arid regions where $\lambda E$ is influenced by the tight coupling between $R_N$, $D_A$ and soil moisture (Gu, Meyers, Pallardy, Hanson, Yang, Heuer, et al., 2006; Small & Kurc, 2003). These might also be the reasons for the relatively higher RMSE over FOG and DEM where such tight coupling tends to dominate $\lambda E$. The present moisture availability retrieval scheme (Eq. 5) does not explicitly consider this coupling, thus causing an overestimation of $\lambda E$ under these conditions. Granger & Gray (1989) have already described a moisture availability function that accounts for the influence of $R_N$, $D_A$ and $\lambda E$ on $M$. Therefore, a hybrid approach of $M$ retrieval may be more useful in this context.

Cloudiness and the associated biased sampling of $T_S$ by the satellite will also affect the performance of STIC. The eight-day $T_S$ products of MODIS are, by definition, a sample of relatively cloud free conditions whilst the tower fluxes are mixtures of clear and cloudy atmospheric conditions. The mismatch in the sampling footprint between MODIS $T_S$ and the flux tower meteorology-micrometeorology we have used will also contribute to the uncertainty of the analysis. Another source of error will be the presence of variable dry-wet patches within a MODIS pixel as well as around the flux towers (LeMone, Tewari, Chen, Alfieri, & Niyogi, 2008). For example, if more than 50% of the area falling within a 1 km MODIS pixel is predominately wet (or dry), the lumped MODIS $T_S$ signal will be biased due to the wetness (or dryness) of the landscape. The degree to which this affects the performance of STIC needs to be assessed, although it is surprising how well the scheme has performed given little is known about within footprint heterogeneity.

The performance of STIC in comparison to other approaches cannot be directly assessed without running the other models with similar datasets, which is beyond the scope of this paper. But some inferences can be drawn by comparing and contrasting the results with other studies that also used $T_S$ observations for computing surface energy fluxes in a single-source or two-source framework. Using hourly measurements of $T_S$ and associated meteorological-micrometeorological variables, Norman et al. (1995) reported RMSE in $\lambda E$ of 60 and 50 $\text{W}\ \text{m}^{-2}$ with Monsoon'90 and FIFE experimental data, while the RMSE in $H$ was 35 (Monsoon'90) and 50 $\text{W}\ \text{m}^{-2}$ (FIFE), respectively. Using a single-source surface energy balance model (SEBS), Su (2002) reported RMSE of 61.34-82.79 $\text{W}\ \text{m}^{-2}$ (in $\lambda E$) and 28.61-36.19 $\text{W}\ \text{m}^{-2}$ (in $H$) over semi-arid shrub and grasses. Interestingly, in all these previous studies the prediction accuracy of $H$ was better than that for $\lambda E$. We note that STIC has a lower RMSE for both $\lambda E$ and $H$ even when applied to a much broader class of land surface types and with the advantage of being independent of any aerodynamic inputs (e.g., wind speed, vegetation height and surface roughness) that are typically required to model $H$. The use of $T_S$ and the temperature-saturation vapor pressure slopes to estimate the near surface moisture and vapor pressure (Fig. 2) provided the information on lower boundary conditions for $\lambda E$ and $H$. Information on surface roughness are also implicit in $T_S$ (high roughness will cause $T_S$ to approach $T_A$), $R_N$ (through albedo and surface emissivity) and $G$ measurements that are direct inputs into STIC.

The spatial scale mismatch between satellite and tower produced varying relationship between $\lambda E_E$ versus $\theta$, $M$ versus $\theta$, and $\lambda E_T$ versus $\theta_d$ among the grasslands and savannas. $\theta$ and $\theta_d$ are the point measurements whereas the satellite retrieved $M$ are the mixed signals of the heterogeneous wetness. Over the forests, $T_S$ signals are mainly contributed by the vegetation and less by the soil. Therefore, the information of canopy wetness and interception evaporation was dominant in the retrieved $M$ and $\lambda E_E$ (Cavanaugh et al., 2011). This led to the poor relationship between $M$, $\theta$, and $\lambda E_E$ in the forest landscapes. The proportion of $\lambda E_E$ appears to be little higher than $\lambda E_T$ since it includes contributions from both the canopy and soil water. However, these fractions are computed at a single snapshot during the MODIS Terra equatorial crossing time and values are likely to vary at the daily time scale. Studies over Mediterranean savanna have shown the proportion of $\lambda E_E$ at the daily scale to be around 25 to 28% (Czikowsky & Fitzjarrald, 2009; David et al., 2006).

The nature of relationship between $\lambda E_T$ and $\theta_d$ (Table 5 and Fig. 11a) supports the hypothesis that transpiration is determined by two limiting conditions in the water-limited ecosystems: (a) uptake limited by available energy when $\theta_d$ is plentiful (Guswa, 2005), (b) by atmospheric $D_A$ and available $\theta_d$ under water-stressed conditions (Denmead & Shaw, 1962). High $D_A$ values result in higher soil-water content thresholds at which transpiration rates are observed to decline. The degree and the direction of the relationship between $\lambda E_T$ and shortwave radiation ($R_G$) also confirmed the fact that the energy content of the radiation absorbed by vegetation influences the transpiration (Pieruschka et al., 2010) when the level of $\theta_d$ crossed a certain limit. Even in the presence of substantial moisture in the root zone and high leaf area, $\lambda E_T$ may still be suppressed due to insufficient radiative load (Small & McConnell, 2008), as seen in Fig. 11b and c.

---

## 6. Conclusions and future research
The analytical method presented here addressed some of the stumbling blocks that previously hindered the direct use of radiometric surface temperature into the PM equation. This also demonstrated a physical integration of $T_S$ and the PM equation to derive a hybrid closure that does not require the specification of surface to atmosphere conductance terms. An initial evaluation of STIC over a variety of biomes and climate types produced encouraging results. Errors between the predicted and observed fluxes appeared reasonable given the uncertainties in the assumptions being made and the data used. The advantage of STIC is that the Priestley-Taylor parameter, $\alpha$, is the only land surface parameterization required and STIC is relatively insensitive to the assumed value of this parameter.

The modest number of input variables required to run STIC suggests that this framework would be suitable for generating spatially explicit latent and sensible heat fields using readily available data from the current generation remote sensing platforms (e.g., MODIS and AIRS). We would also argue that because the surface conductances we produce are non-parametric byproducts of STIC, these estimates could serve as a means of developing and testing land surface parameterizations embedded within the climate and Earth system models.

We envisage further development of STIC. For example, alternative expressions to recover the surface moisture availability from the information of early morning rise in $T_S$ (Anderson et al., 1997, 2007) or by combining $R_N$, $D_A$ and $T_S$ (Granger & Gray, 1989) needs to be explored particularly during the dry down phase. For the surface energy balance mapping with STIC at the regional and continental scale, the challenge is to obtain $R_N$, $G$, $T_A$ and $R_H$ or $e_A$ that are compatible with the radiometric surface temperature data. The advent of improved thermal remote sensing data through future missions like HyspIRI, GeoSTAR and GOES-R may afford an opportunity to extend this method across multiple spatio-temporal scales, allowing for more spatially explicit hydrological and physiological process studies.

## Acknowledgments
K.M. acknowledges the postdoctoral fellowship from Jet Propulsion Laboratory, California Institute of Technology under the contract of National Aeronautics and Space Administration. K.M. and A.J. was earlier supported by NERC grant (NE/E019153/1) during the initial development of the work. The research was initiated at Lancaster University and completed at the Jet Propulsion Laboratory. J.B.F. was supported by JPL Research and Technology Development Grant (01STCR-R.11.118.004) and Earth System Data Record (ESDR) grant (104772-547714.04.16.01.06). We gratefully acknowledge entire FLUXNET site PIs for sharing the eddy covariance data and ORNL for maintaining the harmonized MODIS data. K.M. also acknowledges Dr. Junhak Lee for helping in GIS. The scientific discussions with Dr. Jozsef Szilagyi, University of Nebraska and Dr. Wilfried Brutsaert, Cornell University are also acknowledged. D.N. acknowledges partial support from NSF CAREER, NSF INTEROP, and USDA NIFA/U2U. All the copyrights of 2013 are reserved. K.M., A.J., and E.B. formulated idea, designed research and performed research; everyone contributed to writing the manuscript. The authors declare no conflict of interest.

---

### Table 4
**Error statistics of STIC derived daily average $\lambda E$ and $H$ over the eddy covariance sites of FLUXNET (No daily data was available for DEM).**

| Site | N | $\lambda E$ RMSE ($\text{W}\ \text{m}^{-2}$) | $\lambda E$ MBE ($\text{W}\ \text{m}^{-2}$) | $\lambda E$ r | Mean Obs $\lambda E$ | Mean STIC $\lambda E$ | $\lambda E$ Err (%) | $H$ RMSE ($\text{W}\ \text{m}^{-2}$) | $H$ MBE ($\text{W}\ \text{m}^{-2}$) | $H$ r | Mean Obs $H$ | Mean STIC $H$ | $H$ Err (%) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| AUD | 230 | 11.46 | 0.33 | 0.84 | 21.01 | 22.37 | 6 | 12.65 | -0.65 | 0.89 | 51.43 | 50.93 | 1 |
| FR2 | 138 | 10.63 | 0.12 | 0.95 | 49.47 | 50.96 | 3 | 10.03 | -0.38 | 0.81 | 43.90 | 44.94 | 2 |
| SRM | 138 | 10.39 | -1.02 | 0.91 | 22.96 | 22.04 | 4 | 11.19 | 0.98 | 0.94 | 62.92 | 64.26 | 2 |
| SO2 | 138 | 12.71 | 0.03 | 0.61 | 30.96 | 31.12 | 1 | 12.58 | -0.21 | 0.93 | 71.10 | 71.36 | 0.3 |
| TON | 230 | 16.54 | 0.10 | 0.74 | 29.58 | 31.00 | 5 | 16.70 | -0.62 | 0.92 | 58.38 | 59.87 | 3 |
| WKG | 138 | 9.73 | 0.30 | 0.83 | 17.54 | 16.95 | 3 | 10.03 | -0.35 | 0.93 | 60.55 | 59.95 | 1 |
| FMF | 92 | 13.31 | -5.84 | 0.87 | 35.25 | 29.40 | 17 | 13.89 | 6.03 | 0.93 | 56.78 | 62.62 | 10 |
| FUF | 92 | 13.23 | -2.59 | 0.83 | 30.68 | 28.08 | 8 | 14.33 | 4.03 | 0.94 | 67.27 | 69.87 | 4 |
| FWF | 92 | 10.45 | -5.69 | 0.89 | 26.27 | 20.57 | 22 | 10.02 | 5.35 | 0.96 | 40.24 | 45.94 | 14 |
| SP1 | 92 | 14.66 | -7.27 | 0.95 | 74.77 | 67.16 | 10 | 14.97 | 7.81 | 0.55 | 27.22 | 34.46 | 26 |
| HOW | 230 | 17.13 | 0.35 | 0.85 | 85.29 | 85.47 | 0.02 | 16.75 | -0.43 | 0.67 | 49.29 | 47.10 | 4 |
| FOG | 92 | 26.80 | -8.26 | 0.70 | 114.26 | 111.11 | 3 | 25.72 | 11.50 | 0.48 | 32.14 | 38.44 | 20 |
| LMA | 138 | 16.25 | -7.87 | 0.87 | 40.48 | 32.60 | 19 | 16.43 | 7.32 | 0.89 | 38.34 | 47.43 | 24 |
| MA1 | 92 | 13.15 | 8.08 | 0.88 | 30.89 | 40.64 | 31 | 14.77 | -9.34 | 0.79 | 55.68 | 45.71 | 18 |
| KRU | 138 | 11.44 | 1.40 | 0.82 | 22.57 | 22.68 | 0.5 | 12.76 | -1.36 | 0.87 | 61.05 | 57.54 | 6 |
| WCR | 276 | 15.47 | -1.11 | 0.89 | 33.41 | 32.23 | 4 | 15.12 | 0.07 | 0.70 | 24.61 | 27.35 | 11 |
| HAI | 230 | 10.17 | -0.69 | 0.86 | 24.59 | 22.69 | 8 | 10.23 | 0.21 | 0.86 | 12.65 | 17.71 | 40 |
| UMB | 138 | 13.55 | -4.16 | 0.94 | 50.58 | 46.42 | 8 | 12.79 | 3.02 | 0.88 | 30.71 | 37.55 | 22 |
| HES | 230 | 12.47 | 0.10 | 0.87 | 21.58 | 22.77 | 6 | 12.54 | -1.44 | 0.69 | 25.33 | 17.91 | 29 |
| CAI | 184 | 19.04 | -17.16 | 0.94 | 40.25 | 27.67 | 31 | 19.40 | 17.31 | 0.62 | 5.77 | 20.41 | 253 |
| WET | 230 | 14.69 | -2.53 | 0.82 | 25.97 | 24.89 | 4 | 14.76 | 1.47 | 0.80 | 15.47 | 25.09 | 62 |
| MEH | 138 | 5.59 | -0.74 | 0.98 | 23.29 | 23.14 | 0.6 | 5.64 | 0.67 | 0.90 | 12.21 | 18.30 | 49 |
| GRI | 92 | 10.59 | -8.63 | 0.96 | 24.59 | 18.74 | 24 | 10.52 | 8.27 | 0.61 | 3.59 | 12.50 | 248 |
| OAS | 230 | 16.59 | -0.09 | 0.84 | 23.57 | 23.85 | 1 | 16.45 | -0.22 | 0.85 | 36.10 | 30.77 | 15 |
| LET | 230 | 17.14 | -6.12 | 0.89 | 24.49 | 20.00 | 18 | 17.50 | 6.00 | 0.74 | 17.89 | 27.62 | 54 |
| LOO | 276 | 17.74 | -6.13 | 0.79 | 38.72 | 34.22 | 11 | 14.41 | 0.81 | 0.84 | 22.51 | 27.89 | 24 |
| SF1 | 92 | 11.80 | -1.79 | 0.92 | 33.81 | 31.77 | 6 | 11.84 | 1.44 | 0.93 | 30.33 | 37.52 | 24 |
| TP1 | 92 | 34.89 | -21.86 | 0.81 | 76.16 | 60.15 | 21 | 22.59 | 19.42 | 0.79 | 10.84 | 39.66 | 265 |
| IVO | 138 | 6.26 | -2.15 | 0.93 | 9.36 | 5.76 | 38 | 6.22 | 1.70 | 0.93 | 5.66 | 8.56 | 51 |
| Pooled | 4508 | 15.70 | -3.35 | 0.89 | 41.03 | 37.68 | 8 | 14.52 | 2.48 | 0.88 | 42.84 | 45.32 | 6 |

---

![Fig. 9. Annual E predictions validation.](https://placeholder-for-figure9)  
*Fig. 9. (a) Validation of STIC estimates of annual $E$ against the tower measurements. These are the annual sum of $E$ for multiple years at each of the flux tower sites. (b) Bar chart of mean annual measured $E$ and estimated $E$. Averaging all the individual year $E$ values were done over the eddy covariance sites. The number of years in this averaging varied from at least 1 to maximum 5.*

---

## Appendix A

The surface energy balance equation can be written as:
$$\Phi = H + \lambda E \qquad (A1)$$

where $\Phi =$ net available energy $(\cong R_N - G)$, $H =$ sensible heat flux, $\lambda E =$ latent heat flux, $R_N =$ net radiation, $G =$ conductive surface heat flux or ground heat flux. All the fluxes have units in $\text{W}\ \text{m}^{-2}$. The sensible and latent heat flux can be expressed in the form of aerodynamic transfer equations (Boegh & Soegaard, 2004; Boegh et al., 2002) as follows:
$$H = \rho c_P g_B \Delta T \qquad (A2)$$
$$\lambda E = \frac{\rho c_P}{\gamma} g_B (e_S - e_A) \qquad (A3)$$

In many studies $\Delta T$ is expressed as the difference between $T_S$ and $T_A$ $(\Delta T = T_S - T_A)$ and in such cases, an 'extra conductance' is introduced to compensate for errors arising due to the inequalities between $(T_{SA} - T_A)$ and $(T_{SR} - T_A)$ (Boegh et al., 2002; Su, 2002). But $T_S$ is not the true temperature that is responsible for transferring the sensible heat flux from surface to the atmosphere (Anderson et al., 1997; Norman et al., 1995; Troufleau et al., 1997). Despite the apparent simplicity of Eq. (A2), the main limitation lies in the definition of $T_{SA}$. $T_{SA}$ is theoretically an air temperature at the source/sink height, which is different from the physical temperature ($T_S$) of the surface (Monteith, 1965). It is the temperature of the thin boundary layer in the immediate vicinity of the surface level (Fig. 1) and is responsible for the transfer of heat from surface to the atmosphere. This level is defined as the source height where wind speed is zero and $T_{SA}$ is obtained by extrapolating the logarithmic profile of $T_A$ down to that level (Troufleau et al., 1997).

### Expression for $g_B$
By combining Eqs. (A1), (A2) and (A3) and solving for $g_B$ we get:
$$g_B = \frac{\Phi}{\rho c_P \left(\Delta T + \frac{e_S - e_A}{\gamma}\right)} \qquad (A4)$$

### Expression for $g_S$
According to Boegh et al. (2002), $\lambda E$ can also be expressed as:
$$\lambda E = \frac{\rho c_P}{\gamma} g_S (e_S^* - e_S) \qquad (A5)$$

Combining Eqs. (A3) and (A5) and solving for $g_S$ we can express $g_S$ in terms of $g_B$, $e_S^*$, $e_S$, and $e_A$:
$$g_S = g_B \frac{(e_S - e_A)}{(e_S^* - e_S)} \qquad (A6)$$

Water vapor transfer occurs from within the vegetation (transpiration) and from the immediate vicinity of the vegetation surface (evaporation). The stomatal cavities are assumed to be saturated with respect to water vapor, therefore, $e_S^*$ of dense canopies can always be estimated from $T_S$. For extremely dry bare soil, the evaporating front may sometimes be located little below the dry surface layer and expressing $e_S^*$ in terms of $T_S$ may produce some errors under such circumstances. Recognizing the power of $T_S$ to detect the signal of the surface as well as the subsurface dry-wet regimes (Anderson et al., 2008; Kustas & Anderson, 2009), $e_S^*$ was estimated from $T_S$ in the present case. Given the measurements or estimates of $e_A$, $e_S$, $\Phi$, $\rho$, $c_P$ and $\gamma$, there remains one additional unknown variable ($\Delta T$) in Eqs. (A4) and (A6). Here we have used the Bowen ratio ($\beta$) (Bowen, 1926) equation for expressing $\Delta T$.

### Expression for $\Delta T$
An expression for $\Delta T$ was derived from the Bowen ratio ($\beta$) equation (Bowen, 1926):
$$\beta = \frac{H}{\lambda E} = \gamma \frac{\Delta T}{e_S - e_A} \qquad (A7)$$

Assuming closure of the surface energy balance it is possible to express $\beta$ in terms of evaporative fraction ($\Lambda$) (Shuttleworth et al., 1989):
$$\beta = \frac{1 - \Lambda}{\Lambda} \qquad (A8)$$

$\Lambda$ is defined as the fraction of available energy ($\Phi$) partitioned towards $\lambda E$. Substituting for $\beta$ in Eq. (A8) we get an expression for $\Delta T$ in terms of $\Lambda$:
$$\Delta T = \left(\frac{e_S - e_A}{\gamma}\right)\left(\frac{1 - \Lambda}{\Lambda}\right) \qquad (A9)$$

While expressing $\Delta T$ another extra variable ($\Lambda$) is introduced in Eq. (A9). One more equation is needed to close the system of equations and this equation must introduce the dependence of $\Lambda$ on the conductances. In order to express $\Lambda$ in terms of $g_B$ and $g_S$ we had adopted the advection-aridity hypothesis (Brutsaert & Stricker, 1979).

### Expression for $\Lambda$
According to Brutsaert & Stricker (1979),
$$E_P^* = 2E_{PT}^* - E \qquad (A10)$$

dividing both sides by $E$ we get,
$$\frac{E_P^*}{E} = \frac{2E_{PT}^* - E}{E} \qquad (A11)$$

and dividing the numerator and denominator of the right hand side of Eq. (A11) by $E_{PT}^*$ we get,
$$\frac{E}{E_P^*} = \frac{\frac{E}{E_{PT}^*}}{2 - \frac{E}{E_{PT}^*}} \qquad (A12)$$

Expressing $\Phi$ as $\Phi = E / \Lambda$ gives,
$$\frac{E}{E_{PT}^*} = \frac{s + \gamma}{\alpha s} \Lambda \qquad (A13)$$

Now substituting $E/E_{PT}^*$ from Eq. (A13) into Eq. (A12) and after some algebra gives:
$$\frac{E}{E_P^*} = \frac{\Lambda(s + \gamma)}{2\alpha s - \Lambda(s + \gamma)} \qquad (A14)$$

According to the Penman equation (Penman, 1948) and PM equation (Monteith, 1965),
$$\frac{E}{E_P^*} = \frac{\frac{s\Phi + \rho c_P g_B D_A}{s + \gamma \left(1 + \frac{g_B}{g_S}\right)}}{\frac{s\Phi + \rho c_P g_B D_A}{s + \gamma}} = \frac{s + \gamma}{s + \gamma \left(1 + \frac{g_B}{g_S}\right)} \qquad (A15)$$

Equating Eqs. (A14) and (A15) gives an expression for $\Lambda$ in terms of the conductances,
$$\frac{\Lambda(s + \gamma)}{2\alpha s - \Lambda(s + \gamma)} = \frac{s + \gamma}{s + \gamma \left(1 + \frac{g_B}{g_S}\right)} \qquad (A16)$$

which, after some algebra, gives the final expression of $\Lambda$ in terms of $g_B$ and $g_S$:
$$\Lambda = \frac{2\alpha s}{2s + \gamma \left(2 + \frac{g_B}{g_S}\right)} \qquad (A17)$$

---

![Fig. 10. Partitioned E components and soil moisture scatter.](https://placeholder-for-figure10)  
*Fig. 10. (a) Estimates of the evaporation ($\lambda E_E$) (thick black line) and transpiration ($\lambda E_T$) (thick grey line) components of $\lambda E$ along with the surface moisture availability ($M$) (black dots) and rainfall ($R_F$) (bar) for representative sites of grassland (WKG) and (b) woody savanna (HOW). These two sites have contrasting annual $R_F$ pattern. (c) and (d) Two dimensional scatters between retrieved $M$ (black dots) with the tower measured upper layer soil moisture ($\theta$) for similar representative grassland (WKG) and woody savanna (HOW) sites. (e) and (f) Illustrative examples of the two dimensional scatter between retrieved $\lambda E_E$ and tower measured $\theta$ in grassland (WKG) and woody savanna (HOW).*

---

### Table 5
**Proportion of $\lambda E_E$ and $E_T$ from STIC over all the sites along with the cross correlation between $M$ vs. $\theta$, $\lambda E_E$ vs. $\theta$, $\lambda E_T$ vs. $\theta$, and $M$ vs. $R_F$ (Analysis done with MODIS Terra data only).**

| Biome | Site name | $\Sigma \lambda E_E / \Sigma \lambda E$ | $\Sigma \lambda E_T / \Sigma \lambda E$ | Correlation ($M$ vs. $\theta$) | Correlation ($M$ vs. $R_F$) | Correlation ($\lambda E_E$ vs. $\theta$) | Correlation ($\lambda E_T$ vs. $\theta_d$) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Grassland, shrubland, savanna** | AUD | 0.40 | 0.60 | 0.48 | 0.39 | 0.43 | 0.50 |
| | FR2 | 0.43 | 0.57 | 0.55 | 0.45 | 0.24 | n/a |
| | SRM | 0.34 | 0.66 | 0.76 | 0.47 | 0.69 | 0.33 |
| | SO2 | 0.47 | 0.53 | 0.72 | 0.60 | 0.39 | 0.30 |
| | TON | 0.43 | 0.57 | 0.76 | 0.53 | 0.58 | 0.67 |
| | WKG | 0.34 | 0.66 | 0.76 | 0.47 | 0.79 | 0.62 |
| | FMF | 0.50 | 0.50 | n/a | 0.48 | n/a | n/a |
| | FUF | 0.52 | 0.48 | n/a | 0.33 | n/a | n/a |
| | FWF | 0.50 | 0.50 | n/a | 0.30 | n/a | n/a |
| | SP1 | 0.63 | 0.37 | 0.31 | 0.47 | 0.38 | 0.54 |
| | HOW | 0.52 | 0.48 | 0.78 | 0.62 | 0.53 | n/a |
| | FOG | 0.35 | 0.65 | n/a | 0.41 | n/a | n/a |
| | LMA | 0.51 | 0.49 | 0.63 | 0.37 | 0.50 | 0.34 |
| | DEM | 0.33 | 0.67 | 0.88 | 0.37 | 0.51 | 0.66 |
| | MA1 | 0.45 | 0.55 | 0.58 | 0.38 | 0.53 | 0.45 |
| | KRU | 0.37 | 0.66 | 0.55 | 0.45 | 0.44 | 0.37 |
| | CAI | 0.58 | 0.42 | n/a | 0.18 | 0.0 | n/a |
| | MEH | 0.62 | 0.43 | 0.42 | 0.17 | 0.40 | 0.0 |
| | GRI | 0.61 | 0.43 | n/a | 0.15 | n/a | n/a |
| | LET | 0.44 | 0.60 | n/a | n/a | n/a | n/a |
| **Forests** | WCR | 0.55 | 0.47 | 0.0 | 0.0 | 0.11 | 0.06 |
| | HAI | 0.61 | 0.41 | 0.17 | 0.15 | 0.0 | 0.0 |
| | UMB | 0.51 | 0.49 | n/a | n/a | n/a | n/a |
| | HES | 0.56 | 0.44 | 0.11 | 0.12 | 0.0 | 0.0 |
| | WET | 0.60 | 0.40 | 0.14 | 0.16 | 0.0 | 0.0 |
| | OAS | 0.58 | 0.48 | 0.0 | 0.10 | 0.11 | 0.10 |
| | LOO | 0.57 | 0.44 | 0.18 | 0.17 | 0.0 | 0.0 |
| | SF1 | 0.63 | 0.44 | n/a | n/a | n/a | n/a |
| | TP1 | 0.55 | 0.51 | n/a | 0.14 | n/a | n/a |
| | IVO | 0.60 | 0.52 | n/a | n/a | n/a | n/a |

*n/a: no or spurious measurements of both $\theta$ and $\theta_d$ were available in these sites.*

---

![Fig. 11. Transpiration versus soil moisture and radiation plots.](https://placeholder-for-figure11)  
*Fig. 11. Examples of the scatterplot between (a) retrieved transpiration ($\lambda E_T$) versus measured root zone soil moisture ($\theta_d$) over semi-arid savanna for two different combinations of $D_A$, shortwave radiation ($R_G$) and leaf area index ($L$). The figure in the inset shows good correspondence between $\lambda E_T$ and $\theta_d$ when $D_A$ goes above 10 hPa, (b) $\lambda E_T$ versus $R_G$ for different levels of $\theta_d$, (c) the temporal dynamics of $\lambda E_T$ (blue line) along with the observed $\lambda E$ (black line), $R_G$ (green line), MODIS leaf area index ($L$) (thick red line) and $\theta_d$ (thick grey line) over TON. This shows that $\lambda E_T$ has a direct response to $R_G$ but $\lambda E_T$ does not respond immediately with increase in $\theta_d$ unless $D_A$ rises above some level.*

The advection-aridity hypothesis leads to an assumed link between $g_B$ and radiometric surface temperature ($T_S$). Therefore, the effects of advection and surface moisture are implicit (although not explicit) in the advection-aridity equation.

---

## References
* Amiro, B.D., Orchansky, A. L., Barr, A. G., Black, T. A., Chambers, S. D., Chapin, F. S., et al. (2006). The effect of post-fire stand age on the boreal forest energy balance. *Agricultural and Forest Meteorology*, 140, 41-50.
* Anderson, M. C., Allen, R. G., Morse, A., & Kustas, W. P. (2012). Use of Landsat thermal imagery in monitoring evapotranspiration and managing water resources. *Remote Sensing of Environment*, 122, 50-65.
* Anderson, M. C., Norman, J. M., Diak, G. R., Kustas, W. P., & Mecikalski, J. R. (1997). A two-source time-integrated model for estimating surface fluxes using thermal infrared remote sensing. *Remote Sensing of Environment*, 60, 195-216.
* Anderson, M. C., Norman, J. M., Kustas, W. P., Houborg, R., Starks, P. J., & Agam, N. (2008). A thermal-based remote sensing technique for routine mapping of land-surface carbon, water and energy fluxes from field to regional scales. *Remote Sensing of Environment*, 112, 4227-4241.
* Anderson, M. C., Norman, J. M., Mecikalski, J. R., Otkin, J. A., & Kustas, W. P. (2007). A climatological study of evapotranspiration and moisture stress across the continental United States based on thermal remote sensing. 1: Model formulation. *Journal of Geophysical Research*, 112(D10117), 1-17.
* Ardö, J., Mölder, M., El-Tahir, B.A., & Elkhidir, H. A.M. (2008). Seasonal variation of carbon fluxes in semi-arid Sudan. *Carbon Balance and Management*, 3(7).
* Baier, W. (1966). Studies on dew formation under semi-arid conditions. *Agricultural Meteorology*, 3, 103-112.
* Baldocchi, D.D., Falge, E., Gu, L., Olson, R., Hollinger, D., Running, S., et al. (2001). Fluxnet: a new tool to study the temporal and spatial variability of ecosystem-scale carbon dioxide, water vapor, and energy flux densities. *Bulletin of the American Meteorological Society*, 82(11), 2415-3434.
* Baldocchi, D.D., & Xu, L. (2007). What limits evaporation from Mediterranean oak woodlands-The supply of moisture in the soil, physiological control by plants or the demand by the atmosphere? *Advances in Water Resources*, 30, 2113-2122.
* Ball, J. T., Woodrow, I. E., & Berry, J. A. (1987). A model predicting stomatal conductance and its contribution to the control of photosynthesis under different environmental conditions. In J. Biggins, M. Nijhoff, & Dordrecht (Eds.), *Progress in photosynthesis research*, 4. (pp. 5.221-5.224) Dordrecht, The Netherlands: Martinus-Nijhoff Publishers.
* Barr, A. G., Morgenstern, K., Black, T. A., McCaughey, J. H., & Nesic, Z. (2006). Surface energy balance closure by the eddy covariance method above three boreal forest stands and implications for the measurement of the CO2 flux. *Agricultural and Forest Meteorology*, 140, 322-337.
* Bastiaansssen, W. G. M., Menenti, M., Feddes, R. A., & Holtslag, A. A.M. (1998). A remote sensing surface energy balance algorithm for land (SEBAL) 1. Formulation. *Journal of Hydrology*, 212-213, 198-212.
* Beringer, J., Hutley, L. B., Tapper, N. J., Coutts, A., Kerley, A., & O'Grady, A. P. (2003). Fire impacts on surface heat, moisture and carbon fluxes from a tropical savanna in north Australia. *International Journal of Wildland Fire*, 12, 333-340.
* Boegh, E., & Soegaard, H. (2004). Remote sensing based estimation of evapotranspiration rates. *International Journal of Remote Sensing*, 25(13), 2535-2551.
* Boegh, E., Soegaard, H., & Thomsen, A. (2002). Evaluating evapotranspiration rates and surface conditions using Landsat TM to estimate atmospheric resistance and surface resistance. *Remote Sensing of Environment*, 79, 329-343.
* Bosveld, F. C., & Bouten, W. (2003). Evaluating a model of evaporation and transpiration with observations in a partially wet Douglas-fir forest. *Boundary Layer Meteorology*, 108, 365-396.
* Bouchet, R. J. (1963). Evapotranspiration réelle evapotranspiration potentielle, signification climatique. *International Association of Scientific Hydrology, General Assembly of Berkeley*, Transactions, vol. 2. (pp. 134-142). Berkeley, California: Evaporation.
* Bowen, I. S. (1926). The ratio of heat losses by conduction and by evaporation from any water surface. *Physics Review*, 27, 779-787.
* Brutsaert, W., & Stricker, H. (1979). An advection-aridity approach to estimate actual regional evapotranspiration. *Water Resources Research*, 15(2), 443-450.
* Buck, A. L. (1981). New equations for computing vapor pressure and enhancement factor. *Journal of Applied Meteorology*, 20, 1527-1532.
* Budyko, M. I., Efimova, N. A., Zubenok, L. I., & Strokina, L. A. (1962). The heat balance of the earth's surface. *Izvestya Akademii Nauk SSSR, Seriya Geograficheskaya*, 1, 6-16.
* Cavanaugh, M. L., Kurc, S. A., & Scott, R. L. (2011). Evapotranspiration partitioning in semiarid shrubland ecosystems: A two-site evaluation of soil moisture control on transpiration. *Ecohydrology*, 4(5), 671-681.
* Chavez, J. L., Neale, C. M. U., Hipps, L. E., Prueger, J. H., & Kustas, W. P. (2005). Comparing aircraft-based remotely sensed energy balance fluxes with eddy covariance tower data using heat flux source area functions. *Journal of Hydrometeorology*, 6, 923-940.
* Choudhury, B. J. (1997). Global pattern of potential evaporation calculated from the Penman-Monteith equation using satellite and assimilated data. *Remote Sensing of Environment*, 61, 64-81.
* Cleugh, H. A., Leuning, R., Mu, Q., & Running, S. W. (2007). Regional evaporation estimates from flux tower and MODIS satellite data. *Remote Sensing of Environment*, 106, 285-304.
* Cook, B.D., Davis, K. J., Wang, W. G., Desai, A., Berger, B. W., Teclaw, R. M., et al. (2004). Carbon exchange and venting anomalies in an upland deciduous forest in northern Wisconsin, USA. *Agricultural and Forest Meteorology*, 126, 271-295.
* Curtis, P. S., Hanson, P. J., Bolstad, P., Barford, C., Randolph, J. C., Schmid, H. P., et al. (2002). Biometric and eddy-covariance based estimates of annual carbon storage in five eastern North American deciduous forests. *Agricultural and Forest Meteorology*, 113, 3-19.
* Czikowsky, M. J., & Fitzjarrald, D. R. (2009). Detecting rainfall interception in an Amazonian rain forest with eddy flux measurements. *Journal of Hydrology*, 377, 92-105.
* David, T. S., Gash, J. H. C., Valente, F., Pereira, J. S., Ferreira, M. I., & David, J. S. (2006). Rainfall interception by an isolated evergreen oak tree in a Mediterranean savannah. *Hydrological Processes*, 20, 2713-2726.
* Denmead, O. T., & Shaw, R. H. (1962). Availability of soil water to plants as affected by soil moisture content and meteorological conditions. *Agronomy Journal*, 54, 385-390.
* Dewar, R. C. (1995). Interpretation of an empirical model for stomatal conductance in terms of guard cell function. *Plant Cell and Environment*, 18(4), 365-372.
* Diremeyer, P. A., Cash, B., Kinter, J. L., Stan, C., Jung, T., Marx, L., et al. (2012). Evidence for enhanced land-atmosphere feedback in a warming climate. *Journal of Hydrometeorology*, 13, 981-995.
* Dore, S., Kolb, T. E., Montes-Helu, M., Eckert, S. E., Sullivan, B. W., Hungate, B.A., et al. (2010). Carbon and water fluxes from ponderosa pine forests disturbed by wildfire and thinning. *Ecological Applications*, 20(3), 663-683.
* Droogers, P., & Allen, R. G. (2002). Estimating reference evapotranspiration under inaccurate data conditions. *Irrigation and Drainage Systems*, 16, 33-45.
* Entekhabi, D., Asrar, G. R., Betts, A. K., Beven, K. J., Bras, R. L., Duffy, C. J., et al. (1999). An Agenda for land-surface hydrology research and a call for the second International Hydrological Decade. *Bulletin of the American Meteorological Society*, 80(10), 2043-2058.
* Epstein, H. E., Calef, M. P., Walker, M.D., Chapin, F. S., & Starfield, A.M. (2004). Detecting changes in arctic tundra plant communities in response to warming over decadal time scales. *Global Change Biology*, 10(8), 1325-1334.
* Eslamian, S. S., Gohari, S. A., Zareian, M. J., & Firoozfar, A. (2012). Estimating Penman-Monteith reference evapotranspiration using artificial neural networks and genetic algorithm: a case study. *Arabian Journal for Science and Engineering*, 37(4), 935-944.
* Field, R. T., Fritschen, L. J., Kanemasu, E. T., Smith, E. A., Stewart, J. B., Verma, S. B., et al. (1992). Calibration, comparison, and correction of net radiometer instruments used during FIFE. *Journal of Geophysical Research*, 97(D17), 18681-18695.
* Fisher, J. B., Tu, K. P., & Baldocchi, D.D. (2008). Global estimates of the land-atmosphere water flux based on monthly AVHRR and ISLSCP-II data, validated at 16 FLUXNET sites. *Remote Sensing of Environment*, 112(3), 901-919.
* Flanagan, L. B., Wever, L. A., & Carlson, P. J. (2002). Seasonal and interannual variation in carbon dioxide exchange and carbon balance in a northern temperate grassland. *Global Change Biology*, 8, 599-615.
* French, A. N., Schmugge, T. J., Kustas, W. P., Brubaker, K. L., & Prueger, J. (2003). Surface energy fluxes over El Reno, Oklahoma, using high resolution remotely sensed data. *Water Resources Research*, 39(6), 1164.
* Gao, B.-C., & Kaufman, Y. J. (2003). Water vapor retrievals using Moderate Resolution Imaging Spectroradiometer (MODIS) near-infrared channels. *Journal of Geophysical Research*, 108(D13), 4389.
* Granger, R. J., & Gray, D.M. (1989). Evaporation from natural nonsaturated surfaces. *Journal of Hydrology*, 111, 21-29.
* Granier, A., Biron, P., & Lemoine, D. (2000). Water balance, transpiration and canopy conductance in two beech stands. *Agricultural and Forest Meteorology*, 100, 291-308.
* Grünwald, T., & Bernhofer, C. (2007). A decade of carbon, water and energy flux measurements of an old spruce forest at the Anchor Station, Tharandt. *Tellus*, 59B, 387-396.
* Gu, L., Meyers, T., Pallardy, S. G., Hanson, P. J., Yang, B., Heuer, M., et al. (2006). Direct and indirect effects of atmospheric conditions and soil moisture on surface energy partitioning revealed by a prolonged drought at a temperate forest site. *Journal of Geophysical Research*, 111, D16102.
* Guswa, A. J. (2005). Soil-moisture limits on plant uptake: An upscaled relationship for water-limited ecosystems. *Advances in Water Resources*, 28, 543-552.
* Hall, F. G., Huemmrich, K. F., Goetz, S. J., Sellers, P. J., & Nickeson, J. E. (1992). Satellite remote sensing of surface energy balance: success, failures and unresolved issues in FIFE. *Journal of Geophysical Research*, 97, 19061-19089.
* Heinsch, F. A., Heilman, J. L., McInnes, K. J., Cobos, D. R., Zuberer, D. A., & Roelke, D. L. (2004). Carbon dioxide exchange in a high marsh on the Texas Gulf Coast: Effects of freshwater availability. *Agricultural and Forest Meteorology*, 125, 159-172.
* Hilker, T., Gitelson, A., Coops, N. C., Hall, F. G., & Black, A. T. (2011). Tracking plant physiological properties from multi-angular tower-based remote sensing. *Oecologia*, 165(4), 865-876.
* Hulley, G. C., Hughes, T., & Hook, S. J. (2012). Quantifying uncertainties in land surface temperature (LST) and emissivity retrievals from ASTER and MODIS thermal infrared data. *Geophysical Research Letters*, 117, D23113.
* Huntington, J. L., Szilagyi, J., Tyler, S., & Pohll, G. (2011). Evaluating the complementary relationship for estimating evapotranspiration from arid shrublands. *Water Resources Research*, 47, W05533.
* Jarvis, P. G. (1976). The interpretation of leaf water potential and stomatal conductance found in canopies in the field. *Philosophical Transactions of the Royal Society London, Series B*, 273, 593-610.
* Katul, G. G., Manzoni, S., Palmroth, S., & Oren, R. (2010). A stomatal optimization theory to describe the effects of atmospheric CO2 on leaf photosynthesis and transpiration. *Annals of Botany*, 105(3), 431-442.
* Komatsu, H. (2005). Forest categorization according to dry-canopy evaporation rates in the growing season: comparison of the Priestley-Taylor coefficient values from various observation sites. *Hydrological Processes*, 19, 3873-3996.
* Koster, R. D., et al. (2006). GLACE: The global land-atmosphere coupling experiment. Part 1: overview. *Journal of Hydrometeorology*, 7, 590-610.
* Krishnan, P., Meyers, T. P., Scott, R. L., Kennedy, L., & Heuer, M. (2012). Energy exchange and evapotranspiration over two temperate semi-arid grasslands in North America. *Agricultural and Forest Meteorology*, 153, 31-44.
* Kurc, S. A., & Small, E. E. (2004). Dynamics of evapotranspiration in semiarid grassland and shrubland ecosystems during summer monsoon season, central New Mexico. *Water Resources Research*, 40.
* Kustas, W. P., & Anderson, M. C. (2009). Advances in thermal infrared remote sensing for land surface modeling. *Agricultural and Forest Meteorology*, 149, 2071-2081.
* Kustas, W. P., & Norman, J. M. (1999). Evaluation of soil and vegetation heat flux predictions using simple two-source model with radiometric temperature for partial canopy cover. *Agricultural and Forest Meteorology*, 94, 13-29.
* Landsberg, J. J., Kaufman, M. R., Binkeley, D., Isebrands, J., & Jarvis, P. G. (1991). Evaluating progress towards closed forest models based on fluxes of carbon, water and nutrients. *Tree Physiology*, 9, 1-15.
* Lee, T. J., & Pielke, R. (1992). Estimating the soil surface specific humidity. *Journal of Applied Meteorology*, 31, 480-484.
* LeMone, M.A., Tewari, M., Chen, F., Alfieri, J. G., & Niyogi, D. (2008). Evaluation of the Noah Land Surface Model using data from a Fair-Weather IHOP 2002 day with heterogeneous surface fluxes. *Monthly Weather Review*, 136, 4915-4941.
* Leuning, R., (1995). A critical appraisal of a combined stomatal-photosynthesis model for C3 plants. *Plant, Cell and Environment*, 18, 339-355.
* Lhomme, J. P. (1997). A theoretical basis for the Priestley-Taylor coefficient. *Boundary Layer Meteorology*, 82, 179-191.
* Lipson, D. A., Wilson, R. F., & Oechel, W. C. (2005). Effects of elevated atmospheric CO2 on soil microbial biomass, activity, and diversity in a chaparral ecosystem. *Applied and Environmental Microbiology*, 71(12), 8573-8580.
* Liu, S., Lu, L., Mao, D., & Jia, L. (2007). Evaluating parameterizations of aerodynamic resistance to heat transfer using field measurements. *Hydrology and Earth System Sciences*, 11, 769-783.
* Loescher, H. W., Gholz, H. L., Jacobs, J. M., & Oberbauer, S. F. (2005). Energy dynamics and modeled evapotranspiration from a wet tropical forest in Costa Rica. *Journal of Hydrology*, 315, 274-294.
* Malek, E., McCurdy, G., & Giles, B. (1999). Dew contribution to the annual water balances in semi-arid desert valleys. *Journal of Arid Environments*, 42, 71-80.
* Mallick, K., Jarvis, A., Fisher, J. B., Tu, K. P., Boegh, E., & Niyogi, D. (2013). Latent heat flux and canopy conductance based on Penman-Monteith and Bouchet's complementary hypothesis. *Journal of Hydrometeorology*, 14, 419-442.
* Monteith, J. L. (1965). Evaporation and environment. In G. E. Fogg (Ed.), *Symposium of the Society for Experimental Biology*, The State and Movement of Water in Living Organisms, 19. (pp. 205-234). NY: Academic Press, Inc.
* Moran, M. S., Rahman, A. F., Washburne, J. C., Goodrich, D. C., Weltz, M.A., & Kustas, W. P. (1996). Combining the Penman-Monteith equation with measurements of surface temperature and reflectance to estimate evaporation rates of semiarid grassland. *Agricultural and Forest Meteorology*, 80, 87-109.
* Mu, Q., Heinsch, F. A., Zhao, M., & Running, S. W. (2007). Development of a global evapotranspiration algorithm based on MODIS and global meteorology data. *Remote Sensing of Environment*, 111, 519-536.
* Niyogi, D., Kishtawal, C., Tripathi, S., & Govindaraju, R. S. (2010). Observational evidence that agricultural intensification and land use change may be reducing the Indian summer monsoon rainfall. *Water Resources Research*, 46, W05533.
* Noilhan, J., & Planton, S. (1989). A simple parameterization of land surface processes for meteorological models. *Monthly Weather Review*, 117, 536-549.
* Norman, J. M., Kustas, W. P., & Humes, K. S. (1995). A two-source approach for estimating soil and vegetation energy fluxes in observations of directional radiometric surface temperature. *Agricultural and Forest Meteorology*, 77, 263-293.
* Parlange, M. B., & Katul, G. G. (1992). An advection-aridity evaporation model. *Water Resources Research*, 28(1), 127-132.
* Peichl, M., & Arain, M.A. (2006). Above and belowground ecosystem biomass and carbon pools in an age-sequence of white pine plantations in southern Ontario, Canada. *Agricultural and Forest Meteorology*, 140, 51-63.
* Penman, H. L. (1948). Natural evaporation from open water, bare soil, and grass. *Proceedings of Royal Society London*, A193, 120-146.
* Pielke, R. A., Pitman, A., Niyogi, D., Mahmood, R., McAlpine, C., Hossain, F., et al. (2011). Land use/land cover changes and climate: Modeling analysis and observational evidence. *Wiley Interdisciplinary Reviews: Climate Change*.
* Pieruschka, R., Huber, G., & Berry, J. (2010). Control of transpiration by radiation. *Proceedings of the National Academy of Sciences USA*, 107, 13372-13377.
* Priestley, C. H. B., & Taylor, R. J. (1972). On the assessment of surface heat flux and evaporation using large scale parameters. *Monthly Weather Review*, 100, 81-92.
* Ramirez, J. A., Hobbins, M. T., & Brown, T. C. (2005). Observational evidence of the complementary relationship in regional evaporation lends strong support for Bouchet's hypothesis. *Geophysical Research Letters*, 32, L15401.
* Raupach, M. R., & Finnigan, J. J. (1995). Scale issues in boundary-layer meteorology: Surface energy balance in heterogeneous terrain. *Hydrological Processes*, 9, 589-612.
* Rebmann, C., Zeri, M., Lasslop, G., Mund, M., Kolle, O., Schulze, E. D., et al. (2010). Treatment and assessment of the CO2-exchange at a complex forest site in Thuringia, Germany. *Agricultural and Forest Meteorology*, 150, 684-691.
* Roderick, M. L., & Farquhar, G. D. (2004). Changes in Australian pan evaporation from 1970 to 2002. *International Journal of Climatology*, 25, 1077-1090.
* Santos, A. J. B., Quesada, C. A., Da Silva, G. T., Maia, J. F., Miranda, H. S., Miranda, A.C., et al. (2004). High rates of net ecosystem carbon assimilation by Brachiaria pasture in the Brazilian Cerrado. *Global Change Biology*, 10(5), 877-885.
* Scott, R. L., Hamerlynck, E. P., Jenerette, G. D., Moran, M. S., & Barron-Gafford, G. (2010). Carbon dioxide exchange in a semidesert grassland through drought-induced vegetation change. *Journal of Geophysical Research-Biogeosciences*, 115, G03026.
* Scott, R. L., Jenerette, G. D., Potts, D. L., & Huxman, T. E. (2009). Effects of seasonal drought on net carbon dioxide exchange from a woody-plant-encroached semiarid grassland. *Journal of Geophysical Research-Biogeosciences*, 114, G04004.
* Segal, M., Jia, X., Ye, Z., & Pielke, R. A. (1990). On the effect of daytime surface evaporation on pollution dispersion. *Atmospheric Environment*, 24A(7), 1801-1811.
* Shuttleworth, W. J., Gurney, R. J., Hsu, A. Y., & Ormsby, J. P. (1989). FIFE: The variation in energy partition at surface flux sites. In A. Rango (Ed.), *Remote sensing and large scale processes*, Proceedings of the IAHS third international Assembly, Baltimore, MD, May. IAHS Publication, 186, 67-74.
* Slatyer, R. O., & McIlroy, I. C. (1961). *Practical micrometeorology*. Melbourne, Australia: CSIRO (310 pp.).
* Small, E. E., & Kurc, S. A. (2003). Tight coupling between soil moisture and the surface radiation budget in semiarid environments: Implications for land atmosphere interactions. *Water Resources Research*, 39(10).
* Small, E. E., & McConnell, J. R. (2008). Comparison of soil moisture and meteorological controls on pine and spruce transpiration. *Ecohydrology*, 1, 205-214.
* Strasser, U., & Mauser, W. (2001). Modelling the spatial and temporal variations of the water balance for the Weser catchment 1965-1994. *Journal of Hydrology*, 254, 199-214.
* Su, Z. (2002). The surface energy balance system (SEBS) for estimation of turbulent heat fluxes. *Hydrology and Earth System Sciences*, 6(1), 85-99.
* Sugita, M., Usui, J., Tamagawa, I., & Kaihotsu, I. (2001). Complementary relationship with convective boundary layer model to estimate regional evaporation. *Water Resources Research*, 37(2), 353-365.
* Sulkava, M., Luyssaert, S., Zaehle, S., & Papale, D. (2011). Assessing and improving the representativeness of monitoring networks: The European flux tower network example. *Journal of Geophysical Research*, 116.
* Tobin, D. C., Revercomb, H. E., Knuteson, R. O., Lesht, B.M., Strow, L. L., Hannon, S. E., et al. (2006). Atmospheric Radiation Measurement site atmospheric state best estimates for Atmospheric Infrared Sounder temperature and water vapor retrieval validation. *Journal of Geophysical Research*, 111, D09S14.
* Troufleau, D., Lhomme, J. P., Monteny, B., & Vidal, A. (1997). Sensible heat flux and radiometric surface temperature over sparse Sahelian vegetation. I. An experimental analysis of kB-1 parameter. *Journal of Hydrology*, 188-189, 815-838.
* Twine, T. E., Kustas, W. P., Norman, J. M., Cook, D. R., Houser, P. R., Meyers, T. P., et al. (2000). Correcting eddy-covariance flux underestimates over a grassland. *Agricultural and Forest Meteorology*, 103, 279-300.
* Veenendaal, E. M., Kolle, O., & Lloyd, J. (2004). Seasonal variation in energy fluxes and carbon dioxide exchange for a broad-leaved semi-arid savanna (Mopane woodland) in Southern Africa. *Global Change Biology*, 10(3), 318-328.
* Venturini, V., Islam, S., & Rodriguez, L. (2008). Estimation of evaporative fraction and evapotranspiration from MODIS products using a complementary based model. *Remote Sensing of Environment*, 112(1), 132-141.
* Wan, Z., Zhang, Y., Zhang, Q., & Li, Z.-L. (2004). Quality assessment and validation of the MODIS global land surface temperature. *International Journal of Remote Sensing*, 25(1), 261-274.
* Weaver, H. L. (1990). Temperature and humidity flux-variance relations determined by one dimensional eddy correlation. *Boundary Layer Meteorology*, 53, 77-91.
* Williams, C. A., Hanan, N., Scholes, R. J., & Kutsch, W. (2009). Complexity in water and carbon dioxide fluxes following rain pulses in an African savanna. *Oecologia*, 161(3), 469-480.
* Ye, Z., & Pielke, R. (1993). Atmospheric parameterization of evaporation from non-plant-covered surfaces. *Journal of Applied Meteorology*, 32, 1248-1258.
