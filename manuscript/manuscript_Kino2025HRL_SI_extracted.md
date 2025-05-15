---
title: "Description of the Simple Isotopic Model"
author: "Kanon Kino"
date: "2025-05-10"
---

<!--
To convert to pdf: 
    `pandoc manuscript_Kino2025HRL_SI_extracted.md \\n  --template=template.tex --columns=100\\n  --citeproc \\n  --bibliography=references.bib \\n  -o manuscript_Kino2025HRL_SI_extracted.pdf`
-->


# 1. Description of the Simple Isotopic Model

**The content of this document is same with Text S1 of SI for Kino *et al.* (2025). Please cite the paper if you refer to this document.**

---

We developed a new simple isotopic model (SIM) that calculates isotope fractionation during vapor transport. The model focused on vapor transport from low to high latitudes, assuming a monotonic temperature decrease over the Southern Ocean. It consists of three major processes: (i) isotope fractionation during evaporation from the sea surface ([Section 1.1](#1.1-Evaporation-Process-from-Sea-Surface)), (ii) Rayleigh distillation during vapor transport ([Section 1.2](#1.2-Vapor-Transport-Process)), and (iii) generation of snowfall considering its duration and sublimation above the final precipitation site ([Section 1.3](#1.3-Processes-at-the-Final-Site)). The source of this model can be obtained from [IsoFracPy GitHub repository](https://github.com/kanonundgigue/IsoFracPy).

[Section 1.4](#1.4-Calculation-of-Saturation-Vapor-Pressure) describes the calculation of saturated vapor pressure. [Section 1.5](#1.5-Experimental-settings) describes the experimental settings used in this study.

---

## 1.1 Evaporation Process from Sea Surface

The model begins with the evaporation process from the sea surface. Following the formulation of [Craig and Gordon (1965)](#craig1965) and [Yoshimura *et al.* (2008)](#yoshimura2008), the isotopic ratio of evaporated water ($R_{E}$) depends on the difference between the isotopic ratios of surface seawater ($R_{sea}$) and atmospheric water vapor ($R_{a}$), as well as the effective relative humidity ($h_{eff}$), expressed as:

$$
\begin{array}{r}
R_{E} = \alpha_{kin,evap}\frac{\frac{R_{sea}}{\alpha_{eq}} - h_{eff}R_{a}}{1 - h_{eff}} 
\end{array} \quad (1)
$$ <a id="eq1"></a>

where $\alpha_{eq}$ is the equilibrium state factor (described in Section S1.2.2). The effective relative humidity ($h_{eff}$) is defined as:

$$
\begin{array}{r}
h_{eff} = h_{a}\frac{w_{sat,a}}{w_{sat,sea}} \quad (2)
\end{array}
$$ <a id="eq2"></a>

where $h_{a}$ is the relative humidity of the atmosphere above the sea surface, $w_{sat,a}$ is the saturation mixing ratio of the atmosphere, and $w_{sat,sea}$ is the saturation mixing ratio at the sea surface temperature. In case of $R_{a} = R_{E}$ ([Merlivat and Jouzel, 1979](#merlivat1979); [Dar *et al.*, 2020](#dar2020)),

$$
\begin{array}{r}
R_{a} = R_{E} = \frac{\alpha_{k_{mv}}R_{sea}}{\alpha_{eq}\ (1 - h_{eff} + \alpha_{kin,evap}h_{eff})} \quad {(3)}
\end{array}
$$ <a id="eq3"></a>


The kinetic fractionation factor during evaporation from the sea surface ($\alpha_{kin,evap}$) is expressed as a function of wind speed:

$$
\begin{array}{r}
\alpha_{kin,evap} = \left\{ \begin{matrix}
1 - a_{1}\quad & \left( u < 7\ \text{[m/s]} \right) \\
1 - a_{2}u + a_{3}\quad & \left( u \geq 7\ \text{[m/s]} \right)
\end{matrix} \right.\ \quad (4)
\end{array}
$$ <a id="eq4"></a>

where $u$ is the wind speed at the sea surface. The coefficients are shown in [Table I](#table1).

---

<a id="table1"></a>**Table I. Coefficients for kinetic fractionation factor [(Merlivat and Jouzel, 1979)](#merlivat1979)**

| Species       | $a_1$    | $a_2$      | $a_3$      |
|---------------|----------|------------|------------|
| H$_2^{18}$O   | 0.006    | 0.000285   | 0.00082    |
| HDO           | 0.00528  | 0.0002508  | 0.0007216  |


---

## 1.2 Vapor Transport Process

### 1.2.1 Rayleigh Distillation

This model assumes that water vapor cools adiabatically and excess vapor above saturation condenses. Changes in isotope ratios during vapor transport are expressed based on the basic Rayleigh's equation of [Worden *et al.* (2007)](#worden2007). While Worden's model allows part of the precipitation to reevaporate into the system, this model assumes that formed precipitation is immediately removed from the cloud, as in [Jouzel and Merlivat (1984)](#jouzel1984). Because this model assumes a monotonic decrease of air temperature, $q_{n}$ is determined according to saturated vapor pressure ([Equation 21](#eq21)). Under these conditions, the delta value ($\delta$) of vapor is expressed as:

$$
\begin{array}{r}
\delta_{n + 1} = \frac{\alpha - 1}{q_{n}}\left( q_{n + 1} - q_{n} \right) + \delta_{n} \quad (5)
\end{array}
$$ <a id="eq5"></a>

where $n$ is the calculation step and $q$ is the mixing ratio.

The fractionation factor ($\alpha$) takes values from different equations depending on the temperature range. Above 0 °C, we use the equilibrium fractionation factor, and below −20 °C, we use the effective fractionation factor that accounts for supersaturation effects. Between −20°C and 0°C, these two fraction factors are linearly interpolated.

---

### 1.2.2 Equilibrium Fractionation Factor

The equilibrium fractionation factor ($\alpha_{eq}$) for water isotopes is calculated using empirical equations from Majoube ([1971a](#majoube1971a), [1971b](#majoube1971b)).

$$
\begin{array}{r}
\alpha_{eq} = exp\left( \frac{a_{1}}{T^{2}} + \frac{a_{2}}{T} + a_{3} \right) \quad (6)
\end{array}
$$ <a id="eq6"></a>

where $T$ is temperature. The coefficients are shown in [Table II](#table2).

---

<a id="table2"></a>**Table II. Coefficients for equilibrium Fractionation Factor (Majoube, [1971a](#majoube1971a), [1971b](#majoube1971b))**

| Species                     | $a_1$    | $a_2$     | $a_3$       |
|-----------------------------|----------|-----------|-------------|
| H$_2^{18}$O (vapor-liquid)  | 1137     | -0.4156   | -0.002067   |
| HDO (vapor-liquid)          | 24844    | -76.248   | 0.052612    |
| H$_2^{18}$O (vapor-ice)     | 0        | 11.839    | -0.028224   |
| HDO (vapor-ice)             | 16289    | 0         | -0.0945     |

---

### 1.2.3 Isotope Fractionation during Ice Crystal Formation

Below −20 °C, we use the effective fractionation factor ($\alpha_{eff}$) from [Jouzel and Merlivat (1984)](#jouzel1984) during ice crystal formation:

$$
\begin{array}{r}
\alpha_{eff} = \alpha_{kin,ice}\alpha_{eq} \quad (7)
\end{array}
$$ <a id="eq7"></a>

This kinetic fractionation factor is expressed as:

$$
\begin{array}{r}
\alpha_{kin,ice} = \frac{S}{\alpha_{eq}\left( \frac{D}{D^{'}} \right)(S - 1) + 1} \quad (8)
\end{array}
$$ <a id="eq8"></a>

where $D/D'$ is the ratio of molecular diffusivities, as shown in [Table III](#table3).

---

<a id="table3"></a>**Table III. Ratio of molecular diffusivities [(Merlivat, 1978)](#merlivat1978)**

| Species      | $D/D'$    |
|--------------|-----------|
| H$_2^{18}$O  | 1.02849   |
| HDO          | 1.02512   |


---

The supersaturation ratio over ice ($S$) is defined as:

$$
\begin{array}{r}
S = \left\{ \begin{matrix}
1\quad & \left( T \geq - 20\ \text{[°C]} \right) \\
1 - 0.003T\quad & \left( T < - 20\ \text{[°C]} \right)
\end{matrix}\  \right. \quad (9)
\end{array}
$$ <a id="eq9"></a>

Following [Ciais and Jouzel (1994)](#ciais1994), we linearly interpolate the fractionation factor between −20 and 0 °C.

---

## 1.3 Processes at the Final Site

The amount of vapor in clouds that produce snowfall is expressed using the pressure heights of cloud base and top ($p_{btm}$, $p_{top}$):

$$
\begin{array}{r}
Cld = \frac{\left( p_{btm} - p_{top} \right)}{g} \quad (10)
\end{array}
$$ <a id="eq10"></a>

where $g$ is the gravitational acceleration. Given the observed snowfall flux at the final site ($Sn_{obs}\left\lbrack kg/m^{2}/s \right\rbrack$) and its duration ($d$) as external parameters, the initial snowfall flux generated from the cloud ($Sn_{gen}\lbrack kg/kg/s\rbrack$) is expressed as:

$$
\begin{array}{r}
Sn_{gen} = \frac{\frac{Sn_{obs}}{(1 - f)}}{Cld} \quad (11)
\end{array}
$$ <a id="eq11"></a>

The temporal evolution of isotopic ratios is parameterized as follows. We assume that the amount of cloud vapor ($q_{cld}$) remains constant within each time step. From an Eulerian perspective, this assumes continuous vapor advection, with an amount of $Sn_{gen}$ and isotope ratio of $\delta_{n_{final}}$. In order to preserve isotope ratio, we also assume that the isotopic ratio of cloud vapor at timestep $t$ is updated by a weighting mean.

$$
\begin{array}{r}
\delta_{cld,t} = \frac{\delta_{cld}^{*}\left( Cld - Sn_{gen} \right) + \delta_{n_{final}}Sn_{gen}}{Cld} \quad (12)
\end{array}
$$ <a id="eq12"></a>

where $\delta_{cld}^{*}$ is expressed as

$$
\begin{array}{r}
\delta_{cld}^{*} = \frac{\alpha - 1}{Cld}\left( - Sn_{gen} \right) + \delta_{cld,t - 1} \quad (13)
\end{array}
$$ <a id="eq13"></a>

by using [Equation 5](#eq5). The physical meaning of Equations [12](#eq12) and [13](#eq13) requires further investigation in future.

The final isotopic ratio of generated snowfall can be obtained as a weighted average of snowfall at each time step. In practice, since the snowfall flux is constant, it can be calculated as a simple time average:

$$
\begin{array}{r}
\delta_{Sn} = \frac{\sum\left\lbrack \left( \delta_{cld} + 1 \right)\alpha_{eff} - 1 \right\rbrack}{d} \quad (14)
\end{array}
$$ <a id="eq14"></a>

The total amount of snowfall ($Sn_{tot}\lbrack kg/kg\rbrack$) is given by:

$$
\begin{array}{r}
Sn_{tot} = Sn_{gen}d \quad (15)
\end{array}
$$ <a id="eq15"></a>

When a fraction $f$ of the snowfall amount $Sn_{tot}$ sublimates, the amount of sublimation ($q_{sub}$), the updated water vapor amount near the surface ($q_{surf,updated}$), and its isotopic delta value ($\delta_{surf,updated}$) are expressed as:

$$
\begin{array}{r}
q_{sub} = Sn_{tot}f \quad (16)
\end{array}
$$ <a id="eq16"></a>

$$
\begin{array}{r}
q_{surf,updated} = q_{surf} + q_{sub} \quad (17)
\end{array}
$$ <a id="eq17"></a>

$$
\begin{array}{r}
\delta_{surf,updated} = \frac{q_{surf}\left( \delta_{surf} - 1 \right) + q_{sub}\left( \delta_{sub} - 1 \right)}{q_{surf,updated}} + 1 \quad (18)
\end{array}
$$ <a id="eq18"></a>

Note that $\delta_{sub} = \delta_{Sn}$ because no isotope fractionation is assumed to occur during sublimation.

---

## 1.4 Calculation of Saturation Vapor Pressure

The saturation vapor pressure ($e_{s}$) is calculated using the equation from [Sonntag (1990)](#sonntag1990):

$$
\begin{array}{r}
e_{s} = exp\left( \frac{a_{1}}{T} + a_{2} + a_{3}T + a_{4}T^{2} + a_{5}\ln(T) \right) \quad (19)
\end{array}
$$ <a id="eq19"></a>

where $T$ is temperature. The coefficients are shown in [Table IV](#table4).

---

<a id="table4"></a>**Table IV. Coefficients for saturation vapor pressure [(Sonntag, 1990)](#sonntag1990)**

| Phase   | $a_1$       | $a_2$      | $a_3$         | $a_4$         | $a_5$     |
|---------|-------------|------------|---------------|---------------|-------------|
| Liquid  | -6096.9385  | 21.2409642 | -0.02711193   | 0.00001673952 | 2.433502  |
| Ice     | -6024.5282  | 29.32707   | 0.010613868   | -0.0000131988 | -0.49382577 |


---

The conversion from vapor pressure ($e$) to the mixing ratio ($w$):

$$
\begin{array}{r}
w = \epsilon\frac{e}{P} \quad (20)
\end{array}
$$ <a id="eq20"></a>

where $\epsilon$ is molar mass ratio (water/air) and $P$ is the typical surface air pressure (1013 hPa).

The conversion from mixing ratio ($w$) to specific humidity ($q$):

$$
\begin{array}{r}
q = \frac{w}{1 + w} \quad (21)
\end{array}
$$ <a id="eq21"></a>

---

## 1.5 Experimental settings

This study conducted three different series of experiments.

### 1.5.1 Series 1: Sensitivity experiments for Evaporation Process from Sea Surface near the Coast

Series 1 only used the first part of the SIM. [Table V](#table5) summarizes the parameter settings. Sea surface temperature at evaporation source ($T_{sea}$) was set to 0 °C to fit the actual value of the coastal region. Air temperature ($T$) ranged from −20 to 0 °C to assess multiple situations. For Wind speed ($u$) and relative humidity ($h$), we chose values that enabled us to simulate extreme cases to generate a low isotopic ratio of surface air vapor.

---

<a id="table5"></a>**Table V. Parameters used in Series 1**

| Parameter                               | Symbol            | Value                 | Unit |
|-----------------------------------------|-------------------|------------------------|------|
| Sea surface temperature                 | $T_{\mathrm{sea}}$| 0                     | °C   |
| Air temperature                         | $T$               | -20, -15, -10, -5, 0  | °C   |
| Wind speed                              | $u$               | 30                    | m/s  |
| Relative humidity                       | $h$               | 0.6                   | ND   |

---

### 1.5.2 Series 2: Sensitivity experiments for whole parametercombinations (IsoGSM)

Series 2 used the whole part of SIM and tested the influence of all parameters on the near-surface water vapor isotopic ratio at the final site. [Table VI](#table6) summarizes the parameter settings. Sea surface temperature at evaporation source ($T_{sea}$) ranged from 5 to 15 °C according to typical sea surface temperatures in latitudes where storm track activity is enhanced [(Nakamura *et al.*, 2008)](#nakamura2008). $T$ ranged from 0 to 15 °C to be the same as $T_{sea}$ or slightly lower than $T_{sea}$. Possible values were provided for the cloud temperature at the final site ($T_{final}$) and the cloud-top pressure level at the final site ($p_{top}$). Pressure level of the cloud bottom at the final site ($p_{btm}$) was fixed to 700 hPa as only the difference between $p_{top}$ and $p_{btm}$ affects the sensitivity of [Equation 10](#eq10). For $d$ and $Sn_{tot}(1 - f)$, we set values according to Figures 2e and 2g of Kino *et al.* (2025). For $T_{surf}$ and $R_{surf}$, we selected values consistent with neither the southerly nor northerly (the other) cases, as shown in Figures 2a and 2e of Kino *et al.* (2025).

---

<a id="table6"></a>**Table VI. Parameters used in Series 2**
| Parameter                                | Symbol                  | Value                 | Unit      |
|------------------------------------------|--------------------------|------------------------|-----------|
| Sea surface temperature                  | $T_{\mathrm{sea}}$       | 5, 10, 15              | °C        |
| Air temperature                          | $T$                      | 0, 5, 10, 15           | °C        |
| Wind speed                               | $u$                      | 6.5, 30                | m/s       |
| Relative humidity                        | $h$                      | 0.5, 1                 | ND        |
| Kinetic fractionation considered         | —                        | True, False            | ND (bool) |
| Temperature step size                    | $dT$                     | 0.5, 1                 | °C        |
| Upper-air temperature at final site      | $T_{\mathrm{final}}$     | -20, -10               | °C        |
| Cloud-top pressure at final site         | $p_{\mathrm{top}}$       | 400, 500               | hPa       |
| Cloud-bottom pressure at final site      | $p_{\mathrm{bottom}}$    | 700                    | hPa       |
| Sublimation efficiency                   | $f$                      | 0.1, 0.5, 0.9          | ND        |
| Precipitation duration                   | $d$                      | 0.25, 1, 3             | day       |
| Precipitation flux                       | $Sn_{\mathrm{tot}}(1-f)$ | 1, 2                   | mm/day    |
| Near-surface temperature                 | $T_{\mathrm{surf}}$      | 0                      | °C        |
| Default isotope ratio near the surface   | $R_{\mathrm{surf}}$      | -120                   | ‰         |

---

### S1.5.3 Series 3: Sensitivity experiments for whole parameter combinations (MIROC5-iso)

Series 3 was the same as Series 2 except for precipitation flux ($Sn_{tot}(1 - f$)), which was set to 1 and 10 mm/day according to Figure 3g of Kino *et al.* (2025), to investigate the contribution of precipitation intensity on the near-surface water vapor isotopic ratio at the final site.

---

## SUPPLEMENTS REFERENCES

- <a id="ciais1994"></a>Ciais P, Jouzel J. 1994. Deuterium and oxygen 18 in precipitation: isotopic model, including mixed cloud processes. Journal of Geophysical Research Atmospheres 99: 16793–16803. DOI: 10.1029/94JD00412.
- <a id="craig1965"></a>Craig H, Gordon L I. 1965. Deuterium and oxygen 18 variations in the ocean and the marine atmosphere. In Tongiogi E. (ed). Proceedings of Stable Isotopes in Oceanographic Studies and Paleotemperatures. V. Lishi e F., Pisa, Italy; 9--130.
- <a id="dar2020"></a>Dar SS, Ghosh P, Swaraj A, Kumar A. 2020. Craig-Gordon model validation using stable isotope ratios in water vapor over the Southern Ocean. Atmospheric Chemistry and Physics 20: 11435–11449. DOI: 10.5194/acp-20-11435-2020.
- <a id="jouzel1984"></a>Jouzel J, Merlivat L. 1984. Deuterium and oxygen 18 in precipitation: modeling of the isotopic effects during snow formation. Journal of Geophysical Research Atmospheres 89: 11749–11757. DOI: 10.1029/JD089iD07p11749.
- <a id="merlivat1978"></a>Merlivat L. 1978. Molecular diffusivities of H~2~^16^O, HD^16^O, and H~2~^18^O in gases. *The Journal of Chemical Physics* **69:** 2864--2871. DOI: 10.1063/1.436884.
- <a id="merlivat1979"></a>Merlivat L, Jouzel J. 1979. Global climatic interpretation of the deuterium-oxygen 18 relationship for precipitation. Journal of Geophysical Research Oceans 84: 5029–5033. DOI: 10.1029/JC084iC08p05029.
- <a id="majoube1971a"></a>Majoube M. 1971a. Fractionnement en 180 entre la glace et la vapeur d\'eau. *Journal de Chimie Physique et de Physicochimie Biologique* **68:** 625--636. DOI: 10.1051/jcp/1971680625.
- <a id="majoube1971b"></a>Majoube M. 1971b. Fractionnement en oxygène 18 et en deutérium entre l\'eau et sa vapeur. *Journal de Chimie Physique et de Physicochimie Biologique* **68:** 1423--1436. DOI: 10.1051/jcp/1971681423.
- <a id="nakamura2008"></a>Nakamura H, Sampe T, Goto A, Ohfuchi W, Xie SP. 2008. On the importance of midlatitude oceanic frontal zones for the mean state and dominant variability in the tropospheric circulation. *Geophysical Research Letters* **35:** L15709. DOI: 10.1029/2008GL034010.
- <a id="sonntag1990"></a>Sonntag D. 1990. Important new values of the physical constants of 1986, vapour pressure formulations based on the ITS-90, and psychrometer formulae. *Zeitschrift für Meteorologie* **40:** 340--344.
- <a id="worden2007"></a>Worden J, Noone D, Bowman K, Tropospheric Emission Spectrometer Science Team and Data contributors. 2007. Importance of rain evaporation and continental convection in the tropical water cycle. *Nature* **445:** 528--532. DOI: 10.1038/nature05508.
- <a id="yoshimura2008"></a>Yoshimura K, Kanamitsu M, Noone D, Oki T. 2008. Historical isotope simulation using Reanalysis atmospheric data. Journal of Geophysical Research 113(D19): D19108. DOI: 10.1029/2008JD010074.