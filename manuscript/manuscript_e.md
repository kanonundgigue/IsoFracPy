---
title: Description of a simple isotope model (IsoFracPy module)
date: 2025-02-01
author: Kanon Kino
format:
  pdf:
    toc: false
    toc-depth: 2
    toc-title: Contents
    number-sections: true
    number-depth: 2
    highlight-style: github
---

This document describes a simple isotope model that calculates isotope fractionation during vapor transport. The model focuses on vapor transport from low to high latitudes, assuming monotonical temperature decreases, over the Southern Ocean. It consists of three major processes: (1) isotope fractionation during evaporation from sea surface, (2) Rayleigh distillation during vapor transport, and (3) sublimation of snow at the final site. Source of this model can be obtained by [IsoFracPy GitHub repository](https://github.com/kanonundgigue/IsoFracPy).

# Model Structure

## Evaporation Process from Sea Surface

The model begins with the evaporation process from sea surface, following the formulation of Merlivat and Jouzel (1979). The isotopic ratio of evaporated water ($R_E$) depends on the difference between the isotopic ratios of surface seawater ($R_{sea}$) and atmospheric water vapor ($R_a$), as well as the effective relative humidity ($h_{eff}$), expressed as:

$$
R_E = \alpha_{kin,evap} \frac{R_{sea} / \alpha_{eq} - h_{eff} R_a}{1 - h_{eff}} \tag{1}
$$

where $\alpha_{eq}$ is the equilibrium state factor (described in the next section). The effective relative humidity ($h_{eff}$) is defined as:

$$
h_{eff} = h_a \frac{w_{sat,a}}{w_{sat,sea}} \tag{2}
$$

where $h_a$ is the relative humidity of the atmosphere above the sea surface, $w_{sat,a}$ is the saturation mixing ratio of the atmosphere, and $w_{sat,sea}$ is the saturation mixing ratio at the sea surface temperature.

The kinetic fractionation factor during evaporation from sea surface ($\alpha_{kin,evap}$) is expressed as a function of wind speed:

$$
\alpha_{kin,evap} =
    \begin{cases}
        1 - a_1       \quad& (u <   7 \ \text{[m/s]}) \\
        1 - a_2 u + a_3 \quad& (u \ge 7 \ \text{[m/s]})
    \end{cases} \tag{3}
$$

where $u$ is the wind speed at the sea surface. The coefficients are shown in Table 1.

### Table 1: Coefficients for kinetic fractionation factor (Merlivat, 1978)

|               | $a_1$     | $a_2$       | $a_3$       |
|-----------------------|---------|-----------|-----------|
| $\mathsf{H_2^{18}O}$  | 0.006   | 0.000285  | 0.00082   |
| $\mathsf{HDO}$        | 0.00528 | 0.0002508 | 0.0007216 |

# Vapor Transport Process

## Rayleigh Distillation

This model assumes that water vapor cools adiabatically and excess vapor above saturation condenses. Changes in isotope ratios during vapor transport are expressed based on the basic Rayleigh's equation of Worden et al. (2007). While Worden's model allows part of precipitation to reevaporate into the system, this model assumed that formated precipitation imediately removed from the system, as in Jouzel and Merlivat (1984). Because this model assumes a monotonical decreases of air temperature, $q_n$ is determined according to saturated vapor pressure (eq. 14).
Under these conditions, the delta value ($\delta$) of vapor is expressed as:

$$
    \delta_{n+1} = \frac{\alpha - 1}{q_n} (q_{n+1} - q_n) + \delta_n \tag{4}
$$

where $n$ is the calculation step and $q$ is the mixing ratio. 

The fractionation factor ($\alpha$) takes values from different equations depending on the temperature range. Above 0 °C, we use the equilibrium fractionation factor, and below -20 °C, we use the effective fractionation factor that accounts for iceuration effects. Between -20 and 0 °C, these two fraction factors are linearly interpolated.

## Equilibrium Fractionation Factor

the equilibrium fractionation factor ($\alpha_{eq}$) for water isotopes is calculated using empirical equations from Majoube (1971a; 1971b).

$$
    \alpha_{eq} = \exp \left( \frac{a_1}{T^2} + \frac{a_2}{T} + a_3 \right) \tag{5}
$$


where $T$ is temperature, and subscripts $vl$ and $vi$ represent vapor -> liquid and vapor -> ice phase transitions, respectively. The coefficients are shown in Table 2.

### Table 2： Coefficients for equilibrium Fractionation Factor (Majoube, 1971a; 1971b)
|                                  | $a_1$   | $a_2$     | $a_3$       | 
|-----------------------------------------|---------|-----------|-------------|
| $\mathsf{H_2^{18}O}$ (vapor -> liquid)         | 1137    | -0.4156   | -0.002067   |  
| $\mathsf{HDO}$ (vapor -> liquid)               | 24844   | -76.248   | 0.052612    |  
| $\mathsf{H_2^{18}O}$ (vapor -> ice)         | 0 |11.839  | -0.028224 |             |  
| $\mathsf{HDO}$ (vapor -> ice)               | 16289   |0| -0.0945   |             | 

## Isotope Fractionation during Ice Crystal Formation

Below -20 °C, we use the effective fractionation factor ($\alpha_{eff}$) from Jouzel and Merlivat (1984) during ice crystal formation:

$$
    \alpha_{eff} = \alpha_{kin,ice} \cdot \alpha_{eq} \tag{6}
$$

This kinetic fractionation factor is expressed as：

$$
    \alpha_{kin,ice} = \frac{S}{\alpha_{eq} (D / D') (S - 1) + 1} \tag{7}
$$

where $D/D'$ is the ratio of molecular diffusivities, as shown in Table 3.

### Table 3：Ratio of molecular diffusivities (Merlivat, 1978)

|            | $D/D'$     | 
|-----------------------|---------|
| $\mathsf{H_2^{18}O}$  | 1.02849 |
| $\mathsf{HDO}$        | 1.02512 |


The supersaturation ratio over ice ($S$) is defined as:

$$
    S = 
        \begin{cases}
            1 \quad& (T \ge -20 \ \text{[°C]}) \\
            1 - 0.003 T \quad& (T < -20 \ \text{[°C]})
        \end{cases} \tag{8}
$$

Following Ciais and Jouzel (1994), we linearly interpolate the fractionation factor between -20 and 0 °C.

# Processes at the Final Site

the amount of vapor in clouds that produce snowfall is expressed using the pressure heights of cloud base and top ($p_{btm}$, $p_{top}$):

$$
    Cld = (p_{btm} - p_{top}) / g \tag{9}
$$

where $g$ is the gravitational acceleration. Given the observed snowfall flux at the final site ($Sn_{obs} \mathrm{[kg/m^2/s]}$) and its duration ($d$) as external parameters, the initial snowfall flux generated from the cloud ($Sn_{gen} \mathrm{[kg/kg/s]}$) is expressed as:

$$
    Sn_{gen} = Sn_{obs} / (1 - f) / Cld \tag{10}
$$

The temporal evolution of isotopic ratios is parameterized as follows. We assume that the amount of cloud vapor ($q_{cld}$) remains constant within each time step. From an Eulerian perspective, this assumes continuous vapor advection associated with atmospheric rivers, where the increase in vapor is removed from the cloud as snowfall flux. Changes in the isotopic ratio of vapor during this process are calculated using the temperature-independent Rayleigh distillation equation (eq. 4).

At each time step, the isotopic ratio of cloud vapor is updated by weighting the isotopic ratio of remaining cloud vapor after snow formation ($\delta_{cld}'$) and the isotopic ratio of newly advected vapor ($\delta_{RY}$) according to their respective amounts:

$$
\delta_{cld, updated} = \frac{\delta_{cld}' (q_{cld} - Sn_{gen}) + \delta_{RY} Sn_{gen}}{q_{cld}} \tag{11}
$$

The final isotopic ratio of generated snowfall can be obtained as a weighted average of snowfall at each time step. In practice, since the snowfall flux is constant, it can be calculated as a simple time average:

$$
\delta_{Sn} = \frac{
\sum \left( (\delta_{cld}' + 1) \alpha_{eff} - 1 \right) 
}{d} \tag{12}
$$

The total amount of snowfall ($Sn_{tot} \mathrm{[kg/kg]}$) is given by:

$$
    Sn_{tot} = Sn_{gen} * d \tag{13}
$$

When a fraction $f$ of the snowfall amount $Sn_{tot}$ sublimates, the amount of sublimation ($q_{sub}$), the updated water vapor amount near the surface ($q_{surf, updated}$), and its isotopic delta value ($\delta_{surf, updated}$) are expressed as:

$$
    q_{sub} = Sn_{tot} \cdot f \tag{14}
$$

$$
    q_{surf, updated} = q_{surf} + q_{sub} \tag{15}
$$

$$
    \delta_{surf, updated} = \frac{q_{surf} (\delta_{surf} - 1) + q_{sub} (\delta_{sub} - 1)}{q_{surf, updated}} + 1 \tag{16}
$$


Note that $\delta_{sub} = \delta_{Sn}$ because no isotope fractionation is assumed to occur during sublimation.

# Appendix: Calculation of Saturation Vapor Pressure

The saturation vapor pressure ($e_s$) is calculated using the equation from Sonntag (1990):

$$
    e_s = \exp \left( \frac{a_1}{T} + a_2 + a_3 T + a_4 T^2 + a_5 \ln(T) \right) \tag{17}
$$

where $T$ is temperature. The coefficients are shown in Table 4.

### Table 4： Coefficients for saturation vapor pressure (Sonntag, 1990)
| Phase   | $a_1$        | $a_2$         | $a_3$           | $a_4$           | $a_5$          |
|------|------------|-------------|---------------|---------------|--------------|
| Liquid | -6096.9385 | 21.2409642  | -2.711193e-2  | 1.673952e-5   | 2.433502     |
| Ice | -6024.5282 | 29.32707    | 1.0613868e-2  | -1.3198825e-5 | -0.49382577  |

The conversion from vapor pressure ($e$) to mixing ratio ($w$):

$$
    w = \epsilon \frac{e}{P} \tag{18}
$$

where $\epsilon$ is molar mass ratio (water / air) and $P$ is the typical surface air pressure (1013 hPa).

The conversion from mixing ratio ($w$) to specific humidity ($q$):

$$
    q = \frac{w}{1 + w} \tag{19}
$$

# Calculation Procedure

## Initial Condition Setup

The initial condition is determined by the evaporation process from sea surface. Specifically, the initial isotope ratio of vapor is determined by the following parameters: sea surface temperature ($T_{sea}$) and its isotope ratio ($R_{sea}$), and atmospheric conditions above the sea surface including temperature ($T_{a}$), relative humidity ($h_a$), and wind speed ($u$).

## Changes in Isotope Ratio during Vapor Transport

Changes in isotope ratios during vapor transport are calculated as a Rayleigh distillation process. Following Jouzel and Merivat (1984), we assume that temperature monotonically decreases along the transport path. We also assume that precipitation formed by condensation is immediately removed from the system and is not re-incorpolated into the vapor mass.

## Water Vapor Isotope Ratio at the Final Site

At the final site, a portion of snowfall formed in the upper air (its tempearture is assumed below 0 °C) sublimates in the near-surface layer (assumed to be unsaturated condition). This sublimation process alters the isotope ratio of vapor near the surface. Therefore, the final isotope ratio of vapor near surface is influced by this snow sublimation process. This  can be determined by the following parameters: near-surface conditions before snowfall occurs (relative humidity: $h_{surf}$, specific humidity: $q_{surf}$, isotope ratio of vapor: $R_{surf}$), sublimation rate of snow ($f$), and the virtual snowfall amount (in terms of relative humidity increment) ($dh$).

# Sensitivity Experiments

To understand the model behavior, we conducted comprehensive sensitivity experiments for all parameters (Tables 5 and 6). Based on these results, we identified parameters with low sensitivity to the final isotope ratio of vapor near surface ($R_{surf, updated}$) (Table 5) and assigned standard values to them. We then conducted additional sensitivity expereiments for parameters with high sensitivity (Table 6) to quantitatively evaluate their impacts on $R_{surf, updated}$.

### Table 5: Parameters with low sensitivity to the final isotope ratio of vapor near surface ($R_{surf, updated}$). The assigned standard values are indicated as bold.

| Parameter                          | Symbol       | Value                     | Unit | Note |
|------------------------------|------------|------------------------|------|------|
| Sea surface temperature at evaporation source           | $T_{sea}$  | **[5, 10, 15]**        | °C   | *1   |
| Air temperature at evaporation source              | $T_a$      | [0, 5, 10], **[5, 10, 15]** | °C   |      |
| Wind speed at evaporation source               | $u$        | **6.5**, 35            | m/s  | *2     |
| Realative humidity at evaporation source           | $h$        | 0.5, **1**             | ND   |      |
| Consideration of kinetic fractionation during ice formation         | -          | **True**, False        | bool |      |
| Near-surface temperature at the final site    | $T_{surf}$ | **0**                  | °C   | *3     |
| Near-surface relative humidity at the final site | $h_{surf}$ | 0.7, 0.8               | ND   | *4   |


**Notes:**
1. (*1) Typical sea surface temperatures in latitudes where storm track activity is enhanced range from 5 to 15 °C (Nakamura et al. 2008).
2. (*2) 6.5 m/s was assigned to make $\alpha_{kin,evap}$ constant while realistic values would be higher than that.
3. (*3) Both observation and GCM experiments show that the Antarctic coastal air temperature is approximately 0 °C.
4. (*4) Based on observations and GCM experiments, we assigned 0.75 as a typical value.

   
### Table 6: Parameters with high sensitivity to the final isotope ratio of vapor near surface ($R_{surf, updated}$). 

| Parameter                          | Symbol       | Value                     | Unit | Note |
|----------------------------------|------------|----------------|------|------|
| Upper-air temperature at the final site            | $T_{final}$ | [-30, -10]     | °C   |      |
| Default isotope ratio of vapor near the surface | $R_{surf}$  | [-120, -100]   | ‰    |      |
| Virtual snowfall amount (in terms of crelative humidity increment) ($dh$).    | $h_{surf}$  | [0.1, 0.2]     | ND   |      |
| Sublimation efficiency of snowfall                 | $f$        | [0.1, 0.5, 0.9] | ND   |      |