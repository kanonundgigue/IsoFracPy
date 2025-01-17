# KineticFractionation.py

import numpy as np
import pandas as pd
import warnings
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from .config import *
from .BasicUtility import check_validity
from .EquilibriumFractionation import eq_frac_factor
    
def kin_frac_factor_sea_evap(
    temp_C: float, 
    ISO_TYPE: str = "HDO", 
    u: float  = None
) -> float:
    """
    Calculate the kinetic fractionation factor for sea evaporation.
    
    Parameters:
    - temp_C (float): Temperature (°C).
    - ISO_TYPE (str): Isotope type, either "H218O" or "HDO".
    - u (float): Wind speed [m/s]. If not provided, a default value of 6.5 m/s
      is used for cases where u <7 m/s.

    Methods:
    - Validate the `ISO_TYPE` using the `check_variability` function.
    - Use coefficients for each isotope type based on Merlivat (1978),
    - If `u` < 7 m/s, apply a simplified formula; otherwise, use a linear
      relationship based on `u`.

    Returns:
    - float: Kinetic fractionation factor for sea evaporation.

    Raises:
    - UserWarning: If `u` is not provided, a default value is used with a warning.

    References:
    -  Merlivat and Jouzel (1979): Isotopic fractionation during evaporation.
    """

    check_validity(ISO_TYPE, iso_type_list, "ISO_TYPE")
    
    params_dict = {
        "H218O": [0.006, 0.000285, 0.00082],
        "HDO": [0.00528, 0.0002508, 0.0007216] 
    } # Coefficients from Merlivat (1978)


    if u is None:
        u = 6.5
        warnings.warn(
            f"Wind speed (u) is not provided. A default value (u < 7 [m/s])"
            "is used.", UserWarning
        )
        
    a, b, c = params_dict[ISO_TYPE]

    if u < 7:
        return 1 - a
    else:
        return 1 - (b * u + c)
        
def kin_frac_factor_supersat_ice(
    temp_C: float, 
    supersat_ice_param: float = 0.003, 
    temp_thres: float = -20, 
    ISO_TYPE: str = "HDO", 
    DIFFUSION_REF: str = "M78"
) -> float:
    """
    Calculate the kinetic fractionation factor for supersaturated ice formation.

    Parameters:
    - temp_C (float): Temperature (°C).
    - supersat_ice_param (float): Supersaturation parameter. Default is 0.003.
    - temp_thres (float): Temperature threshold (°C) for supersaturation.
    - ISO_TYPE (str): Isotope type, either "H218O" or "HDO". Default is "HDO".
    - DIFFUSION_REF (str): Diffusion ratio reference, either "M78"
      (Merlivat, 1978) or "C03" (Cappa et al., 2003). Default is "M78".

    Methods:
    - Validate `ISO_TYPE`: and `DIFFUSION_REF` using `check_validity`.
    - Calculate the supersaturation factor (`S`) based on temperature.
    - Calculate the equilibrium fractionation factor (`alpha_eql`) using
      `eq_frac_factor`.
    - Use the diffusion ratio (`DD`) of molecules to calculate the kinetic factor.

    Returns:
    - float: Kinetic fractionation factor for supersaturated ice formation.
    
    References:
    - Jouzel and Merivat (1984): Equations of kinetic fractionation factors.
    """
    
    check_validity(ISO_TYPE, iso_type_list, "ISO_TYPE")
    check_validity(DIFFUSION_REF, diffusion_ref_list, "DIFFUSION_REF")
    
    params_dict = {
        "H218O": {"M78": 1.02849, "C03": 1.03189},
        "HDO": {"M78": 1.02512, "C03": 1.01636}
        }

    DD = params_dict[ISO_TYPE][DIFFUSION_REF]
    
    S = 1 if temp_C >= temp_thres else 1 - supersat_ice_param * temp_C

    alpha_eql = eq_frac_factor(
        temp_C, ISO_TYPE = ISO_TYPE, PHASE_TYPE = "vi"
    )

    alpha_kin = S / (alpha_eql * DD * (S - 1) + 1)

    return alpha_kin

def plot_kin_frac_factor(
    temp_C_list: list,
    supersat_ice_param: float = 0.003, 
    temp_supersat_ice_thres: float = -20, 
    surface_wind: float = 6.5,
):
    """
    Plot kinetic fractionation factors for various isotope types.

    Parameters:
    - temp_C_list (list): List of temperature (°C) for calculating
      fractionation factors.
    - supersat_ice_param (float): Supersaturation parameter for ice. Default is 0.003.
    - temp_supersat_ice_thres (float): Temperature threshold (°C) for
      supersaturation for ice. Default is -20 °C.
    - surface_wind (float): Wind speed (m/s) for calculating sea evaporation
      factors. Default is 6.5 m/s.

    Methods:
    - Split the temperature range into positive and negative subsets.
    - For each isotope type:
        - Caluclate and plot supersaturation factors for negative temperatures
          using `kin_frac_factor_supersat_ice`.
        - Calculate and plot sea evaporation factors for positive temperatures
          using `kin_frac_factor_sea_evap`.

    Returns:
    - None: The function directly generates and displays the plot.

    Notes:
    - The kinetic fractionation factors of supersaturation for ice and sea surface 
      evaporation are displayed in separate regions of the plot.
    """
    # Process temperature ranges
    if temp_C_list[-1] > 0: 
        temp_negative_list = temp_C_list[temp_C_list<0]
        temp_negative_list = np.append(temp_negative_list, 0)
    else:
        temp_negative_list = temp_C_list
    if temp_C_list[0] < -2:
        temp_positive_list = temp_C_list[temp_C_list>=-2]
    else:
        temp_positive_list = temp_C_list


    fig = plt.figure(layout="tight", figsize=(10, 5))
    fig.suptitle("Kinetic fractionation factors")
    gs = GridSpec(1, 2, figure=fig)

    for i, ISO_TYPE in enumerate(iso_type_list):
        ax = fig.add_subplot(gs[0, i])

        subplot_title = f"({'a' if i == 0 else 'b'}) {ISO_TYPE}"
        ax.set_title(subplot_title, loc="left")

        for temp in [0, temp_supersat_ice_thres]:
            ax.axvline(x=temp + temp0_K, 
                   color="black", linewidth=0.7, linestyle="--")

        for j, DIFFUSION_REF in enumerate(diffusion_ref_list):
            alpha_supersat_ice_list = [
                kin_frac_factor_supersat_ice(
                    temp_C,
                    supersat_ice_param,
                    temp_supersat_ice_thres,
                    ISO_TYPE,
                    DIFFUSION_REF
                ) for temp_C in temp_negative_list
            ]
            
            color = "blue" if j == 0 else "red"
            ax.plot(
                [temp_C + temp0_K for temp_C in temp_negative_list], 
                alpha_supersat_ice_list,
                color=color,
                label=DIFFUSION_REF
            )

        alpha_sea_evap_list = [
            kin_frac_factor_sea_evap(
                temp_C, u=surface_wind
            ) for temp_C in temp_positive_list
        ]
        ax.plot(
                [temp_C + temp0_K for temp_C in temp_positive_list], 
                alpha_sea_evap_list,
                color="green", ls="dashed",
                label="MJ79"
            )
        ax.set_xlabel("Temperature [K]")
        ax.set_ylabel("$\\alpha_{kin}$ [ND]")

        plt.legend()
    plt.show()

    