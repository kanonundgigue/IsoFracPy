# InitialConditions.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.gridspec import GridSpec
from .config import *
from .BasicUtility import (
    check_validity,
    delta_from_ratio,
    sat_mixing_ratio,
    mixing_ratio_to_specific_humidity,
)
from .EquilibriumFractionation import eq_frac_factor
from .KineticFractionation import kin_frac_factor_sea_evap
from .SeaEvaporationIsotopeCalculation import (
    initial_sea_evap_fractionation,
    prepare_mixing_ratio
)

def get_init_conditions(
    temp_sea_list: list, 
    temp_air_list:list , 
    ISO_TYPE: str = "HDO",
    pressure: float = 1013.25, 
    SAT_EQ_REF: str = "Sonntag",
    R_sea: float = 1, 
    h_air: float = 1,
) -> dict:
    """
    Calculate initial conditions for Q-dQ analysis based on sea surface evaporation.

    Parameters:
    - temp_sea_list (list): List of sea surface temperatures (°C).
    - temp_air_list (list): List of air temperatures (°C).
    - ISO_TYPE (str): Isotope type, either "H218O" or "HDO". Default is "HDO".
    - pressure (float): Standard atmospheric pressure [hPa]. Default is 1013.25 hPa.
    - SAT_EQ_REF (str): Type of saturated vapor equation, either "Sonntag" or "CC". 
      Default is "Sonntag".
    - R_sea (float): Standard isotopic ratio of sea surface water. Default is 1.
    - h_air (float): Relative humidity of air. Default is 1, assuming saturated
      conditions.

    Methods:
    - Validate the input `ISO_TYPE` and ensure temperature lists are consistent.
    - Calculate saturated mixing ratios for air and sea.
    - compute equilibrium and kinetic fractionation factors for sea evaporation.
    - Derive isotopic ratios and delta values for air after evaporation.
    - Add specific humidity values to the results dictionary.
    
    Returns:
    - results_dict (dict): Dictionary containin calculated variables:
        - `temp_air`: Air temprature list (°C).
        - `temp_sea`: Sea surface temprature list (°C).
        - `w_sat_air`: Saturated mixing ratio of vapor (g/kg).
        - `w_sat_sea`: Saturated mixing ratio of seawater (g/kg).
        - `alpha_eql_sea`: Equilibrium fractionation factors 
           for liquid -> vapor (ND).
        - `alpha_kin_sea`: Kinetic fractionation factors (ND).
        - `R_air`: Isotopic ratio of vapor (ND).
        - `delta_air`: Delta values of vapor (‰).
        - `q_sat_air`: Specific humidity (g/kg).

        Raises:
        - ValueError: If input temperature lists have inconsistent lenghts or
          invalid values (e.g., air temperature exceeding sea temperature).
    """

    check_validity(ISO_TYPE, iso_type_list, "ISO_TYPE")

    temp_list_len = len(temp_air_list)
    if temp_list_len != len(temp_sea_list):
        raise ValueError(f"Lengths of temp_sea_list and temp_air_list are different.")
    if not np.all(temp_air_list <= temp_sea_list):
        raise ValueError(f"Air temperature should be equal or lower than sea temperature.")

    w_sat_air_list, w_sat_sea_list = prepare_mixing_ratio(
        temp_air_list, temp_sea_list
    )

    alpha_eql_sea_list = [
        eq_frac_factor(temp, ISO_TYPE = ISO_TYPE)
        for temp in temp_sea_list
    ]
    alpha_kin_sea_list = [
        kin_frac_factor_sea_evap(temp, ISO_TYPE = ISO_TYPE)
        for temp in temp_sea_list
    ]

    R_air_list = [
        initial_sea_evap_fractionation(
            w_sat_air_list[i], w_sat_sea_list[i], 
            alpha_kin_sea_list[i], alpha_eql_sea_list[i], 
            h_air=h_air, R_sea = R_sea
        ) for i in range(temp_list_len)
    ]

    delta_air_list = [delta_from_ratio(R) for R in R_air_list]

    # Prepare results dictionary
    results_dict = {
        "temp_air": temp_air_list,
        "temp_sea": temp_sea_list,
        "w_sat_air": w_sat_air_list,
        "w_sat_sea": w_sat_sea_list,
        "alpha_eql_sea": alpha_eql_sea_list,
        "alpha_kin_sea": alpha_kin_sea_list,
        "R_air": R_air_list,
        "delta_air": delta_air_list,      
    }

    # Add specific humidity to results dictionary
    results_dict["q_sat_air"] = [
        mixing_ratio_to_specific_humidity(w_sat_air_list[i])
        for i in range(temp_list_len)
    ]

    return results_dict

def plot_init_conditions(
    results_dict: dict, 
    ISO_TYPE: str = "HDO",
    cmap: str = "winter"
):
    """
    Plot initial conditions for Q-dQ analysis.

    Parameters:
    - results_dict (dict): dictionary containing the results of
      `get_init_conditions`.
    - ISO_TYPE (str): Isotope type, either "H218O" or "HDO". Default is "HDO".
    - cmap (str): Colormap for plotting. Default is "winter".

    Methods:
    - Validate the input `ISO_TYPE`.
    - Use the colormap to assing colors to data points.
    - Plot specific humidity (`q_sat_air`) against delta values (`delta_air`)
      for each temperature pair.

    Returns:
    - None: The function directly generates and displays the plot.

    Notes:
    - Each data point is labeled with its corresponding air and sea 
      temperatures.
    """
    
    check_validity(ISO_TYPE, iso_type_list, "ISO_TYPE")
    
    colors = colormaps[cmap]    
    
    fig = plt.figure(layout="tight", figsize=(5, 5))
    fig.suptitle("Initial conditions (determined by sea surface evaporation process)")
    gs = GridSpec(1, 1, figure=fig)

    ax = fig.add_subplot(gs[0, 0])

    for i, temp_air in enumerate(results_dict["temp_air"]):
        color = colors(i / len(results_dict))
        label= (
            "$T_{air}=$" + f"{temp_air} °C, " 
            + "$T_{sea}=$" + f"{results_dict['temp_sea'][i]} °C"
        )
        
        ax.scatter(
            results_dict["q_sat_air"][i], results_dict["delta_air"][i],
            color = color,
            label = label
        )    
    ax.set_xlabel("Specific humidity [g/kg]")
    ax.set_ylabel(f"{iso_label_dict[ISO_TYPE]} [‰]")
    plt.legend()
    plt.show()