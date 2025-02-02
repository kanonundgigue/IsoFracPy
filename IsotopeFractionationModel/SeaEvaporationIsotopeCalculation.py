# SeaEvaporationIsotopeCalculation.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from .config import *
from .BasicUtility import (
    check_validity,
    delta_from_ratio,
    sat_mixing_ratio,
    mixing_ratio_to_specific_humidity,
)

def effective_humidity(
    w_sat_air: float, 
    w_sat_sea: float, 
    h_air: float = 1
) -> float:
    """
    Calculate effective relative humidity.

    Parameters:
    - w_sat_air (float): Saturated vapor mixing ratio in the air (g/kg).
    - w_sat_sea (float): Saturated vapor mixing ratio at the sea surface (g/kg).
    - h_air (float): Relative humidity of air. Default is 1, assuming　saturated 
      conditions.

    Methods:
    - Calculate the ratio of saturated mixing ratios for air and sea, adjusted
      by the relative humidity of air.
    
    Returns:
    - float: Effective relative humidity (ND).
    """

    return h_air * w_sat_air / w_sat_sea

def prepare_mixing_ratio(
    temp_air_list: list, 
    temp_sea_list: list
)-> tuple:
    """
    Prepare saturated mixing ratio lists for air and sea temperatures.

    Parameters:
    - temp_air_list (list): List of air temperatures (°C).
    - temp_sea_list (list): List of sea surface temperatures (°C).

    Methods:
    - Use `sat_mixing_ratio` to calculate saturated mixing ratios for given
      temperature lists.

    Returns:
    - tuple: Two lists containing saturated mixing ratios for air and sea
      (both in g/kg).
    """
    def _get_w_list(temp_list):
        return [sat_mixing_ratio(temp) for temp in temp_list]

    w_sat_air_list = _get_w_list(temp_air_list)
    w_sat_sea_list = _get_w_list(temp_sea_list)

    return w_sat_air_list, w_sat_sea_list

def plot_effective_humidity(
    temp_air_list: list, 
    temp_sea_list: list,
    h_air: float = 1
):
    """
    Plot effective relative humidity and saturated mixing ratios.

    Parameters:
    - temp_air_list (list): List of air temperatures (°C).
    - temp_sea_list (list): List of sea surface temperatures (°C).
    - h_air (float): Relative humidity of air. Default is 1.

    Methods:
    - Calculate saturated mixing ratios for air and sea using
      `prepare_mixing_ratio`.
    - Compute effective relative humidity for all temperature pairs.
    - Plot saturated mixing ratios and effective humidity on separate subplots.

    Returns:
    - None: The function directly generates and displays the plots.
    """

    w_sat_air_list, w_sat_sea_list = prepare_mixing_ratio(
        temp_air_list, temp_sea_list
    ) 

    h_eff_list = [
        100 * effective_humidity(w_air, w_sea, h_air)
        for w_air, w_sea in zip(w_sat_air_list, w_sat_sea_list)
    ]

    plot_data = {
        "air": {
            "temp": temp_air_list, "w": w_sat_air_list, "color": "red", "marker": "o"
        },
        "sea": {
            "temp": temp_sea_list, "w": w_sat_sea_list, "color": "blue", "marker": "x"
        },
    }

    fig = plt.figure(layout="tight", figsize=(10, 5))
    fig.suptitle("Effective humidity")
    gs = GridSpec(1, 2, figure=fig)

    ax1 = fig.add_subplot(gs[0, 0])
    _plot_saturated_mixing_ratio(ax1, plot_data)

    ax2 = fig.add_subplot(gs[0, 1])
    _plot_effective_humidity(ax2, temp_air_list, h_eff_list)

    plt.show()
    
def _plot_saturated_mixing_ratio(ax, plot_data: dict):
    """
    Helper function to plot saturated mixing ratios.

    Parameters:
    - ax (matplotlib.axes.Axes): Axes object to plot on.
    - plot_data (dict): Data for plotting, including temperature, mixing ratio,
      color, and marker for air and sea.

    Returns:
    - None: The function modifies the given Axes object.
    """
    ax.set_title("(a) Saturated mixing ratio", loc="left")
    ax.axvline(x=0, color="black", linewidth=0.7, linestyle="--")
    for key, data in plot_data.items():
        ax.scatter(
            data["temp"], data["w"],
            color=data["color"], marker=data["marker"], alpha=0.5,
            label=key
        )
    ax.set_xlabel("Temperature [°C]")
    ax.set_ylabel("Mixing ratio [g/kg]")
    ax.legend()

def _plot_effective_humidity(ax, temp_air_list: list, h_eff_list: list):
    """
    Helper function to plot effective humidity.

    Parameters:
    - ax (matplotlib.axes.Axes): Axes object to plot on.
    - temp_air_list (list): List of air temperature (°C).
    - h_eff_list (list): List of effective humidity values (%).

    Returns:
    - None: The function modifies the given Axes object.
    """
    ax.set_title("(b) Effective humidity", loc="left")
    ax.axvline(x=0, color="black", linewidth=0.7, linestyle="--")
    ax.axhline(y=100, color="black", linewidth=0.7, linestyle="--")
    ax.scatter(temp_air_list, h_eff_list, color="green")        
    ax.set_xlabel("Temperature [°C]")
    ax.set_ylabel("Humidity [%]")
    
def initial_sea_evap_fractionation(
    w_sat_air: float, 
    w_sat_sea: float,    
    alpha_kin: float, 
    alpha_eql: float, 
    h_air: float = 1,
    R_sea: float = 1,
) -> float:
    """
    Calculate the initial isotope ratio of evaporated air.

    Parameters:
    - w_sat_air (float): Saturated vapor mixing ratio in the air (g/kg).
    - w_sat_sea (float): Saturated vapor mixing ratio at the sea surface (g/kg).
    - alpha_kin (float): Kinetic fractionation factor.
    - alpha_eql (float): Equilibrium fractionation factor.
    - h_air (float): Relative humidity of air. Default is 1, assuming saturated 
      conditions.    
    - R_sea (float): Standard isotopic ratio of sea surface water. Default is 1.
    
    Returns:
    - float: Isotope ratio of air.

    References:
    - The main equation was modified after Merlivat and Jouzel (1979).
    """
    h_eff = effective_humidity(w_sat_air, w_sat_sea, h_air)

    return alpha_kin * R_sea / alpha_eql / (1 - h_eff + h_eff * alpha_kin)
