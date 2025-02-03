# RayleighDistillation.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import colormaps
from .config import *
from .BasicUtility import (
    check_validity,
    sat_specific_humidity,
    mixing_ratio_to_specific_humidity,
)
from .KineticFractionation import (
    kin_frac_factor_sea_evap,
    kin_frac_factor_ice,
)
from .SeaEvaporationIsotopeCalculation import initial_sea_evap_fractionation
from .InitialCondition import get_init_conditions

def adjust_alpha_raindrop_evap(
    alpha_con: float, 
    alpha_evap: float = 1.098, 
    reevap_factor: float = 0.5,
) -> float:
    """
    Adjust fractionation factor considering raindrop evaporation.

    Parameters:
    - alpha_con (float): Equilibrium fractionation factor for condensation(ND).
    - alpha_evap (float): Effective fractionation factor for re-evaporation(ND). Default value is 1.098 (following after Worden et al. 2007).
    - reevap_factor (float): Re-evaporation factor (0–1). Default value is 0.5

    Methods:
    - Adjust the fractionation factor to account for isotopic effects due to
      re-evaporation of raindrops.

    Returns:
    - float: Adjusted fractionation factor.

    References
    - Worden et al. 2007, Eq. 4.14.
    """    
    return alpha_con * (1 - reevap_factor / alpha_evap) / (1 - reevap_factor)

def rayleigh_process(
    temp_init: float, 
    temp_fin: float,
    q_init: float,
    delta_init: float,
    alpha_eff_list: list, 
    alpha_kin_list: list, 
    temp_default_list: list, 
    reevap_factor: float = 0.5,
    BOOL_REEVAP: bool = True, 
    dt: float = 0.5
) -> dict:
    """
    Perform Rayleigh distillation calculations.

    Parameters:
    - temp_init (float): Initial temperature (°C).
    - temp_fin (float): Final temperature (°C).
    - q_init (float): Initial specific humidity (g/kg).
    - delta_init (float): Initial delta value in ratio (ex., X'/X - 1 [ND]).
    - alpha_eff_list (list): List of effective fractionation factors for condensation.
    - temp_default_list (list): List of default temperatures for interpolation.
    - reevap_factor (float): Re-evaporation factor. Default is 0.5.
    - BOOL_REEVAP (bool): Whether to consider re-evaporation. Default is True.
    - dt (float): Temperature step size (°C). Default is 0.5.

    Methods:
    - Validate inputs for initial and final conditions.
    - Create a temperature array for rayleigh process.
    - Perform Rayleigh calculations, including adjustments for re-evaporation.

    Returns:
    - dict: Dictionary containing specific humidity ('q'), delta values ('delta' [‰]), 
      and temperatures ('temp').
        
    Notes:
    - Ensure `temp_init > temp_fin` for condensation.
    - Assume monotonic cooling with condensation.
    
    References:
    - Worden et al. (2007).

    """
    validate_rayleigh_inputs(temp_init, temp_fin, q_init)
    
    temp_rayleigh_list = np.arange(temp_init, temp_fin, -np.abs(dt))
    q_list, delta_list = initialize_rayleigh_arrays(temp_rayleigh_list, q_init, delta_init)
    
    for i, temp in enumerate(temp_rayleigh_list):
        if i == 0:
            continue

        q_list[i] = sat_specific_humidity(temp)
        dq = q_list[i] - q_list[i - 1]

        idx = np.nanargmin(np.abs(temp_default_list - temp))
        alpha = alpha_eff_list[idx]

        if temp >= 0 and BOOL_REEVAP == True:
            # Re-evaporation from raindrop            
            alpha = adjust_alpha_raindrop_evap(alpha, reevap_factor=reevap_factor)
            
        delta_list[i] = rayleigh_step(alpha, q_list[i - 1], dq, delta_list[i - 1])
        

    return {
        "q": q_list,
        "delta": delta_list * 1000, # ratio -> permil
        "temp": temp_rayleigh_list,
    }

def validate_rayleigh_inputs(temp_init: float, temp_fin: float, q_init: float):
    """
    Validate inputs for Rayleigh distilation calculations.

    Paramters:
    - temp_init (float): Initial temperature (°C).
    - temp_fin (float): Final temperature (°C).
    - q_init (float): Initial specific humidity (g/kg).

    Methods:
    - Check that the initial temperature is greater than the final temperature.
    - Ensure the initial specific humidity is positive.
    
    Raises:
    - ValueError: If inputs do not meet the required conditions.
    """

    if temp_init <= temp_fin:
        raise ValueError(
            f"temp_init must be greater than temp_fin. "
            f"Given temp_init: {temp_init}, temp_fin:{temp_fin}"
        )

    if q_init <= 0:
        raise ValueError(f"q_init must be positive.")

def initialize_rayleigh_arrays(
    temp_rayleigh_list: np.ndarray,
    q_init: float, 
    delta_init: float
) -> tuple[list, list]:
    """
    Initialize arrays for Rayleigh distillation.

    Parameters:
    - temp_rayleigh_list (ndarray): Temperature array.
    - q_init (float): Initial specific humidity.
    - delta_init (float): Initial delta value.

    Returns:
    - tuple: Arrays for specific humidity (`q_list`) and delta values (`delta_list`).
    """
    q_list = np.zeros_like(temp_rayleigh_list)
    delta_list = np.zeros_like(temp_rayleigh_list)
    q_list[0] = q_init
    delta_list[0] = delta_init
    return q_list, delta_list

def rayleigh_step(alpha: float, q: float, dq: float, delta: float) -> float:
    """
    Perform a single step of the Rayleigh distillation process.

    Paramters:
    - alpha (float): Fractionation factor (ND).
    - q (float): Specific humidity (g/kg)
    - dq (float): Change in specific humidity (g/kg).
    - delta (float): Previous delta value in ratio (ex., X'/X - 1 [ND]).

    Methods:
    - Update the delta value using the Rayleigh distillation equation.
    
    Returns:
    - float: Updated delta value.

    References
    - Worden et al. 2007, Eq. 4.13.
    """
    return (alpha - 1) / q * dq + delta


def plot_q_dq(
    config: dict,
    rayleigh_dict: dict, 
    surface_dict: dict,
    ISO_TYPE: str = "HDO", 
    title: str = "Rayleigh distillation model",
    cmap: str = "winter",
    num_of_subplot: int = 1,
    ax = None,
    xlim: tuple = (None, None),
    ylim: tuple = (None, None),
):
    """
    Plot Q-dQ relationship for rayleigh distillation.

    Parameters:
    - config (dict): Configuration for the model.
    - rayleigh_dict (dict): Rayleigh distillation results.
    - surface_dict (dict): Surface conditions.
    - ISO_TYPE (str): Isotope type, either "H218O" or "HDO". Default is "HDO".
    - title (str): Plot title. Default is "Rayleigh distillation model".
    - cmap (str): Colormap for plotting. Default is "winter".

    Returns:
    - None: Displays the plot.    
    """

    check_validity(ISO_TYPE, iso_type_list, "ISO_TYPE")

    ylabel = "$\mathsf{\delta D}$" if ISO_TYPE == "HDO" else "$\mathsf{\delta ^{18}O}$"

    colors = colormaps[cmap]
    
    if num_of_subplot == 1:
        fig = plt.figure(layout="tight", figsize=(5, 5))
        fig.suptitle(title)
        gs = GridSpec(1, 1, figure=fig)
        ax = fig.add_subplot(gs[0, 0])
    else:
        if ax == None:
            raise ValueError(
                f"Axes object is not provided while number of figure is more than 1."
            )
    
    for i, temp in enumerate(rayleigh_dict.keys()):
        color = colors(i / len(rayleigh_dict))
        _plot_data(ax, config, rayleigh_dict, surface_dict, temp, color)

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel("Specific humidity [g/kg]")
    ax.set_ylabel(f"{ylabel} [‰]")
    plt.grid()
    
    if num_of_subplot == 1:
        plt.legend()
        plt.show()

def _plot_data(
    ax, 
    config: dict, 
    rayleigh_dict: dict, 
    surface_dict: dict, 
    temp: float, 
    color: tuple
):
    """
    Helper function to plot data on the plot.

    Parameters:
    - ax (matplotlib.axes.Axes): Axes object to plot on.
    - config (dict): Configurations for the model.
    - rayleigh_dict (dict): Rayleigh distillation results.
    - surface_dict (dict): Surface conditions.
    - temp (float): Current temperature.
    - color (tuple): Color for the plot.

    Returns:
    - None: Modifies the given plot object.
    """
    ax.plot(
        rayleigh_dict[temp]["q"],
        rayleigh_dict[temp]["delta"],
        color=color, 
        label=temp, 
        zorder=1
    )
    # ax.scatter(
    #     rayleigh_dict[temp]["q"],
    #     rayleigh_dict[temp]["delta"],
    #     marker="o",
    #     color=color,
    #     zorder=1
    # )
    ax.scatter(
        surface_dict[temp]["q"],
        surface_dict[temp]["delta"],
        marker="^",
        color=color,
        zorder=2
    )
    ax.scatter(
        surface_dict[temp]["snow"],
        surface_dict[temp]["delta_snow"],
        marker="x",
        color=color,
        zorder=3
    )
    ax.scatter(
        config["q_surf"],
        config["delta_q_surf"],
        marker="o",
        color="black",
        zorder=3
    )
    