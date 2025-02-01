# EquilibriumFractionation.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from .config import *
from .BasicUtility import check_validity

def combine_alpha_eq(
    temp_list: list,
    alpha_eqi_list: list,
    alpha_eql_list: list,
    temp_thres_min: float = -20,
    temp_thres_max: float = 0,
) -> list:
    """
    Combine alpha_eqi and alpha_eql values based on temperature conditions 
    with linear interpolation.

    Parameters:
    - temp_list (list): A list of temperature values (°C).
    - alpha_eqi_list (list): A list of fractionation factors for "vi" (vapor -> ice).
    - alpha_eql_list (list): A list of fractionation factors for "vl" (vapor -> liquid).
    - temp_thres_min (float): Minimum temperature threshold for the "vi" phase.
      Default is -20 °C.
    - temp_thres_max (float): Maximum temperature threshold for the "vi" phase.
      Default is 0 °C.

    Methods:
    - Assign alpha_eqi values for temperatures <= temp_thres_min.
    - Assign alpha_eql values for temperatures <= temp_thres_max.
    - Use linear interpolation for intermediate temperature values.

    Returns:
    - list: A list of combined equilibrium fractionation factors.

    Notes:
    - Input lists must have the same length.
    - Linear interpolation is performed using Pandas' `interpolate` method.
    """

    # Convert inputl list to NumPy arrays for efficient processing
    temp_array = np.array(temp_list)
    alpha_eqi_array = np.array(alpha_eqi_list)
    alpha_eql_array = np.array(alpha_eql_list)

    # Initialize alpha_eq_array with NaN values
    alpha_eq_array = np.full_like(temp_array, np.nan, dtype=np.float64)

    # Assign values based on temperature thresholds
    alpha_eq_array[temp_array <= temp_thres_min] = alpha_eqi_array[
        temp_array <= temp_thres_min
    ]
    alpha_eq_array[temp_array > temp_thres_max] = alpha_eql_array[
        temp_array > temp_thres_max
    ]

    # Fill NaN values using linear interpolation
    alpha_eq_array = pd.Series(alpha_eq_array).interpolate(method="linear").to_numpy()

    # Return the interpolated values as a list
    return alpha_eq_array.tolist()
    
def eq_frac_factor(
    temp_C: float,                                             
    ISO_TYPE: str = "HDO", 
    PHASE_TYPE: str = "vl"
) -> float:
    """
    Calculate equilibrium fractionation factor for a given isotope type.

    Parameters:
    - temp_C (float): Temperature (°C).
    - ISO_TYPE (str): Isotope type, either "H218O" or "HDO". Default is "HDO".
    - PHASE_TYPE (str): Phase type, either "vl" (vapor -> liquid) or "vi" (vapor -> ice). 
      Default is "vl".

    Methods:
    - Validate `ISO_TYPE` and `PHASE_TYPE` using the `check_validity` function.
    - Convert the input temperature to Kelvin.
    - Calculate the equilibrium fractionation factor for specified `ISO_TYPE` and `PHASE_TYPE`.

    Returns:
    - float: Equilibrium fracionation factor (ND).

    Raises:
    - ValueError: If `PHASE_TYPE` is not "vl" or "vi".
    
    References:
    - This function uses the empirical equations of Majoube (1971a, 1971b).
    """

    check_validity(ISO_TYPE, iso_type_list, "ISO_TYPE")
    check_validity(PHASE_TYPE, phase_type_list, "PHASE_TYPE")

    temp = temp_C + temp0_K
    
    params_dict = {
        "H218O": {
            "vl": [1137, -0.4156, -0.002067],
            "vi": [0, 11.839, -0.028224],
        },
        "HDO": {
            "vl": [24844, -76.248, 0.052612],        
            "vi": [16289, 0, -0.0945]
        },
    }    
   
    a1, a2, a3 = params_dict[ISO_TYPE][PHASE_TYPE]        
    alpha =  np.exp(a1 / temp**2 + a2 / temp + a3)
    

    return alpha

def prepare_combined_alpha_eq(
    temp_list: list, 
    temp_thres_min: float = -20,
    temp_thres_max: float = 0,    
    ISO_TYPE: str = "HDO"
) -> list:
    """
    Prepare a combined list of equilibrium fractionation factors.

    Parameters:
    - temp_list (list): A list of temperature values (°C).
    - temp_thres_min (float): Minimum temperature threshold for the "vi" phase.
      Default is -20 °C.
    - temp_thres_max (float): Maximum temperature threshold for the "vl" phase.
      Default is 0 °C.
    - ISO_TYPE (str): Isotope type, either "H218O" or "HDO". Default is "HDO".

    Methods:
    - Validate `ISO_TYPE` using the `check_validity` function.
    - Calculate `alpha_eqi` (vapor -> ice) and `alpha_eql` (vapor -> liquid)
      for the given temperatures.
    - Combine the calculated values using linear interpolation between 
      `temp_thres_min` and `temp_thres_max`.

    Returns:
    - list: A list of combined equilibrium fractionation factors.

    Notes:
    - This function calls `eq_frac_factor` and `combine_alpha_eq`.
    """
    check_validity(ISO_TYPE, iso_type_list, "ISO_TYPE")
    
    alpha_eqi_list = [
        eq_frac_factor(temp, ISO_TYPE=ISO_TYPE, PHASE_TYPE="vi")
        for temp in temp_list
    ]
    alpha_eql_list = [
        eq_frac_factor(temp, ISO_TYPE=ISO_TYPE, PHASE_TYPE="vl")
        for temp in temp_list
    ]
    
    alpha_eq_list = combine_alpha_eq(
        temp_list, alpha_eqi_list, alpha_eql_list, 
        temp_thres_min=temp_thres_min, temp_thres_max=temp_thres_max,
    )

    return alpha_eq_list
    
def plot_eq_frac_factor(
    temp_C_list: list,
    temp_thres_min: float = -20,
    temp_thres_max: float = 0,    
):
    """
    Plot equilibrium fractionation factors for different isotope and phase types.

    Parameters:
    - temp_C_list (list): List of temperatures (°C) for which the 
      equilibrium fractionation factors are calculated.
    - temp_thres_min (float): Minimum temperature threshold for the "vi" phase.
      Default is -20 °C.
    - temp_thres_max (float): Maximum temperature threshold for the "vi" phase.
      Default is 0 °C.

    Methods:
    - For each isotope type, calculate combined equilibrium fractionation factors
      using `prepare_combined_alpha_eq`.
    - Plot results with temperature (converted to Kelvin).

    Returns:
    - None: The function directly generates and displays the plot.
    """
    fig = plt.figure(layout="tight", figsize=(10, 5))
    fig.suptitle("Equilibrium fractionation factors")
    gs = GridSpec(1, 2, figure=fig)

    for i, ISO_TYPE in enumerate(iso_type_list):
        ax = fig.add_subplot(gs[0, i])

        subplot_title = f"({'a' if i == 0 else 'b'}) {ISO_TYPE}"
        ax.set_title(subplot_title, loc="left")

        for temp in temp_thres_min, temp_thres_max:
            ax.axvline(
                x=temp + temp0_K, color="black", linewidth=0.7, linestyle="--"
            )
            
        alpha_list = prepare_combined_alpha_eq(
            temp_C_list,  
            temp_thres_min=temp_thres_min, 
            temp_thres_max=temp_thres_max,
            ISO_TYPE=ISO_TYPE
        )

        ax.plot(
            [temp_C + temp0_K for temp_C in temp_C_list], 
            alpha_list, 
        )
        ax.set_xlabel("Temperature [K]")
        ax.set_ylabel("$\\alpha_{eq}$ [ND]")
    
    plt.show()
    

            
            
        