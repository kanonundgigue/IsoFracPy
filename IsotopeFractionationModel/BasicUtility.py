# BasicUtility.py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from .config import *
    
def check_validity(value, valid_values: list, value_name: str):
    """
    Check if a given value is within a list of valid values.

    Parameters:
    - value: The value to be checked.
    - valid_values (list): A list of valid values that the value should match.
    - value_name (str): A string representing the name of the value being checked,
                        used for a generating error messages.
    Methods:
    - Check if the value is in the list of valid_values. 
    
    Raises:
    - ValueError: If the value is not in the list of valid_values.
    """
    
    if value not in valid_values:
        raise ValueError(f"Invalid {value_name}: {value}. Must be one of {valid_values}")

def delta_from_ratio(R: float) -> float:
    """
    Convert isotope ratio to delta value.

    Parameters:
    - R (float): Isotope ratio.

    Methods:
    - Calculate delta value using the formula: (R - 1) * 1000.

    Returns: 
    - float: Isotope delta value [‰].
    """

    return (R - 1) * 1000
                  
def sat_specific_humidity(temp_C: float, SAT_EQ_REF: str = "Sonntag") -> float:
    """
    Calculate saturated specific humidity from temperature.

    Parameters:
    - temp_C (float): Temperature (°C).
    - SAT_EQ_REF (str): Type of saturated vapor equation, either "Sonntag" or "CC". Default is "Sonntag".

    Methods:
    - Calculate saturated mixing ratio.
    - Convert it to specific humidity.

    Returns:
    - float: Saturated specific humidity (g/kg).
    """
    w = sat_mixing_ratio(temp_C, SAT_EQ_REF = SAT_EQ_REF)
    return mixing_ratio_to_specific_humidity(w)
    
def sat_mixing_ratio(
    temp_C: float, 
    pressure: float = 1013.25,
    SAT_EQ_REF: str = "Sonntag"
) -> float:
    """
    Calculate saturated vapor mixing ratio from temperature.

    Parameters:
    - temp_C (float): Temperature (°C).
    - pressure (float): Standard atmospheric pressure [hPa]. Default is 1013.25 hPa.
    - SAT_EQ_REF (str): Type of saturated vapor equation, either "Sonntag" or "CC". Default is "Sonntag".

    Methods:
    - Use the specified saturation equation to calculate saturated vapor pressure.
    - Convert the saturated vapor pressure to mixing ratio using the molar mass ratio.
    
    Returns:
    - float: Saturated vapor mixing ratio (g/kg).
    """ 
    temp = temp_C + temp0_K
    method = SaturationFunc.get_method(SAT_EQ_REF)
    es = method(temp)
    epsilon = Mw / Ma # Molar mass ratio
    ws = epsilon * es / pressure

    return 1000 * ws

def mixing_ratio_to_specific_humidity(w1000: float) -> float:
    """
    Convert mixing ratio to specific humidity.

    Parameters:
    - w1000 (float): Mixing ratio (g/kg).

    Methods:
    - Convert mixing ratio to specific humidity using the formula: q = w (1 + w).
    
    Returns
    - float: Specific humidity (g/kg).
    """""
    
    w = w1000 / 1000
    q = w / (1 + w)
    return q * 1000

class SaturationFunc:
    """
    A class for handling saturation equations and methods.

    Attributes:
    - method_map (dict): Mapping of aliases to saturation equation methods.

    Methods:
    - clausius_clapeyron: Calculate saturated vapor pressure using Clausius-Clapeyron equation.
    - Sonntag: Calculate saturated vapor pressure using Sonntag equation.
    - get_method: Retrieve the appropriate method based on alias.
    """

    method_map ={
        "cc": "clausius_clapeyron",
        "clausius_clapeyron": "clausius_clapeyron",
        "clausius clapeyron": "clausius_clapeyron",
        "sonntag":  "sonntag",
        "sonntag1990": "sonntag",
        "sonntag (1990)": "sonntag",
    }
    
    @staticmethod
    def clausius_clapeyron(temp: float) -> float:
        """
        Calculate saturated vapor pressure using Clausius-Clapeyron equation.
        
        Parameters:
        - temp (float): Temperature (K).

        Methods:
        - Calculate saturated vapor pressure using the formula with constants for Clausius-Clapeyron.
    
        Returns:
        - float: Saturated vapor pressure (hPa).
        """
        return e0 * np.exp(Lv / (R_const / (Mw / 1000)) * (1 / temp0_K - 1 /temp))    
        
    @staticmethod
    def sonntag(temp: float) -> float:
        """
        Calculate saturated vapor pressure.
        
        Parameters:
        - temp (float): Temperature (K).

        Methods:
        - Determine whther to use the liquid or ice parameter set based on the temperature.
        - Calculate saturated vapor pressure using Sonntag equation.
        
        Returns:
        - float: Saturated vapor pressure (hPa).

        Notes:
        - Switches between liquid and ice formulations at 0 °C.
        
        References:
        - Sonntag (1990): Sonntag, D, “Important new values of the physical constants of 1986, 
          vapour pressure formulations based on the ITS-90 and psychrometer formulae”, 
          Zeitschrift fur Meteorologie, 1990, 40(5), 340-344.
        """
    
        params ={
            "liquid": [-6096.9385, 21.2409642, -2.711193e-2, 1.673952e-5, 2.433502],
            "ice": [-6024.5282, 29.32707, 1.0613868e-2, -1.3198825e-5, -0.49382577],
        }
            
        key = "liquid" if temp >= temp0_K else "ice"
        a1, a2, a3, a4, a5 = params[key]
    
        return np.exp(a1/temp + a2 + a3 * temp + a4 * temp**2  + a5 * np.log(temp)) / 100  

    @staticmethod
    def default(temp):
        """
        Raise an error for unsupported saturation equation references.

        Parameters:
        - temp: This parameter is not used in the function 
                but is required for consistency with other methods.

        Methods:
        - Raise a ValueError indicating that the specified saturation equation references is unsupported.

        Raises:
        - ValueError: If an unsupported SAT_EQ_REF value is encountered.
        """
        raise ValueError(f"Unsupported SAT_EQ_REF value.")
        
    @classmethod
    def get_method(cls, SAT_EQ_REF: str):
        """
        Retrieve the appropriate saturation equation method based on the alias.

        Parameters:
        - SAT_EQ_REF (str): Saturation equation alias.

        Methods:
        - Map the alias to the corresponding method in the class.

        Returns:
        - Callable: The saturation equation method.
        """
        method_name = cls.method_map.get(SAT_EQ_REF.lower(), "default")
        return getattr(cls, method_name)
        

def plot_sat_vapor_pressure(temp_C_list):
    """
    Plot saturated vapor pressures against temperature.

    Parameters:
    - temp_C_list (list): List of temperatures (°C).

    Methods:
    - Calculate saturated vapor pressures for each temperature using different methods.
    - Plot the results with a logarithmic scale for the pressure axis.
    
    Returns:
    - None: The function directly generates and displays the plot.    
    """
    fig = plt.figure(layout="tight", figsize=(5, 5))
    fig.suptitle("Saturated vapor pressures")
    gs = GridSpec(1, 1, figure = fig)
    ax = fig.add_subplot(gs[0, 0])
    
    ax.axhline(
        y=1, color="black", linewidth=0.7, linestyle="--"
    )
    ax.axvline(
        x=temp0_K, color="black", linewidth=0.7, linestyle="--"
    )
    
    for SAT_EQ_REF in sat_eq_ref_list:
        es_list = []
        for temp_C in temp_C_list:
            temp = temp_C + temp0_K
            method = SaturationFunc.get_method(SAT_EQ_REF)
            es_list.append(method(temp))
    
        ax.plot(
            [temp_C + temp0_K for temp_C in temp_C_list], 
            es_list, label=SAT_EQ_REF
        )
        
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("Pressure [hPa]")
    ax.set_yscale("log")
    plt.legend()
    plt.show()