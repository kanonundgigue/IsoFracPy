# PrepareIsotopeFractionationFactors.py
import numpy as np
import pandas as pd
from .config import *

from .BasicUtility import check_validity
from .EquilibriumFractionation import (
    combine_alpha_eq, 
    eq_frac_factor,
    prepare_combined_alpha_eq,
)
from .KineticFractionation import (
    kin_frac_factor_sea_evap,
    kin_frac_factor_ice,
)

def interp_alpha_kin_ice(alpha_kin_list, invalid_thres: float = 1) -> list: 
    """
    Interpolate alpha_kin values with linear interpolation.

    Parameters:
    - alpha_kin_list (list): A list of alpha_kin values, where values >= 1 are need to be replaced.

    Method:
    - Convert the input list to a NumPy array.
    - Replace values >= 1 with NaN to mark them as invalid.
    - Ensure the last value in the array is set to 1.
    - Use pandas' linear interpolation to fill NaN values based on surrounding values.
    - Convert the interpolated array back to a Python list and return it.

    Returns:
    - list: A list of interpolated alpha_kin values.

    Reference:
    - The interpolation follows the approach described by Ciais and Jouzel (1994) and Yoshimura (2009).
    """
    alpha_kin_array = np.array(alpha_kin_list)
    alpha_kin_array[alpha_kin_array >= invalid_thres] = np.nan
    alpha_kin_array[-1] = invalid_thres  # Assuming the ascending order of temperature
    alpha_kin_array = pd.Series(alpha_kin_array).interpolate(method="linear").to_numpy()
    return alpha_kin_array.tolist()

def prepare_frac_factors(
    temp_list: list, 
    temp_thres_min: float = -20,
    temp_thres_max: float = 0,        
    ISO_TYPE: str = "HDO"
) -> dict:
    """
    Compute fractionation factors.
    """
    
    check_validity(ISO_TYPE, iso_type_list, "ISO_TYPE")
    
    alpha_eq_list = prepare_combined_alpha_eq(
        temp_list, 
        temp_thres_min=temp_thres_min, 
        temp_thres_max=temp_thres_max,
        ISO_TYPE=ISO_TYPE
    )
    
    alpha_kin_list = [
        kin_frac_factor_ice(
            temp, ISO_TYPE=ISO_TYPE
        )
        for temp in temp_list
    ]
    alpha_kin_list = interp_alpha_kin_ice(alpha_kin_list)
    
    alpha_eff_list = [
        alpha_eq_list[i] * alpha_kin_list[i]
        for i, temp in enumerate(temp_list)
    ]

    return {
        "alpha_eq": alpha_eq_list,
        "alpha_kin": alpha_kin_list,
        "alpha_eff": alpha_eff_list,
    }
    

    