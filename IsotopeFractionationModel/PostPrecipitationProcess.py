# PostPrecipitationProcesses.py
"""
Post-precipitation processes module.

This module includes functions to calculate the effects of sublimation from 
falling snow into surface air and to update the isotopic composition and 
specific humidity of the air as a result of these processes.
"""
import numpy as np
from .config import *

def post_snowfall(
    q_surf: float, 
    delta_q_surf: float,
    snow: float,
    delta_snow: float,
    resub_factor: float,
    rh_surf: float = 0.7,
    drh: float = 0.1,
    # snow_duration_factor: float = 100,
) -> dict:
    """
    Calculate the effects of sublimation from falling snow on surface air.

    Parameters:
    - q_surf (float): Initial specific humidity of the surface air (g/kg).
    - delta_q_surf (float): Initial isotopic ratio of surface vapor (ex., X'/X - 1 [ND]).
    - snow (float): Specific humidity of snowfall (g/kg).
    - delta_snow (float): Isotopic ratio of the snowfall (ex., X'/X - 1 [ND]).
    - resub_factor (float): Fractionation of snow mass that sublimates into the air.
    - rh_surf (float): Initial relative humidity of the surface air (ND).
    - drh (float): Change in relative humidity of the surface air (ND).
    # - snow_duration factor (float): Tuning parameter for snowfall duration.
    #   Default is 100.

    Methods:
    - Calculate the effective amount of sublimated snowfall (`snow_eff`) by scaling 
      `snow` using `snow_duration_factor`.
    - Update the surface air's specific humidity and isotopic composition using
      the `resublimation` function.

    Returns:
    - dict: Updated conditions including:
        - `delta_snow` (float): Isotopic ratio of snowfall (ex., X'/X - 1 [ND]).
        - `snow` (float): Effective sublimated snow (g/kg).
        - `delta` (float): Updated isotopic ratio of surface vapor (ex., X'/X - 1 [ND]).
        - `q` (float): Updated specific humidity of surface air (g/kg).
    """

    # snow_eff = snow * snow_duration_factor
    snow_eff = q_surf * drh / rh_surf / resub_factor

    delta_q_surf_updated, q_surf_updated = resublimation(
        delta_q_surf, delta_snow, q_surf, snow_eff, resub_factor
    )
    return {
        "delta_snow": delta_snow * 1000,
        "snow": snow_eff,
        "delta":delta_q_surf_updated * 1000, 
        "q": q_surf_updated
    }

def resublimation(
    delta_vapor: float, 
    delta_snow: float, 
    qv: float, 
    snow: float, 
    f: float
) -> tuple:
    """
    Calculate the updated specific humidity and isotopic composition after snow sublimation.
    
    Parameters:
    - delta_vapor (float): Initial isotopic ratio of surrface vapor (ex., X'/X - 1 [ND]).
    - delta_snow (float): Isotopic ratio of the snowfall (ex., X'/X - 1 [ND]).
    - qv (float): Initial specific humidity fo the surface vapor (g/kg).
    - snow (float): Specific humidity of snowfall (g/kg).
    - f (float): Fraction of snow mass that sublimates into the air.

    Methods:
    - Calculate the sublimated snow amount (`qi`) as the product of `snow` and `f`.
    - Compute the updated specific humidity (`q_updated`) by adding `qi` and `qv`.
    - Update the isotopic ratio of vapor (`R_updated`) by weighting the initial
      ratio (`Rv`) and the snow ratio (`Ri`) based on their respective masses.

    Returns:
    - tuple:
        - float: Updated isotopic ratio of surface vapor (ex., X'/X - 1 [ND]).
        - float: Updated specific humidity of surface vapor (g/kg).
    """
    
    qi = np.abs(snow) * f # 降雪のうち昇華により大気に戻る量
    
    Rv = delta_vapor + 1 # 水蒸気の同位体比（レイリー蒸留後)
    Ri = delta_snow + 1 # 降雪の同位体比
    print()
    q_updated = qv + qi # 水蒸気量の更新
    R_updated = (qv * Rv + qi * Ri) / q_updated # 同位体比の更新

    return R_updated - 1, q_updated    