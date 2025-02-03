# PostPrecipitationProcesses.py
"""
Post-precipitation processes module.

This module includes functions to calculate the effects of sublimation from 
falling snow into surface air and to update the isotopic composition and 
specific humidity of the air as a result of these processes.
"""
import numpy as np
from .config import *
from .RayleighDistillation import rayleigh_step


def snowfall_time_integration(delta_final, q_final, alpha_final, snow_dt, duration):
    delta_cloud_vapor = delta_final
    delta_snow = 0
    for loop in range(duration):
        hoge = rayleigh_step(alpha_final, q_final, snow_dt, delta_cloud_vapor)

        delta_cloud_vapor = (hoge * (q_final - snow_dt) + delta_final * snow_dt) / q_final
        delta_snow = delta_snow + ( (hoge + 1) * alpha_final - 1) / duration
    return delta_cloud_vapor, delta_snow

def calc_snow_dt(config):
    prcp_dt = config["prcp_perday"] / day_per_sec # kg/m2/s
    airmass = (config["p_btm"] - config["p_top"]) * 100 / grav # kg/m2        
    if config["BOOL_RESUB"]:
        return prcp_dt / (1 - config["resub_factor"]) / airmass * 1000 # dq/dt; g/kg/s
    else:
        return prcp_dt / airmass * 1000 # dq/dt; g/kg/s

def generate_snowfall(delta_final, q_final, alpha_final, config):
    duration = int(config["prcp_duration"] * day_per_sec)
    snow_dt = calc_snow_dt(config)
    delta_cloud_vapor, delta_snow = snowfall_time_integration(delta_final, q_final, alpha_final, snow_dt, duration)
    snow = snow_dt * duration

    return snow, delta_snow

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
    
    qsub = np.abs(snow) * f # 降雪のうち昇華により大気に戻る量
    
    Rv = delta_vapor + 1 # 水蒸気の同位体比（レイリー蒸留後)
    Ri = delta_snow + 1 # 降雪の同位体比

    q_updated = qv + qsub # 水蒸気量の更新
    R_updated = (qv * Rv + qsub * Ri) / q_updated # 同位体比の更新

    return R_updated - 1, q_updated    