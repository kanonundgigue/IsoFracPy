# PostPrecipitationProcesses.py
"""
Post-precipitation processes module.

This module includes functions to compute the impact of snowfall and sublimation
on the isotopic composition and specific humidity of surface air.
"""
import numpy as np
from .config import *
from .RayleighDistillation import rayleigh_step


def snowfall_time_integration(
    delta_ry: float, 
    q_ry: float, 
    alpha: float, 
    snow_dt: float, 
    duration: int
) -> tuple[float, float]:
    """
    Perform time integration to compute the evolution of isotopic composition during 
    snowfall.
    
    Parameters:
    - delta_ry (float): Initial isotopic ratio of cloud cover (ND).
    - q_ry (float): Initial specific humidity of cloud vapor (g/kg).
    - alpha (float): Equilibrium fractionation factor for snow formation.
    - snow_dt (float): Snowfall rate per unit time (g/kg/s).
    - duration (int): Duration of snowfall in seconds.
        
    Returns:
    - tuple:
      - float: Updated isotopic ratio of cloud vapor (ND).
      - float: Mean isotopic ratio of snowfall (ND).
    """
    delta_cloud_vapor = delta_ry
    delta_snow_sum = 0
    
    for _ in range(duration):
        delta_hoge = rayleigh_step(alpha, q_ry, -snow_dt, delta_cloud_vapor)

        # Update cloud vapor isotopic ratio iteratively
        delta_cloud_vapor = (
            delta_hoge * (q_ry - snow_dt) + delta_ry * snow_dt
        ) / q_ry
        
        # Accumulate snowfall isotopic ratio
        delta_snow_sum += ((delta_hoge + 1) * alpha - 1) 
    
    # Compute mean isotopic ratio for snowfall
    delta_snow = delta_snow_sum / duration
    
    return delta_cloud_vapor, delta_snow

def calc_snow_dt(config: dict) -> float:
    """
    Compute snowfall rate per unit time (snow_dt) based on precipitation and air mass.
    
    Parameters:
    - config (dict): Configuration dictionary containing:
      - "prcp_perday" (float): Daily precipitation amount (kg/m2/day).
      - "p_btm" (float): Bottom pressure level (hPa).
      - "p_top" (float): Top pressure level (hPa).
      - "BOOL_RESUB" (bool): Whether resublimation is considered.
      - "resub_factor" (float): Fraction of snow mass that sublimates.
      
    Returns:
    - float: Snowfall rate per unit time (g/kg/s).
    """
    try:
        # Convert precipitation rate to kg/m2/s
        prcp_dt = config["prcp_perday"] / day_per_sec 
        # Air mass in kg/m2
        Q_cld = (config["p_btm"] - config["p_top"]) * 100 / grav
        if config["BOOL_RESUB"]:
            return prcp_dt / (1 - config["resub_factor"]) / Q_cld * 1000 # g/kg/s
        else:
            return prcp_dt / Q_cld * 1000 # g/kg/s
    except KeyError as e:
        raise ValueError(f"Missing required config parameter: {e}")
        
def generate_snowfall(
    delta_ry: float, 
    q_ry: float, 
    alpha: float, 
    config: dict
) -> tuple[float, float]:
    """
    Generate snowfall and compute its isotopic composition.
    
    Parameters:
    - delta_ry (float): Initial isotopic ratio of cloud vapor (ND).
    - q_ry (float): Initial specific humidity of cloud vapor (g/kg).
    - alpha (float): Equilibrium fractionation factor for snow formation.
    - config (dict): Configuration dictionary containing:
      - "prcp_duration" (float): Duration of snowfall in days.
      - "prcp_perday" (float): Daily precipitation amount (kg/m2/day).
      - "p_btm" (float): Bottom pressure level (hPa).
      - "p_top" (float): Top pressure level (hPa).
      - "BOOL_RESUB" (bool): Whether resublimation is considered.
      - "resub_factor" (float): Fraction of snow mass that sublimates.
      
    Returns:
    - tuple:
      - float: Total snowfall (g/kg).
      - float: Mean isotopic ratio of snowfall (ND).
    """
    try:
        duration = int(config["prcp_duration"] * day_per_sec)
    except KeyError as e:
        raise ValueError(f"Missing required prcp_duration parameter in config: {e}")
        
    snow_dt = calc_snow_dt(config)
    snow = snow_dt * duration
    
    delta_cloud_vapor, delta_snow = snowfall_time_integration(
        delta_ry, q_ry, alpha, snow_dt, duration
    )

    return snow, delta_snow

def resublimation(
    delta_vapor: float, 
    delta_snow: float, 
    qv: float, 
    snow: float, 
    f: float
) -> tuple[float, float]:
    """
    Calculate the updated specific humidity and isotopic composition after snow sublimation.
    
    Parameters:
    - delta_vapor (float): Initial isotopic ratio of surrface vapor (ex., X'/X - 1 [ND]).
    - delta_snow (float): Isotopic ratio of the snowfall (ex., X'/X - 1 [ND]).
    - qv (float): Initial specific humidity fo the surface vapor (g/kg).
    - snow (float): Specific humidity of snowfall (g/kg).
    - f (float): Fraction of snow mass that sublimates into the air.

    Returns:
    - tuple:
        - float: Updated isotopic ratio of surface vapor (ex., X'/X - 1 [ND]).
        - float: Updated specific humidity of surface vapor (g/kg).
    """
    
    qsub = np.abs(snow) * f # Snow mass converted to vapor
    
    Rv = delta_vapor + 1 # Initial vapor isotopic ratio
    Ri = delta_snow + 1 # Snow isotopic ratio

    q_updated = qv + qsub # Updated specific humidity
    R_updated = (qv * Rv + qsub * Ri) / q_updated # Updated isotopic ratio

    return R_updated - 1, q_updated    