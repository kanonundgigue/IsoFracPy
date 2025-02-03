import numpy as np
from IsotopeFractionationModel import (
    get_init_conditions,
    rayleigh_process,
    plot_q_dq,
    generate_snowfall,
    # post_snowfall,
    resublimation,
    prepare_frac_factors,
    sat_specific_humidity
)

def configure_model():
    """
    Configure model parameters.

    Returns:
    - dict: Dictionary containing model parameters.
    """
    config_dict = {
        "ISO_TYPE": "HDO",  # Isotope type to simulate. Options are "HDO" or "H218O".
        "dt": 0.5,  # Temperature step size for the simulation [°C].
        "h_air": 1,
        "surface_wind": 6.5,
        "temp_min": -40,  # Minimum temperature for the simulation range [°C].
        "temp_max": 10,  # Maximum temperature for the simulation range [°C].
        "temp_thres_min": -20,  # Lower temperature threshold for ice-vapor fractionation ("vi" phase) [°C].
        "temp_thres_max": 0,  # Upper temperature threshold for liquid-vapor fractionation ("vl" phase) [°C].
        #
        ## Tuning parameters
        "temp_sea_init_list": [ 0, 0, 5, 10, 15],  # List of initial sea surface temperatures [°C].
        "temp_air_init_list": [-5, 0, 5, 10, 15],  # List of initial air temperatures above the sea surface [°C].
        "temp_air_fin": -30,  # Final air temperature after Rayleigh distillation [°C].
        "BOOL_REEVAP": False,  # Whether to consider raindrop re-evaporation during Rayleigh distillation.
                              # (Ref. Worden et al., 2007).
        "reevap_factor": 0.7,  # Re-evaporation factor for raindrops. Affects the fractionation adjupsstment.
        "BOOL_RESUB": True,  # Whether to consider snow sublimation as a post-precipitation process.
        "resub_factor": 0.5,  # Sublimation fraction factor. Controls the effect of snow sublimation on isotope ratios.
        "prcp_perday": 2, # Daily precipitation amount 
        "prcp_duration": 1, # Duration of precipitation (Number of day).
        "p_btm": 700, # Air pressure of cloud bottom [hPa].
        "p_top": 400, # Air pressure of cloud top [hPa].
        
        # "snow_duration_factor": 100,  # Parameter tuning the duration of snowfall effects.
        "delta_q_surf": -120,  # Isotopic ratio (δ-value) of surface vapor [‰].
        "temp_surf": 0,  # Temperature of the surface air layer [°C].
        "rh_surf": 0.75,  # Relative humidity of the surface air layer [ND].
        "drh": 0.2,  # Change in relative humidity of the surface air layer [ND].
        "ALPHA_RY_MODE": "eff",  # Mode for fractionation factors: "eff" for effective fractionation 
                             # (including kinetic effects during supercooling), 
                             # "eq" for equilibrium fractionation only.
    }

    # Generate a list of default temperatures based on the specified range and step size
    config_dict["temp_default_list"] = np.arange(
        config_dict["temp_min"], config_dict["temp_max"], config_dict["dt"]
    )
    q_sat = sat_specific_humidity(config_dict["temp_surf"])
    config_dict["q_surf"] = config_dict["rh_surf"] * q_sat
    # Specific humidity at the final surface air [g/kg].

    return config_dict
    
def main():
    """
    Main function for running the Rayleigh distillation model.

    Methods:
    - Configure model parameters.
    - Calculate initial conditions for sea-air interaction.
    - Compute fractionation factors across temperature ranges.
    - Perform Rayleigh distillation and post-precipitation processes.
    - Plot the results of the model simulation.

    Returns:
    - None: Displays the plot.
    """

    config = configure_model()
    
    # Calculate initial conditions
    initial_dict = initialization(config)

    # Calculate fractionation factors for the given temperature ranges
    frac_factors_dict, alpha_ry_mode_list = get_fractionation_factors(config)

    rayleigh_results_dict, post_precipitation_results_dict = process_vapor_isotopes(
        config, initial_dict, frac_factors_dict, alpha_ry_mode_list
    )
        
    # Plot results
    title = f"Rayleigh distillation model" 
    title += f" (f={config['resub_factor']})" if config["BOOL_RESUB"] else f" (No resublimation.)"
    
    plot_q_dq(
        config, rayleigh_results_dict, post_precipitation_results_dict,
        ISO_TYPE="HDO", title=title
    )

def initialization(config: dict):
    """
    Calculate initial conditions for sea-air interaction.

    Parameters:
    - config (dict): Configuration dictionary containing model parameters.

    Returns:
    - dict: Initial conditions including specific humidity and isotope ratios.
    """
    return get_init_conditions(
        config["temp_sea_init_list"],
        config["temp_air_init_list"],
        ISO_TYPE=config["ISO_TYPE"],
        h_air=config["h_air"],
        surface_wind=config["surface_wind"]
    )

def get_fractionation_factors(config: dict):
    """
    Calculate fractionation factors for the temprature range.

    Parameters:
    - config (dict): Configuration dictionary containing model parameters.

    Returns:
    - tuple:
        - dict: Fractionation factors dictionary.
        - list: Selected alpha_ry_mode_list based on ALPHA_RY_MODE.
    """
    frac_factors_dict = prepare_frac_factors(
        config["temp_default_list"], 
        temp_thres_min=config["temp_thres_min"],
        temp_thres_max=config["temp_thres_max"],
        ISO_TYPE=config["ISO_TYPE"]
    )
    alpha_ry_modes = {
        "eq": frac_factors_dict["alpha_eq"],
        "eff": frac_factors_dict["alpha_eff"]
    }
    alpha_ry_mode_list = alpha_ry_modes.get(
        config["ALPHA_RY_MODE"], frac_factors_dict["alpha_eff"]
    )

    return frac_factors_dict, alpha_ry_mode_list

def process_vapor_isotopes(
    config: dict, 
    initial_dict: dict,
    frac_factors_dict: dict, 
    alpha_ry_mode_list: list
) -> tuple:
    """
    Perform Rayleigh and post precipitation processes.

    Parameters:
    - config (dict): Configuration dictionary containing model parameters.
    - initial_dict (dict): Initial conditions.
    - frac_factors_dict (dict): Fractionation factors.
    - alpha_ry_mode_list (list): Selected fractionation factors.

    Returns:
    - tuple:
        - dict: Results from Rayleigh distillation.
        - dict: Results from post-precipitation process.
    """
    rayleigh_results_dict = {
        temp_air_init: rayleigh_process(
            temp_air_init, 
            config["temp_air_fin"], 
            initial_dict["q_sat_air"][i], 
            initial_dict["delta_air"][i] / 1000,  # permil -> ratio
            alpha_ry_mode_list,            
            frac_factors_dict["alpha_kin"], 
            config["temp_default_list"], 
            BOOL_REEVAP=config["BOOL_REEVAP"],
            reevap_factor=config["reevap_factor"], 
            dt=config["dt"],
            )
        for i, temp_air_init in enumerate(config["temp_air_init_list"])
    }

    post_precipitation_results_dict = {
        temp_air_init: perform_post_precipitation(
            config, rayleigh_results, alpha_ry_mode_list,
        ) for temp_air_init, rayleigh_results in rayleigh_results_dict.items()
    }
        
    return rayleigh_results_dict, post_precipitation_results_dict

def perform_post_precipitation(config, rayleigh_results, alpha_ry_mode_list):
    """
    Perform the post-precipitation process (e.g., snow sublimation).

    Parameters:
    - config (dict): Configuration dictionary containing model parameters.
    - rayleigh_results (dict): Results from Rayleigh distillation.
    - alpha_ry_mode_list (list): Selected fractionation factors.

    Returns:
    - dict: Results from the post-precipitation process.
    """
    alpha_final_idx = np.nanargmin(
        np.abs(config["temp_default_list"] - rayleigh_results["temp"][-1])
    )
    alpha_final = alpha_ry_mode_list[alpha_final_idx]
    q_final = rayleigh_results["q"][-1] # g/kg
    delta_final = rayleigh_results["delta"][-1] / 1000 # permil -> ND
    
    snow, delta_snow = generate_snowfall(delta_final, q_final, alpha_final, config)
    
    if not config["BOOL_RESUB"]:
        return {
            "delta_snow": delta_snow * 1000,
            "snow": snow,
            "delta": config["delta_q_surf"], 
            "q": config["q_surf"]
        }
    
    # Resublimation process
    delta_q_surf_updated, q_surf_updated = resublimation(
        config["delta_q_surf"] / 1000,
        delta_snow,
        config["q_surf"],
        snow,
        config["resub_factor"]
    )        
        
    return {
        "delta_snow": delta_snow * 1000,
        "snow": snow,
        "delta":delta_q_surf_updated * 1000, 
        "q": q_surf_updated
    }
        
if __name__ == "__main__":
    main()        

        