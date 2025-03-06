def set_param(ISO_TYPE: str = "HDO", ALL_PARAM=True):
    sens_dict = set_sensitive_param(ISO_TYPE)
    if ALL_PARAM:
        test_dict = set_test_param(ISO_TYPE)
    else:        
        test_dict = set_fix_param(ISO_TYPE)
        
    return sens_dict | test_dict
    
def set_sensitive_param(ISO_TYPE: str = "HDO"):
    return{
            "resub_factor": [0.2, 0.5],    
            "prcp_duration":[0.25, 1, 3], 
            # "CLOUD_AGING": [True, False],
            "temp_air_fin": [-20, -10],
        "temp_air_init_list": [[5, 10, 15], [0, 5, 10]], 
    }
    
def set_test_param(ISO_TYPE: str = "HDO"):
    return{
            # Post precipitation
            "prcp_perday":[1, 2],
            "p_top": [400, 500],
            # "CLOUD_AGING": [True, False],
           #  # Background surface condition
            "rh_surf": [0.7, 0.8],
            # Rayleigh distillation
            "dt": [0.5, 1],
            "ALPHA_RY_MODE": ["eq", "eff"], 
           #  # Initial evaporation
            # "delta_q_surf": [-120, -100] if ISO_TYPE == "HDO" else [-15, -12.5],
            "h_air": [1, 0.5],
            "surface_wind": [6.5, 35],
            
        }  
    
def set_fix_param(ISO_TYPE: str = "HDO"):
            return {
            # Post precipitation
            "prcp_perday":[2],
            "p_top": [500],
            # "CLOUD_AGING": [False],
           #  # Background surface condition
            "rh_surf": [0.75],
            # Rayleigh distillation
            "dt": [0.5],
            "ALPHA_RY_MODE": ["eff"], 
           #  # Initial evaporation
            # "delta_q_surf": [-120] if ISO_TYPE == "HDO" else [-15],
            "h_air": [1],
            "surface_wind": [6.5],
            # "temp_air_init_list": [[5, 10, 15]],
        }