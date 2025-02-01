# config.py

# Constants
temp0_K = 273.15  # Freezing point of water [K]
e0 = 6.112  # Reference saturation vapor pressure [hPa]
Lv = 2.5e6  # Latent heat of vaporization [J/kg]
R_const = 8.314  # Universal gas constant [kJ/(mol·K)]
Mw = 18.015  # Molar mass of water [g/mol]
Ma = 28.964  # Molar mass of air [g/mol]

# Universal parameters
iso_type_list = ["HDO", "H218O"]
phase_type_list = ["vl","vi"]
diffusion_ref_list = ["M78", "C03"]
sat_eq_ref_list = ["Clausius Clapeyron", "Sonntag (1990)"]

iso_label_dict = {
    "HDO":"$\mathsf{\delta D}$",
    "H218O": "$\mathsf{\delta ^{18}O}$"
}

# Export all constants and parameters
__all__ = [
    "temp0_K",
    "e0",
    "Lv",
    "R_const",
    "Mw",
    "Ma",
    "iso_type_list",
    "phase_type_list", 
    "diffusion_ref_list",
    "sat_eq_ref_list",
    "iso_label_dict"
]