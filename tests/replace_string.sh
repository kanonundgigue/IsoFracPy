#!/bin/bash
# Script to replace strings
ext="py"
directory="$HOME/modern_model_comparison/scripts/IsoFracPy"

replacements=(
    # "compute_isotope_delta,delta_from_ratio"
    # "compute_saturated_specific_humidity,sat_specific_humidity"
    # "compute_saturated_vapor_mixing_ratio,sat_vapor_mixing_ratio"
    # "compute_equilibrium_fractionation_factor,eq_frac_factor"
    # "compute_kinetic_fractionation_factor_for_seasurface_evaporation,kin_frac_factor_sea_evap"
    # "compute_kinetic_fractionation_factor_for_supersaturation_condensation,kin_frac_factor_supersat"
    # "compute_effective_relative_humidity,effective_humidity"
    # "get_initial_conditions,get_init_conditions"
    # "plot_initial_conditions,plot_init_conditions"
    # "combine_alpha_eqi_and_eql,combine_alpha_eq"
    # "compute_fractionation_factors,prepare_frac_factors"
    # "interpolate_alpha_kin,interp_alpha_kin"
    # "compute_alpha_for_raindrop_evaporation,adjust_alpha_raindrop_evap"
    # "plot_rayleigh_distillation,plot_rayleigh"
    # "compute_post_snowfall,post_snowfall"
#    "compute_isotope_ratio_of_air,initial_sea_evap_fractionation"
#    "kin_frac_factor_supersat,kin_frac_factor_supersat_ice"
    # "interp_alpha_kin,interp_alpha_kin_supersat_ice"
    "compute_rayleigh_process,rayleigh_step"
)


find "$directory" -type f -name "*.${ext}"| while read -r file; do
    for replacement in "${replacements[@]}"; do
        IFS="," read -r search_string replace_string <<< "$replacement"
        # Replace file content
        if grep -q "$search_string" "$file"; then
            sed -i "s/${search_string}/${replace_string}/g" "$file"
            echo "Replaced '$search_string' with '$replace_string' in: $file"
        fi
    done
done

echo "Replacement completed."
