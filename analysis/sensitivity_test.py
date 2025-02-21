import os, pathlib, sys, warnings
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import colormaps
import pandas as pd
import itertools
import concurrent.futures

# Original modules
from save_figure_with_confirmation import save_figure_with_confirmation as save_fig

# path_for_repo=f"/data37/kanon/modern_model_comparison/scripts/IsoFracPy" # server
path_for_repo=f"/Users/kanon/Work/IsoFracPy" # local
iso_model_dir = pathlib.Path(f"{path_for_repo}").resolve()
sys.path.append(str(iso_model_dir))

from main import (
    configure_model, 
    initialization, 
    get_fractionation_factors,
    process_vapor_isotopes
)
from IsotopeFractionationModel import (
    plot_q_dq,
    sat_specific_humidity,
)

def process_combination(
    index: int, 
    combination: tuple, 
    param_test_dict: dict, 
    param_fix_dict: dict, 
    results_dir: str, 
    title: str, 
    GENERATE_FIG: bool, 
    fig=None, 
    gs=None, 
    subplot_hnum_max=None, 
    xlim=None, 
    ylim=None, 
    num_of_analysis=None,
    GENERATE_CSV = True
):
    # try:
    # Configure the model
    config = configure_model_with_params(param_fix_dict, param_test_dict, combination)    
    # Execute the model
    rayleigh_results_dict, post_precipitation_results_dict = execute_model(config)    
    
    if GENERATE_CSV:
        # Save results to CSV
        write_results_csv(
            config, list(param_test_dict.keys()), combination,
            rayleigh_results_dict, post_precipitation_results_dict, 
            results_dir, title
        )


    # except Exception as e:
    #     print(f"Error processing combination {index}: {e}")

    if GENERATE_FIG:
        # Plot results
        plot_results_for_combination(
            config, rayleigh_results_dict, post_precipitation_results_dict,
            index, combination, param_test_dict.keys(), 
            fig, gs, subplot_hnum_max, xlim, ylim, num_of_analysis
        ) 
def run_sensitivity_analysis(
    param_test_dict: dict,
    param_fix_dict: dict = None,
    xlim: tuple = (0, 10),
    ylim: tuple = (-300, -50),
    fig_size_unit: tuple = (5, 5),
    subplot_hnum_max: int = 4,
    title: str = "tuning",
    iso_model_dir: str = ".",
    results_dir: str = "analysis", 
    GENERATE_FIG: bool = True,
    BOOL_SAVE_FIG: bool = True,    
    GENERATE_CSV: bool = True,
    notebook_name: str = None,
    max_workers: int = 16,
):    
    """
    Do sensitivity analysis by varing parameters and plotting results.

    Parameters:
    - param_test_dict (dict): Dictionary of parameters to vary with their values.
    - param_fix_dict (dict, optional): Dictionary of fixed parameters. Default is None.
    - xlim (tuple): x-axis limits for the plots. Default is (0, 10).
    - ylim (tuple): y-axis limits for the plots. Default is (-300, -50).
    - fig_size_unit (tuple): Unit size for each subplot (width, height). Default is (5, 5).
    - subplot_hnum_max (int): Maximum number of subplots in a row. Default is 4.
    - title (str): Title for the entire figure and table. Default is "tuning".
    - iso_model_dir (str): Base directory for the isotope model. Default is ".".
    - results_dir (str): Subdirectory to save analysis results. Default is "analysis".
    - BOOL_SAVE_FIG (bool): Whether save figure or not. Default is True.
    - notebook_name (str): Name of executing notebook name. Default is "".

    Returns:
    - None: The function directly displays the plots and saves the figure.
    """
    param_combinations, num_of_analysis = prepare_parameter_combinations(
        param_test_dict
    )

    if GENERATE_FIG:
        num_of_plot_x, num_of_plot_y = calculate_subplot_grid(num_of_analysis, subplot_hnum_max)
        fig, gs = initialize_figure(title, fig_size_unit, num_of_plot_x, num_of_plot_y)
    else:
        fig, gs = None, None

    if GENERATE_CSV:
        _check_existing_tables(results_dir, title)

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                process_combination, index, combination, param_test_dict, param_fix_dict, results_dir, title, GENERATE_FIG, 
                fig, gs, subplot_hnum_max, xlim, ylim, num_of_analysis, GENERATE_CSV
            )
            for index, combination in enumerate(param_combinations)
        ]
        for future in concurrent.futures.as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Error in combination {future}: {e}")
    if GENERATE_FIG:
        # Finalize and save figure
        finalize_and_save_figure(
            fig, title, iso_model_dir, results_dir, BOOL_SAVE_FIG, notebook_name
        )

    print("Sensitivity analysis completed.")

def _check_existing_tables(results_dir, title):
    # Set results directory path
    results_path = set_results_path(iso_model_dir, results_dir, results_type="tables")
    # Remove existing CSV files in results directory.
    csv_dir = os.path.join(results_path, "tables")
    existing_csv_files = [f for f in os.listdir(csv_dir) if f.endswith(".csv")]

    if existing_csv_files:
        print("The following CSV files already exist:")
        for file in existing_csv_files:            
            if title in file:
                print(f"- {file}")
                confirm = input ("Delete? If not, it may be overwritten. (y/n)").strip().lower()
                if confirm in ["y", "yes"]:
                    os.remove(os.path.join(csv_dir, file))
                else:
                    print("Keep this file.")


def get_simulation_results_dict(
    config: dict, 
    rayleigh_results_dict: dict, 
    post_precipitation_results_dict: dict, 
    temp: float
) -> dict:
    """
    Retrieve simulation results for a given temperature.

    Parameters:
    - config (dict): Configuration.
    - rayleigh_results_dict (dict): Dictionary storing Rayleigh model results.
    - post_precipitation_results_dict (dict): Dictionary storing post-precipitation model results.
    - temp (float): Temperature key for retrieving results.
    
    Returns:
    - dict: Dictionary containing simulation results for the given temperature.
    """
    return {
        "ry_ini_q": rayleigh_results_dict[temp]["q"][0],
        "ry_ini_delta": rayleigh_results_dict[temp]["delta"][0],
        "ry_fin_q": rayleigh_results_dict[temp]["q"][-1],
        "ry_fin_delta": rayleigh_results_dict[temp]["delta"][-1],
        "snow_q": post_precipitation_results_dict[temp]["snow"],
        "snow_delta": post_precipitation_results_dict[temp]["delta_snow"],
        "obs_fin_q": post_precipitation_results_dict[temp]["q"],
        "obs_fin_delta": post_precipitation_results_dict[temp]["delta"],
        "obs_ini_q": config["q_surf"],
        "obs_ini_delta": config["delta_q_surf"]        
    }

    return results_dict
    
def write_results_csv(
    config: dict, 
    param_list: list, 
    combination: list, 
    rayleigh_results_dict: dict,
    post_precipitation_results_dict: dict, 
    results_dir: str,
    title: str,
) -> None:
    """
    Write simulation results to CSV files.

    Parameters:
    - config (dict): Configuration containing surface initial conditions.
    - param_list (list): List of parameter names.
    - combination (list): List of parameter values.
    - rayleigh_results_dict (dict): Dictionary storing Rayleigh model results.
    - post_precipitation_results_dict (dict): Dictionary storing post-precipitation model results.
    - results_dir (str): Directory path to save ressults.
    - title (str):

    Returns:
    - None
    """
    base_columns = [
        "ry_ini_q", "ry_ini_delta", "ry_fin_q", "ry_fin_delta", 
        "snow_q", "snow_delta", "obs_fin_q", "obs_fin_delta"
    ]    

    for i, temp in enumerate(rayleigh_results_dict.keys()):
        results_dict = get_simulation_results_dict(
            config, rayleigh_results_dict, post_precipitation_results_dict, temp
        )
        # Add parameters from combination
        additional_columns = []
        for param, val in zip(param_list, combination):
            if param == "temp_air_init_list":
                results_dict["temp_air_init"] = val[i]
                additional_columns.append("temp_air_init")
            elif param == "temp_sea_init_list":
                results_dict["temp_sea_init"] = val[i]
                additional_columns.append("temp_sea_init")
            # elif param not in ["q_surf", "delta_q_surf"]:
            else:
                results_dict[param] = val
                additional_columns.append(param)  
        if "temp_sea_list" not in additional_columns:
            results_dict["temp_sea_init"] =  config["temp_sea_init_list"][i]
            additional_columns.append("temp_sea_init")        
        
        # Construct DataFrame 
        columns = base_columns + additional_columns
        results_df = pd.DataFrame([results_dict], columns=columns)

        # Set output file path
        results_path = set_results_path(iso_model_dir, results_dir, results_type="tables")
        csv_filename = f"{title}_{i}.csv"
        csv_path = f"{results_path}/tables/{csv_filename}"

        # Save CSV file
        mode, header = ("a", False) if os.path.exists(csv_path) else ("w", True)
        if mode == "w":
            print(f"New file is created:")
            print(f"- {csv_path}") 
        results_df.to_csv(csv_path, index=False, header=header, mode=mode)            

        
def execute_model(config):
    initial_dict = initialization(config)
    frac_factors_dict, alpha_mode_list = get_fractionation_factors(config)

    rayleigh_results_dict, post_precipitation_results_dict = process_vapor_isotopes(
        config, initial_dict, frac_factors_dict, alpha_mode_list
    )
    return rayleigh_results_dict, post_precipitation_results_dict
    
def prepare_parameter_combinations(param_test_dict: dict) -> tuple:
    """
    Prepare parameter combinations for sensitivity analysis.

    Parameters:
    - param_test_dict (dict): Dictionary of parameters to vary with their values.

    Returns:
    - tuple:
        - list: List of parameter combinations.
        - int: Number of combinations.
    """    
    param_combinations = list(itertools.product(*param_test_dict.values()))
    return param_combinations, len(param_combinations)

def calculate_subplot_grid(num_of_analysis: int, subplot_hnum_max: int) -> tuple:
    """
    Calculate the grid size for subplots.

    Parameters:
    - num_of_analysis (int): Total number of analyses.
    - subplot_hnum_max (int): Maximum number of subplots in a row.

    Returns:
    - tuple: Number of plots along x and y axes.
    """
    num_of_plot_x = subplot_hnum_max if num_of_analysis >=subplot_hnum_max else num_of_analysis 
    num_of_plot_y = num_of_analysis//subplot_hnum_max + 1 if num_of_analysis >subplot_hnum_max else 1
    return num_of_plot_x, num_of_plot_y

def initialize_figure(
    fig_title: str, 
    fig_size_unit: tuple, 
    num_of_plot_x: int,
    num_of_plot_y: int,
):
    """
    Initialize the figure for sensitivity analysis.

    Parameters:
    - title (str): Title for the entire figure.
    - fig_size_unit (tuple): Unit size for each subplot (width, height).
    - num_of_plot_x (int): Number of subplots in a row.
    - num_of_plot_y (int): Number of subplots in a column.

    Returns:
    - tuple:
        - matplotlib.figure.Figure: Initialized figure object.
        - matplotlib.gridspec.GridSpec: GridSpec object for subplots.
    """        
    figsize = calculate_figure_size(fig_size_unit, num_of_plot_x, num_of_plot_y)
    fig = plt.figure(layout="tight", figsize=figsize)
    fig.suptitle(fig_title, y=1.02)
    gs = GridSpec(num_of_plot_y, num_of_plot_x, figure=fig)
    return fig, gs

def calculate_figure_size(fig_size_unit, num_of_plot_x, num_of_plot_y):
    """
    Calculate figure size dynamically.
    """
    return (fig_size_unit[0] * num_of_plot_x, fig_size_unit[1] * num_of_plot_y)
    
def configure_model_with_params(param_fix_dict: dict, param_test_dict: dict, combination: tuple) -> dict:
    """
    Configure the model with the given parameter combination.

    Parameters:
    - param_fix_dict (dict): Fixed parameters.
    - param_test_dict (dict): Test parameters to vary.
    - combination (tuple): Specific parameter combination.

    Returns:
    - dict: Configuration dictionary for the model.    
    """
    config = configure_model()
    if param_fix_dict:
        config.update(param_fix_dict)
        
    config.update(dict(zip(param_test_dict.keys(), combination)))

    return config

def generate_subplot_label(index):
    """
    Generate subplot labels like (a1).
    """
    alphabet = "abcdefghijklmnopqrstuvwxyz"
    prefix = alphabet[index % 26]
    suffix = str(index // 26 + 1) if index >= 26 else ""
    return f"{prefix}{suffix}"

def plot_results_for_combination(
    config: dict, 
    rayleigh_results_dict: dict, 
    post_precipitation_results_dict: dict,
    index: int, 
    combination: tuple, 
    param_names, 
    fig, 
    gs, 
    subplot_hnum_max: int, 
    xlim: tuple, 
    ylim: tuple, 
    num_of_analysis: int
):
    """
    Plot results for a specific parameter combination.

    Parameters:
    - config (dict): Configuration dictionary for the model.
    - index (int): Index of the current combination.
    - combination (tuple): Current parameter combination.
    - param_names (iterable): Names of parameters being tested.
    - fig (matplotlib.figure.Figure): Figure object.
    - gs (matplotlib.gridspec.GridSpec): GridSpec object for subplots.
    - subplot_hnum_max (int): Maximum number of subplots in a row.
    - xlim (tuple): x-axis limits for the plot.
    - ylim (tuple): y-axis limits for the plot.
    - num_of_analysis (int): Total number of analyses.

    Returns:
    - None: Modifies the plot directly.
    """      



    title = f"({generate_subplot_label(index)})\n"+ "\n".join(
            [f"{param}={val}" for param, val in zip(param_names, combination)]
        )

    ax = fig.add_subplot(gs[index // subplot_hnum_max, index % subplot_hnum_max])
    ax.set_title(title, loc="left")                
    plot_q_dq(
        config,
        rayleigh_results_dict, 
        post_precipitation_results_dict, 
        ISO_TYPE=config["ISO_TYPE"], 
        title=title, 
        xlim=xlim, 
        ylim=ylim, 
        num_of_subplot=num_of_analysis, 
        ax=ax)

    plt.legend()

def set_results_path(iso_model_dir: str, results_dir: str, results_type: str = "figures") -> str:
    results_path = os.path.join(iso_model_dir, results_dir)
    os.makedirs(os.path.join(results_path, results_type), exist_ok=True)

    return results_path

def finalize_and_save_figure(fig, fig_title, iso_model_dir, results_dir, BOOL_SAVE_FIG, notebook_name):    
    """
    Finalize the figure layout and save it to the results directory.

    Parameters:
    - fig (matplotlib.figure.Figure): Figure object.
    - fig_title (str): Title for the figure.
    - iso_model_dir (str): Base directory for the isotope model.
    - results_dir (str): Subdirectory to save analysis results.
    - BOOL_SAVE_FIG (bool): Whether save figure or not. Default is True.
    - notebook_name (str): Name of this notebook.

    Returns:
    - None: Saves the figure to disk if `BOOL_SAVE_FIG` is True.
    """     
    # plt.tight_layout(rect=[0, 0, 1, 0.98]) # [left, bottom, right, top]
    # plt.tight_layout()
    plt.show()

    results_path = set_results_path(iso_model_dir, results_dir, results_type="figures")
    
    if BOOL_SAVE_FIG:
        save_fig(
            fig, 
            results_dir=results_path,
            figure_name=fig_title.replace("}", "").replace("{", "").replace("$", "").replace(",", "").replace(" ", "_"),
            notebook_path=f"{os.path.abspath(notebook_name)}.ipynb",
            check_overwrite=False
        )