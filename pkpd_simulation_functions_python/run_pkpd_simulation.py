import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.colors as mcolacors
import pickle
from tqdm import tqdm  # For progress bar
from multiprocessing import Pool, cpu_count
import os
import sys
import time
from datetime import datetime 
from datetime import timedelta
import math
import seaborn as sns
import concurrent.futures
from scipy.optimize import minimize
from functools import partial

#******************************************************************************************************************************************#
# Artemether-Lumefantrine PK-PD Simulation Function: Efficacy Calculation
#******************************************************************************************************************************************#

def run_simulation_AL(idx, weight, n_patients, pmax_artemether, ec50_artemether, hill_artemether, pmax_lum, ec50_lum, hill_lum, folder):
    
    """

    Run PK-PD simulation for artemether-lumefantrine (AL).
    
    Executes C++ simulation and calculates the corresponding drug efficacy based on weight dependent pharmacokinetic/pharmacodynamic parameters.
    
    Parameters
    ----------
    idx : int
        Index identifier for this simulation run. 
    weight : float
        Patient weight in kg.
    n_patients : int
        Number of patients of the specified weight to be simulated.
    hill_artemether : float
        Hill coefficient for artemether dose-response curve.
    pmax_artemether : float
        Maximum parasite killing rate for artemether. The value must be between 0 and 1.
    ec50_artemether : float
        EC50 value for artemether in ng/microliter.
    pmax_lum : float
        Maximum parasite killing rate for lumefantrine. The value must be between 0 and 1.
    ec50_lum : float
        EC50 value for lumefantrine in ng/microliter.
    hill_lum : float
        Hill coefficient for lumefantrine dose-response curve.
    folder : str
        Name of the folder where results are saved.

    Returns
    -------
    tuple
        (idx, weight, hill_artemether, pmax_artemether, ec50_artemether, 
         pmax_lum, ec50_lum, hill_lum, efficacy, df_weight)
        
        - idx : int
            Input index identifier.
        - weight : float
            Input patient weight.
        - efficacy : float
            Calculated treatment efficacy (0-1).
        - df_weight : pandas.DataFrame
            Output dataframe containing simulation results of format:
            | PID | HOUR | COMP2CONC_ARTEMETHER | COMP2CONC_LUM | PARASITEDENSITY |
    
    Notes
    -----
    Calls external C++ executable via subprocess. Requires compiled binary at specified path = /home/venitha_b/Projects/2019-test-ppq-pk/build/run_ppq_pk
    
    """

    command_weight = [
        "/home/venitha_b/Projects/2019-test-ppq-pk/build/run_ppq_pk",
        "--AL",
        "-n", str(n_patients),
        "--pmax_artemether", str(pmax_artemether),
        "--ec50_artemether", str(ec50_artemether), # ng/microliter
        "--hill_artemether", str(hill_artemether),
        "--pmax_lum",str(pmax_lum),
        "--ec50_lum", str(ec50_lum), #ng/ml
        "--hill_lum", str(hill_lum), 
        "--weight", str(weight),
        "-o", "1"
    ]

    if folder is None:
                folder = os.getcwd()
    else:
        os.makedirs(folder, exist_ok=True)

    try:
        # Run the command and capture output
        output_weight = subprocess.run(command_weight, capture_output=True, text=True, check=True, cwd=folder)

        # Process the output
        lines_weight = output_weight.stdout.splitlines()
        df_weight = pd.DataFrame(
            [line.split() for line in lines_weight],
            columns=["PID", "HOUR", "COMP2CONC_ARTEMETHER", "COMP2CONC_LUM", "PARASITEDENSITY"]
        )

        # Clean the DataFrame
        df_weight = df_weight.iloc[1:].apply(pd.to_numeric, errors='coerce')

        # Calculate the efficacy
        failed_treatment_count = df_weight[df_weight['PARASITEDENSITY'] >= 10].shape[0]
        total_patients = df_weight['PID'].nunique()
        efficacy = (100 - ((failed_treatment_count / total_patients) * 100)) if total_patients > 0 else 0

        return (idx, weight, pmax_artemether, ec50_artemether, hill_artemether, pmax_lum, ec50_lum, hill_lum, efficacy, df_weight)

    except subprocess.CalledProcessError as e:
        print(f"Error in subprocess for weight: {weight}Kg: {e} and pmax_artemether: {pmax_artemether}")
        return (idx, weight, pmax_artemether, ec50_artemether, hill_artemether, pmax_lum, ec50_lum, hill_lum, np.nan, None)
    except Exception as e:
        print(f"An unexpected error occurred for weight: {weight}Kg: {e} and pmax_artemether: {pmax_artemether}")
        return (idx, weight, pmax_artemether, ec50_artemether, hill_artemether, pmax_lum, ec50_lum, hill_lum, np.nan, None)

#******************************************************************************************************************************************#
# DHA-PPQ PK-PD Simulation Function: Efficacy Calculation
#******************************************************************************************************************************************#

def run_simulation_DHA_PPQ(idx, weight, age, n_patients, pmax_dha, ec50_dha, hill_dha, pmax_ppq, ec50_ppq, hill_ppq, folder):
    
    """
    Run PK-PD simulation for dihydroartemisinin-piperaquine (DHA-PPQ).
    
    Executes C++ simulation and calculates the corresponding drug efficacy based on weight and age-dependent pharmacokinetic/pharmacodynamic parameters.
    
    Parameters
    ----------
    idx : int
        Index identifier for this simulation run.
    weight : float
        Patient weight in kg.
    age : int
        Patient age in years. 
    n_patients : int
        Number of patients of the specified weight and age to be simulated.
    hill_dha : float
        Hill coefficient for DHA dose-response curve.
    pmax_dha : float
        Maximum parasite killing rate for DHA. The value must be between 0 and 1.
    ec50_dha : float
        EC50 value for DHA in ng/microliter.
    pmax_ppq : float
        Maximum parasite killing rate for piperaquine. The value must be between 0 and 1.
    ec50_ppq : float
        EC50 value for piperaquine in ng/microliter.
    hill_ppq : float
        Hill coefficient for piperaquine dose-response curve.
    folder : str
        Name of the folder where results are saved. 
    
    Returns
    -------
    tuple
        (idx, weight, age, hill_dha, pmax_dha, ec50_dha, pmax_ppq, ec50_ppq, 
         hill_ppq, efficacy, df_weight)
        
        - idx : int
            Input index identifier.
        - weight : float
            Input patient weight.
        - age : int
            Input patient age.
        - efficacy : float
            Calculated treatment efficacy (0-1).
        - df_weight : pandas.DataFrame
            Output dataframe containing simulation results of format:
            | PID | HOUR | COMP2CONC_DHA | COMP2CONC_PPQ | PARASITEDENSITY |
    
    Notes
    -----
    Calls external C++ executable via subprocess. Requires compiled binary at specified path = /home/venitha_b/Projects/2019-test-ppq-pk/build/run_ppq_pk
    Age parameter influences piperaquine pharmacokinetics in the simulation.
    
    """
    
    command_weight = [
        "/home/venitha_b/Projects/2019-test-ppq-pk/build/run_ppq_pk",
        "--DHA-PPQ",
        "-n", str(n_patients),
        "--pmax_dha", str(pmax_dha),
        "--ec50_dha", str(ec50_dha), # ng/microliter
        "--hill_dha", str(hill_dha),
        "--pmax_ppq",str(pmax_ppq),
        "--ec50_ppq", str(ec50_ppq), #ng/ml
        "--hill_ppq", str(hill_ppq),
        "--age", str(age), 
        "--weight", str(weight),
        "-o", "1"
    ]

    if folder is None:
                folder = os.getcwd()
    else:
        os.makedirs(folder, exist_ok=True)

    try:
        # Run the command and capture output
        output_weight = subprocess.run(command_weight, capture_output=True, text=True, check=True, cwd=folder)

        # Process the output
        lines_weight = output_weight.stdout.splitlines()
        df_weight = pd.DataFrame(
            [line.split() for line in lines_weight],
            columns=["PID", "HOUR", "COMP2CONC_DHA", "COMP2CONC_PPQ", "PARASITEDENSITY"]
        )

        # Clean the DataFrame
        df_weight = df_weight.iloc[1:].apply(pd.to_numeric, errors='coerce')

        # Calculate the efficacy
        failed_treatment_count = df_weight[df_weight['PARASITEDENSITY'] >= 10].shape[0]
        total_patients = df_weight['PID'].nunique()
        efficacy = (100 - ((failed_treatment_count / total_patients) * 100)) if total_patients > 0 else 0

        return (idx, weight, age, pmax_dha, ec50_dha, hill_dha, pmax_ppq, ec50_ppq, hill_ppq, efficacy, df_weight)

    except subprocess.CalledProcessError as e:
        print(f"Error in subprocess for weight: {weight}Kg: {e} and pmax_dha: {pmax_dha}")
        return (idx, weight, age, pmax_dha, ec50_dha, hill_dha, pmax_ppq, ec50_ppq, hill_ppq, np.nan, None)
    except Exception as e:
        print(f"An unexpected error occurred for weight: {weight}Kg: {e} and pmax_dha: {pmax_dha}")
        return (idx, weight, age, pmax_dha, ec50_dha, hill_dha, pmax_ppq, ec50_ppq, hill_ppq, np.nan, None)

#******************************************************************************************************************************************#
# AS-AQ PK-PD Simulation Function: Efficacy Calculation
#******************************************************************************************************************************************#

def run_simulation_AS_AQ(idx, weight, n_patients, pmax_artesunate, ec50_artesunate, hill_artesunate, pmax_aq, ec50_aq, hill_aq, folder):
    
    """
    Run PK-PD simulation for artesunate-amodiaquine (AS-AQ).
    
    Executes C++ simulation and calculates the corresponding drug efficacy.
    
    Parameters
    ----------
    idx : int
        Index identifier for this simulation run.
    weight : float
        Patient weight in kg.
    n_patients : int
        Number of patients of the specified weight to be simulated.
    hill_artesunate : float
        Hill coefficient for artesunate dose-response curve.
    pmax_artesunate : float
        Maximum parasite killing rate for artesunate. The value must be between 0 and 1.
    ec50_artesunate : float
        EC50 value for artesunate in ng/microliter.
    pmax_amodiaquine : float
        Maximum parasite killing rate for amodiaquine. The value must be between 0 and 1.
    ec50_amodiaquine : float
        EC50 value for amodiaquine in ng/microliter.
    hill_amodiaquine : float
        Hill coefficient for amodiaquine dose-response curve.
    folder : str
        Name of the folder where results are saved. 
    
    Returns
    -------
    tuple
        (idx, weight, hill_artesunate, pmax_artesunate, ec50_artesunate, 
         pmax_amodiaquine, ec50_amodiaquine, hill_amodiaquine, efficacy, df_weight)
        
        - idx : int
            Input index identifier.
        - weight : float
            Input patient weight.
        - efficacy : float
            Calculated treatment efficacy (0-1).
        - df_weight : pandas.DataFrame
            Output dataframe containing simulation results of format:
            | PID | HOUR | COMP2CONC_AS | COMP2CONC_AQ | COMP2CONCDEAQ | PARASITEDENSITY |
    
    Notes
    -----
    Calls external C++ executable via subprocess. Requires compiled binary at specified path = /home/venitha_b/Projects/2019-test-ppq-pk/build/run_ppq_pk
    """
    
    command_weight = [
        "/home/venitha_b/Projects/2019-test-ppq-pk/build/run_ppq_pk",
        "--AS-AQ",
        "-n", str(n_patients),
        "--pmax_artesunate", str(pmax_artesunate),
        "--ec50_artesunate", str(ec50_artesunate), # ng/microliter
        "--hill_artesunate", str(hill_artesunate),
        "--pmax_aq",str(pmax_aq),
        "--ec50_aq", str(ec50_aq), #ng/ml
        "--hill_aq", str(hill_aq), 
        "--weight", str(weight),
        "-o", "1"
    ]

    if folder is None:
                folder = os.getcwd()
    else:
        os.makedirs(folder, exist_ok=True)

    try:
        # Run the command and capture output
        output_weight = subprocess.run(command_weight, capture_output=True, text=True, check=True, cwd=folder)

        # Process the output
        lines_weight = output_weight.stdout.splitlines()
        df_weight = pd.DataFrame(
            [line.split() for line in lines_weight],
        columns=["PID", "HOUR", "COMP2CONC_AS", "COMP2CONC_AQ", "COMP2CONCDEAQ", "PARASITEDENSITY"]
        )

        # Clean the DataFrame
        df_weight = df_weight.iloc[1:].apply(pd.to_numeric, errors='coerce')

        # Calculate the efficacy
        failed_treatment_count = df_weight[df_weight['PARASITEDENSITY'] >= 10].shape[0]
        total_patients = df_weight['PID'].nunique()
        efficacy = (100 - ((failed_treatment_count / total_patients) * 100)) if total_patients > 0 else 0

        return (idx, weight, pmax_artesunate, ec50_artesunate, hill_artesunate, pmax_aq, ec50_aq, hill_aq, efficacy, df_weight)

    except subprocess.CalledProcessError as e:
        print(f"Error in subprocess for weight: {weight}Kg: {e} and hill_artesunate: {hill_artesunate}")
        return (idx, weight, pmax_artesunate, ec50_artesunate, hill_artesunate, pmax_aq, ec50_aq, hill_aq, np.nan, None)
    except Exception as e:
        print(f"An unexpected error occurred for weight: {weight}Kg: {e} and hill_artesunate: {hill_artesunate}")
        return (idx, weight, pmax_artesunate, ec50_artesunate, hill_artesunate, pmax_aq, ec50_aq, hill_aq, np.nan, None)

#******************************************************************************************************************************************#

def run_simulation(act_combo, params_list, folder=None, num_workers=None):

    """
    Run PK-PD simulation in parallel for single or multiple parameter combinations.
    
    Distributes simulation tasks across multiple CPU cores using multiprocessing and aggregates results from individual drug-combination-specific simulation functions.
    
    Parameters
    ----------
    act_combo : str
        Choice of ACT to simulate. One of: `AL` (artemether-lumefantrine), `DHA-PPQ` (dihydroartemisinin-piperaquine) or `AS-AQ` (artesunate-amodiaquine).
    params_list : list
        List of parameter combinations to simulate.
        (see `prepare_pkpd_parameters.prepare_parameters` for format).
    folder : str, optional
        Name of the folder where results are saved. If None (default), uses the current working directory.
    num_workers : int, optional
        Number of CPU cores to use for parallel processing. If None (default), uses all available CPU cores.
      
    Returns
    -------
    A tuple (results, dfs) where:
    - results : list of tuples
        Each tuple contains:
        - Input parameters used for the simulation
        - Calculated efficacy % (0-100)
    - dfs : list of pandas.DataFrame
        Aggregated results from all simulations with columns:
        - `weight` : Patient weight in kg
        - `efficacy` : Treatment efficacy % (0-100)
        - PK/PD parameters (Hill, pMax, EC50) for each drug component
    
    Notes
    -----
    Uses `multiprocessing.Pool` for parallel execution. Progress is displayed via `tqdm` progress bar. 
    
    Examples
    --------
    >>> params = [(0, 20, 2, [0.99], [0.01], [1], [0.99], [0.1], [2])]
    >>> results = run_simulation('AL', params, num_workers=4)
    >>> results = run_simulation('DHA-PPQ', params)
    >>> results = run_simulation('AS-AQ', params, folder='simulation_results')
    """

    # Define number of workers (use all available CPU cores) if not provided
    if num_workers is None:
        num_workers = cpu_count()
    results = []
    dfs = []
    if folder is None:
                folder = os.getcwd()
    else:
        os.makedirs(folder, exist_ok=True)

    if act_combo == "AL":
        with Pool(processes=num_workers) as pool:
            simulation_func = partial(run_simulation_AL, folder=folder)
            with tqdm(total=len(params_list), desc="Running simulations") as pbar:
                for result in pool.starmap(simulation_func, params_list):
                    idx, weight, pmax_artemether, ec50_artemether, hill_artemether, pmax_lum, ec50_lum, hill_lum, efficacy, df_weight = result
                    results.append((idx, weight, hill_artemether, pmax_artemether, ec50_artemether, pmax_lum, ec50_lum, hill_lum, efficacy, df_weight))
                    if df_weight is not None:
                        dfs.append(df_weight)
                        pbar.update(1)

    elif act_combo == "DHA-PPQ":
        with Pool(processes=num_workers) as pool:
            simulation_func = partial(run_simulation_DHA_PPQ, folder=folder)
            with tqdm(total=len(params_list), desc="Running simulations") as pbar:
                for result in pool.starmap(simulation_func, params_list):
                    idx, weight, age, pmax_dha, ec50_dha, hill_dha, pmax_ppq, ec50_ppq, hill_ppq, efficacy, df_weight = result
                    results.append((idx, weight, age, pmax_dha, ec50_dha, hill_dha, pmax_ppq, ec50_ppq, hill_ppq, efficacy, df_weight))
                    if df_weight is not None:
                        dfs.append(df_weight)
                        pbar.update(1)
                        
    elif act_combo == "AS-AQ":
        with Pool(processes=num_workers) as pool:
            simulation_func = partial(run_simulation_AS_AQ, folder=folder)
            with tqdm(total=len(params_list), desc="Running simulations") as pbar:
                for result in pool.starmap(simulation_func, params_list):
                    idx, weight, pmax_artesunate, ec50_artesunate, hill_artesunate, pmax_aq, ec50_aq, hill_aq, efficacy, df_weight = result
                    results.append((idx, weight, pmax_artesunate, ec50_artesunate, hill_artesunate, pmax_aq, ec50_aq, hill_aq, efficacy, df_weight))
                    if df_weight is not None:
                        dfs.append(df_weight)
                        pbar.update(1)
    else:
        raise ValueError(f"Unsupported ACT combination: {act_combo}")
    return results, dfs


#**********************************************************************************************************************************************#

if __name__ == "__main__":
    # Empty for now - add test code here if needed later
    pass  # pass is required when there is no other code; python throws an error if you give it an empty code block

#******************************************************************************************************************************************#

