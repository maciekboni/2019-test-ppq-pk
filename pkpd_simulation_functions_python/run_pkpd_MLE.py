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

#******************************************************************************************************************************************#
# Artemether-Lumefantrine PK-PD Simulation Function: Log-Likelihood for MLE
#******************************************************************************************************************************************#

def run_simulation_AL_MLE(params, N_obs, k_obs, folder):

    """
    Calculate the log-likelihood for artemether-lumefantrine (AL) efficacy via PK-PD simulation.
    
    Executes C++ simulation and computes log-likelihood based on observed parasitemia clearance data (N_obs, k_obs) for maximum likelihood estimation (MLE).
    
    Parameters
    ----------
    params : list
        List of parameter combinations to simulate.
        (see `prepare_pkpd_parameters.prepare_parameters` for format).
    N_obs : int
        Observed initial parasitemia count from TES data.
    k_obs : float
        Observed parasitemia clearance rate or time-to-clearance metric.
    folder : str
        Name of the folder where results will be stored.
    
    Returns
    -------
    float
        Log-likelihood value for a parameter combination. Higher values indicate better fit to the observed data.
    
    Notes
    -----
    Used internally by MLE optimization routines. Calls external C++ executable via subprocess.
    Requires C++ simulation binary to be compiled and accessible at specified path = `/home/venitha_b/Projects/2019-test-ppq-pk/build/run_ppq_pk`.
    """
     
    N_obs_data = N_obs
    k_obs_data = k_obs
    weight, pmax_art, ec50_art, hill_art, pmax_partner, ec50_partner, hill_partner = params
    
    command_weight = [
        "/home/venitha_b/Projects/2019-test-ppq-pk/build/run_ppq_pk",
        "--AL",
        "-n", '1000',
        "--pmax_artemether", str(pmax_art),
        "--ec50_artemether", str(ec50_art),
        "--hill_artemether", str(hill_art), 
        "--pmax_lum",str(pmax_partner),
        "--ec50_lum", str(ec50_partner), #ng/ml
        "--hill_lum", str(hill_partner), 
        "--weight", str(weight),
        "-o", "1"
    ]
        
    try:
        # Run the command, k and capture output
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
        successful_treatment_count = total_patients - failed_treatment_count
        efficacy = float(successful_treatment_count / total_patients) if total_patients > 0 else 0
        #efficacy_percentage = (100 - ((failed_treatment_count / total_patients) * 100)) if total_patients > 0 else 0

        #log_likelihood = (successful_treatment_count * np.log(efficacy + 1e-9)) + ((total_patients - successful_treatment_count) * np.log((1 - (efficacy + 1e-9))))

        log_likelihood = (k_obs_data * np.log(efficacy + 1e-9)) + ((N_obs_data - k_obs_data) * np.log((1 - (efficacy) + 1e-9)))

        return log_likelihood
    
    except subprocess.CalledProcessError as e:
        print(f"Error in subprocess for weight: {weight}Kg: {e}")
        log_likelihood = -1e9  # Penalize this parameter set
    except Exception as e:
        print(f"An unexpected error occurred for weight: {weight}Kg: {e}")
        log_likelihood = -1e6

#******************************************************************************************************************************************#
# DHA-PPQ PK-PD Simulation Function: Log-Likelihood for MLE
#******************************************************************************************************************************************#

def run_simulation_DHA_PPQ_MLE(params, N_obs, k_obs, folder):
    
    N_obs_data = N_obs
    k_obs_data = k_obs
    weight, age, pmax_art, ec50_art, hill_art, pmax_partner, ec50_partner, hill_partner = params
    
    command_weight = [
        "/home/venitha_b/Projects/2019-test-ppq-pk/build/run_ppq_pk",
        "--DHA-PPQ",
        "-n", '1000',
        "--pmax_dha", str(pmax_art),
        "--ec50_dha", str(ec50_art),
        "--hill_dha", str(hill_art), 
        "--pmax_ppq",str(pmax_partner),
        "--ec50_ppq", str(ec50_partner), #ng/ml
        "--hill_ppq", str(hill_partner), 
        "--weight", str(weight),
        "--age", str(age),
        "-o", "1"
    ]
        
    try:
        # Run the command, k and capture output
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
        successful_treatment_count = total_patients - failed_treatment_count
        efficacy = float(successful_treatment_count / total_patients) if total_patients > 0 else 0
        #efficacy_percentage = (100 - ((failed_treatment_count / total_patients) * 100)) if total_patients > 0 else 0

        #log_likelihood = (successful_treatment_count * np.log(efficacy + 1e-9)) + ((total_patients - successful_treatment_count) * np.log((1 - (efficacy + 1e-9))))

        log_likelihood = (k_obs_data * np.log(efficacy + 1e-9)) + ((N_obs_data - k_obs_data) * np.log((1 - (efficacy) + 1e-9)))

        return log_likelihood
    
    except subprocess.CalledProcessError as e:
        print(f"Error in subprocess for weight: {weight}Kg: {e}")
        log_likelihood = -1e9  # Penalize this parameter set
    except Exception as e:
        print(f"An unexpected error occurred for weight: {weight}Kg: {e}")
        log_likelihood = -1e6

#******************************************************************************************************************************************#
# AS-AQ PK-PD Simulation Function: Log-Likelihood for MLE
#******************************************************************************************************************************************#

def run_simulation_AS_AQ_MLE(params, N_obs, k_obs, folder):
    
    N_obs_data = N_obs
    k_obs_data = k_obs
    weight, pmax_art, ec50_art, hill_art, pmax_partner, ec50_partner, hill_partner = params
    
    command_weight = [
        "/home/venitha_b/Projects/2019-test-ppq-pk/build/run_ppq_pk",
        "--AS-AQ",
        "-n", '1000',
        "--pmax_artesunate", str(pmax_art),
        "--ec50_artesunate", str(ec50_art),
        "--hill_artesunate", str(hill_art), 
        "--pmax_aq",str(pmax_partner),
        "--ec50_aq", str(ec50_partner), #ng/ml
        "--hill_aq", str(hill_partner), 
        "--weight", str(weight),
        "-o", "1"
    ]
        
    try:
        # Run the command, k and capture output
        output_weight = subprocess.run(command_weight, capture_output=True, text=True, check=True, cwd=folder)

    # Process the output
        lines_weight = output_weight.stdout.splitlines()
        df_weight = pd.DataFrame(
            [line.split() for line in lines_weight],
            columns=["PID", "HOUR", "COMP2CONC_AS", "COMP2CONC_AQ", "COMP2CONC_DEAQ", "PARASITEDENSITY"]
        )

        # Clean the DataFrame
        df_weight = df_weight.iloc[1:].apply(pd.to_numeric, errors='coerce')

        # Calculate the efficacy
        failed_treatment_count = df_weight[df_weight['PARASITEDENSITY'] >= 10].shape[0]
        total_patients = df_weight['PID'].nunique()
        successful_treatment_count = total_patients - failed_treatment_count
        efficacy = float(successful_treatment_count / total_patients) if total_patients > 0 else 0
        #efficacy_percentage = (100 - ((failed_treatment_count / total_patients) * 100)) if total_patients > 0 else 0

        #log_likelihood = (successful_treatment_count * np.log(efficacy + 1e-9)) + ((total_patients - successful_treatment_count) * np.log((1 - (efficacy + 1e-9))))

        log_likelihood = (k_obs_data * np.log(efficacy + 1e-9)) + ((N_obs_data - k_obs_data) * np.log((1 - (efficacy) + 1e-9)))

        return log_likelihood
    
    except subprocess.CalledProcessError as e:
        print(f"Error in subprocess for weight: {weight}Kg: {e}")
        log_likelihood = -1e9  # Penalize this parameter set
    except Exception as e:
        print(f"An unexpected error occurred for weight: {weight}Kg: {e}")
        log_likelihood = -1e6

#******************************************************************************************************************************************#
# Negative Log-Likelihood Function for MLE
#******************************************************************************************************************************************#

def neg_log_likelihood(params, n_obs, k_obs, act_combo, folder):
    if folder is None:
        folder = os.getcwd()  # Use current working directory if no folder is specified
    else:
        os.makedirs(folder, exist_ok=True)  # Create folder if it doesn't exist
    if act_combo == "AL":
        log_likelihood = run_simulation_AL_MLE(params, n_obs, k_obs, folder)
    elif act_combo == "DHA-PPQ":
        log_likelihood = run_simulation_DHA_PPQ_MLE(params, n_obs, k_obs, folder)
    elif act_combo == "AS-AQ":
        log_likelihood = run_simulation_AS_AQ_MLE(params, n_obs, k_obs, folder)
    else:
            raise ValueError(f"Unsupported ACT combination: {act_combo}")
    return -log_likelihood

#******************************************************************************************************************************************#

if __name__ == "__main__":
    # Empty for now - add test code here if needed later
    pass  # pass is required when there is no other code; python throws an error if you give it an empty code block

#******************************************************************************************************************************************#