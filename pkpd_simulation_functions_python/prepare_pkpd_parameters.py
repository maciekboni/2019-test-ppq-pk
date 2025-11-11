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

def prepare_parameters(num_patients, act_combo, experiment_num, pmax_art_list, ec50_art_list, hill_art_list, pmax_partner_list, ec50_partner_list, hill_partner_list, min_weight = 15, max_weight = 100, weight_step = 5, patient_age = None):

    """
    Prepare parameter combinations for the PKPD simulation.

    Parameters
    ----------
    num_patients : int
        Number of patients for each simulation.
    act_combo : str
        Active drug combination. One of "AL", "DHA-PPQ", "AS-AQ".
    experiment_num : int
        Experiment number for folder naming.
    pmax_art_list : list
        List of pMax values for the ART component. E.g. [0.99, 0.99998]
    ec50_art_list : list
        List of EC50 values for the ART component in ng/microliter. E.g. [0.01, 0.1]
    hill_art_list : list
        List of Hill coefficients for the ART component. E.g. [1, 2]
    pmax_partner_list : list
        List of pMax values for the partner drug. E.g. [0.99, 0.99998]
    ec50_partner_list : list
        List of EC50 values for the partner drug in ng/microliter. E.g. [0.01, 0.1]
    hill_partner_list : list
        List of Hill coefficients for the partner drug. E.g. [1, 2]
    min_weight : int, optional
        Minimum patient weight in kg (default: 15).
    max_weight : int, optional
        Maximum patient weight in kg (default: 100).
    weight_step : int, optional
        Step size for patient weight in kg (default: 5).
    patient_age : list, optional
        List of patient ages corresponding to patient weights in years (only relevant for DHA-PPQ simulations). Default is None. 
        The length of this list should match the number of weight values.
        E.g. if weights are [15, 20, 25], ages could be [1, 2, 3]. 

    Returns
    -------
    tuple
        (weight, params_list, folder, file_suffix)
        - weight : numpy array
            Array of patient weights.
        - params_list : list
            List of parameter combinations for simulations.
        - folder : str
            Name of the folder where results are saved.
        - file_suffix : str
            Suffix for the parameter file.
    """

    # Number of patients for each simulation
    n_patients = num_patients

    # Specify PD parameters for the ART component
    # Make sure they are stored as lists
    pmax_art = list(pmax_art_list)
    print(f"PMax for ART component: {pmax_art}")
    ec50_art = list(ec50_art_list)
    print(f"EC50 for ART component: {ec50_art}")
    hill_art = list(hill_art_list)
    print("Hill Coefficient for ART component: ", hill_art)

    # Specify PD parameters for the partner drug
    # Make sure they are stored as lists
    pmax_partner = list(pmax_partner_list)
    print(f"PMax for partner drug: {pmax_partner}")
    ec50_partner = list(ec50_partner_list)
    print(f"EC50 for partner drug: {ec50_partner}")
    hill_partner = list(hill_partner_list)
    print("Hill Coefficient for partner drug: ", hill_partner)

    # Patient weights
    weight = np.arange(min_weight, max_weight + 1, weight_step)  
    print("Patient weight: ", weight)

    # Prepare output folder 
    experiment = str(experiment_num)
    drug = str(act_combo).replace("-", "_")
    file_suffix = f"model_run_{drug}_{experiment}"
    folder = f"{file_suffix}"
    # Create folder in cwd if it doesn't exist
    os.makedirs(folder, exist_ok=True)  

    if act_combo != "DHA-PPQ":
        patient_age = None  # Age is only relevant for DHA-PPQ simulations
        params_list = []
        for i in range(len(weight)):
            for j in range(len(pmax_art)):
                    for k in range(len(ec50_art)):
                        for l in range(len(hill_art)):
                            for m in range(len(pmax_partner)):
                                for n in range(len(ec50_partner)):
                                    for o in range(len(hill_partner)):
                                        params_list.append([i, weight[i], n_patients, pmax_art[j], ec50_art[k], hill_art[l], pmax_partner[m], ec50_partner[n], hill_partner[o]])
    else:
        patient_age = list(patient_age)
        params_list = []
        for i in range(len(weight)):
            age = patient_age[i]
            for j in range(len(pmax_art)):
                for k in range(len(ec50_art)):
                    for l in range(len(hill_art)):
                        for m in range(len(pmax_partner)):
                            for n in range(len(ec50_partner)):
                                for o in range(len(hill_partner)):
                                    params_list.append([i, weight[i], age, n_patients, pmax_art[j], ec50_art[k], hill_art[l], pmax_partner[m], ec50_partner[n], hill_partner[o]])
    print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Parameter list: {params_list}")

    # Save a copy of the values tested
    with open(f"{folder}/{file_suffix}_parameters.txt", 'w') as file:
        file.write(f"Date and Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        file.write(f"Number of patient(s): {str(n_patients)}\n")
        file.write(f"Weight Value(s) Tested (kg): {str(weight.tolist())}\n")
        file.write(f"pmax_artemether Value(s): {str(pmax_art)}\n")
        file.write(f"ec50_artemether Value(s): {str(ec50_art)}\n")
        file.write(f"Hill Co-efficient(s) Artemether: {str(hill_art)}\n")
        file.write(f"pmax_lum Value(s): {str(pmax_partner)}\n")
        file.write(f"ec50_lum Value(s): {str(ec50_partner)}\n")
        file.write(f"Hill Co-efficient(s) Lumefantrine: {str(hill_partner)}\n")
        if act_combo == "DHA-PPQ":
            file.write(f"Patient Age(s) Tested (years): {str(patient_age)}\n")
    print(f"Parameters saved to {folder}/{file_suffix}_parameters.txt")

    return weight, params_list, folder, file_suffix

#******************************************************************************************************************************************#

if __name__ == "__main__":
    # Empty for now - add test code here if needed later
    pass  # pass is required when there is no other code; python throws an error if you give it an empty code block

#******************************************************************************************************************************************#