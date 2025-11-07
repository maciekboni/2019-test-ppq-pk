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

def prepare_parameters(num_patients, act_combo, experiment_num, pmax_art_list, ec50_art_list, hill_art_list, pmax_partner_list, ec50_partner_list, hill_partner_list, min_weight=15, max_weight=100, weight_step=5):

    # Define number of workers (use all available CPU cores)
    num_workers = cpu_count()

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

    # Create parameter combinations
    params_list = []
    for i in range(len(weight)):
        for j in range(len(pmax_art)):
                for k in range(len(ec50_art)):
                    for l in range(len(hill_art)):
                        for m in range(len(pmax_partner)):
                            for n in range(len(ec50_partner)):
                                for o in range(len(hill_partner)):
                                    params_list.append([i, weight[i], n_patients, pmax_art[j], ec50_art[k], hill_art[l], pmax_partner[m], ec50_partner[n], hill_partner[o]])

    print(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Parameter list: {params_list}")

    # Prepare output folder 
    experiment = str(experiment_num)
    drug = str(act_combo).replace("-", "_")
    file_suffix = f"model_run_{drug}_{experiment}"
    folder = f"{file_suffix}"

    # Create folder in cwd if it doesn't exist
    os.makedirs(folder, exist_ok=True)  

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
    print(f"Parameters saved to {folder}/{file_suffix}_parameters.txt")

    return weight, params_list, folder, file_suffix