#Date: 2025-04-25

#Objective:
#Let's try to figure out why the model is not producing the same efficacy across all weight groups

#Ideas for how to do this:
# 1. Output the state variables for each weight group and graph them; perhaps I will be able to catch something
# 2. Output the rates that determine how the drug moves across the compartments for each weight group - maybe the issue is with the rates
# 3. Output the killing constant a for each weight group
# 4. :think What is the EMax in our model? Its the pmax right? 

# Pasting these here for easier access

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
import datetime 
from math import *

# Trying out Idea 1: Output the state variables for each weight group and graph them
# First editing main.cpp/artemether.cpp/lum.cpp to output the state variables

