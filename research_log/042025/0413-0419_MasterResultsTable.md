# 04/13 - 04/19: Master Results Table

### Description: Implementing changes one at a time and noting down the parameters and corresponding efficacy

| Change | `--art` efficacy for n = 1000 | `--AL --pmax_lum = 0` efficacy for n = 1000 | `--AL` efficacy for n = 1000 | Global Parameters | Artemether specific parameters | Dosing Schedule | Lumefantrine specific parameters | Dosing Schedule |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| None from 09/09 |  | 64.2 | 96.9 | PMF = 1.0
AGE = 25
WEIGHT = 54

 | median_weight  =  48.5
patient_blood_volume = 5500000.0
pmax_art = 0.983 | Erroneous, outdated | median_weight  =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Outdated |
| 1. Change output results function into output_art_monotherapy and output_combination_AL

2. Changed G_OUTPUT_TYPE TO G_CLO_OUTPUT_TYPE

3. Edited the output function calls in `--art` and `--AL` | 1.4000000000000057  | 67.5 | 98.2  | PMF = 1.0
AGE = 25
WEIGHT = 54

 | median_weight  =  48.5
patient_blood_volume = 5500000.0
pmax_art = 0.983 | Erroneous, outdated | median_weight  =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Outdated |
| Changed `--pmax_art` = 0.99997 | 100 | - | - | PMF = 1.0
AGE = 25
WEIGHT = 54 | Patient_blood_volume = 5500000.0
pmax_art = 0.99997 | Erroneous, outdated | median_weight  =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Outdated |
| Changing the formula to calculate the PMF in main.cpp (implemented density dependence correctly) |  | 1.9000000000000057 | 92.3 | PMF = 1.0
AGE = 25
WEIGHT = 54 | median_weight  =  48.5
patient_blood_volume = 5500000.0
pmax_art = 0.983 | Erroneous, outdated | median_weight  =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Outdated |
| G_CLO_PMF = 10.0 (from 1.0) |  | 0 | 45.49999999999999  | PMF = 10
AGE = 25
WEIGHT = 54 | median_weight  =  48.5
patient_blood_volume = 5500000.0
pmax_art = 0.983 | Erroneous, outdated | median_weight  =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Outdated |
| Changing artemether dosing schedule, fixing dosing bug |  | 0 | 60.8 | PMF = 10
AGE = 25
WEIGHT = 54 | median_weight  =  48.5
patient_blood_volume = 5500000.0
pmax_art = 0.983 | Fixed! | median_weight  =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Outdated |
| Changing lumefantrine dosing schedule |  | 0 | 58.900000000000006  | PMF = 10
AGE = 25
WEIGHT = 54 | median_weight  =  48.5
patient_blood_volume = 5500000.0
pmax_art = 0.983 | Fixed! | median_weight  =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Fixed! |
| Implementing set_age_and_weight() |  | 0 | 61.7 | PMF = 10
AGE = 25
WEIGHT = 54 | median_weight  =  48.5
patient_blood_volume = 5500000.0* (48.5/54) = 
pmax_art = 0.983 | Fixed | median_weight  =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Fixed! |
| Scaling Hill Equation by patient blood volume. NOTE: EC50 SCALED BY PATIENT BLOOD VOLUME |  | 0 | 56.6 | PMF = 10
AGE = 25
WEIGHT = 54 | median_weight  =  48.5
patient_blood_volume = 5500000.0* (48.5/54) =
pmax_art = 0.983 | Fixed | median_weight  =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Fixed |
| 1) Making median_weight the same in `pkpd_artemether.cpp` (48.5) and `pkpd_lum.cpp` (54) to 54 Kg

2) Editing age and weight implementation |  | 0 | 57.6 
 | PMF = 10
AGE = 25
WEIGHT = 54 | median_weight  =  48.5
patient_blood_volume = 5500000.0
pmax_art = 0.983 | Fixed | median_weight  =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Fixed |
| Moving `initial_log10_parasitemia` from dha constructor to `initialize_params()` |  | 0 | 59.300000000000004  | PMF = 10
AGE = 25
WEIGHT = 54 | median_weight  =  48.5
patient_blood_volume = 5500000.0
pmax_art = 0.983 | Fixed | median_weight  =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Fixed |
| Not scaling EC50 in Hill Equation by `patient_blood_volume` |  | 0 | 0 | PMF = 10
AGE = 25
WEIGHT = 54 | median_weight  =  48.5
patient_blood_volume = 5500000.0
pmax_art = 0.983 | Fixed | =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Fixed |
| Copying class and functions in `pkpd_dha.cpp` to `pkpd_artemether.cpp` |  | 0 | 57.6  | PMF = 10
AGE = 25
WEIGHT = 54 | median_weight  =  48.5
patient_blood_volume = 5500000.0
pmax_art = 0.983 | Fixed | =  54.0;
age = 25.0;
patient_blood_volume = 5500000.0 | Fixed |
|  |  |  |  |  |  |  |  |  |
|  |  |  |  |  |  |  |  |  |