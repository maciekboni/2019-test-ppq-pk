# 04/13 - 04/19

### Sunday, 04/13

### **Main goal:**

Revert back changes to when model was calibrated to AL 96.5% efficacy/ DHA 69.5% efficacy for patients = 54 kg to see what new parameter/function is affecting the efficacy to drop

### **Objective 1: Reverted changes made to main.cpp, pkpd_artemether.cpp and pkpd_lum.cpp to 09092024**

- [ ]  Created folder `research_log` to store research logs
    - Created weekly e-mail reminder to update research logs
- [ ]  Created a copy of main.cpp, pkpd_artemether.cpp, pkpd_lum.cpp and corresponding header files. Tracking changes via Github
    - Location: /home/venitha_b/Projects/2019-test-ppq-pk/ pkpd_model_changes/04132025_backup
- [ ]  Reverted changes to `main.cpp` , `pkpd_artemether.cpp` and `pkpd_lum.cpp` and corresponding header files to 09/09/2024
    - Copy of the 09092024 cpp files are in /home/venitha_b/Projects/2019-test-ppq-pk/ pkpd_model_changes/09092024_backup

### **Objective 2: Testing if `--art` and `--AL` modules give 69.5% / 96.5% efficacy for weight = 54 kg after revision**

- [ ]  Making sure that the `--art` and `--AL` modules do indeed give 69.5%/96.5% efficacies as previously calibrated
    - Analysis is in /home/venitha_b/Projects/2019-test-ppq-pk/plots/AL/parameter_analysis/04132025
    - Results:
        - AL: n = 10k, pmax_artemether = 0.983
            
            `Length of simulation:  0:04:37.493341 
            The failed treatment count for weight 54kg is: 249
            The percentage efficacy for weight 54kg is: 97.51`
            
        - Artemether monotherapy: n = 10k, pmax_artemether = 0.983
            - With `--AL` and `--pmax_lum = 0`
                
                `Length of simulation:  0:04:28.233104 
                The failed treatment count for weight 54kg is: 3171
                The percentage efficacy for weight 54kg is: 68.28999999999999`
                
            - With `--art` : TO DO

### Objective 3: Debugging the dosage schedule of `pkpd_artemether.cpp` and observing `--art` and `--AL` efficacy for weight = 54 kg

- Results:
    - Artemether monotherapy: n = 10k, pmax_artemether = 0.983
        - With `--AL` and `--pmax_lum = 0`
        - `Today's date:  2025-04-13
        Length of simulation:  0:03:48.042578
        The failed treatment count for weight 54kg is: 0
        The percentage efficacy for weight 54kg is: 100.0`
    - AL: n = 10k, pmax_artemether = 0.983
        - Weight = 54 Kg
            - `Today's date:  2025-04-13
            Length of simulation:  0:03:54.101743
            The failed treatment count for weight 54kg is: 0
            The percentage efficacy for weight 54kg is: 100.0`
        - Weight = 12 Kg
            - `Today's date:  2025-04-13 
            Length of simulation:  0:03:47.829280 
            The failed treatment count for weight 12kg is: 0
            The percentage efficacy for weight 12kg is: 100.0`

### Objective 4: Made the following changes:

- `pkpd_artemether.cpp` :
    - Changed constructor to `pkpd_artemether::pkpd_artemether(patient age, patient weight)`
    - Moving log10 parasitaemia calculation from constructor to `initialize_params()`
    - Adding `set_age_and_weight()`
    - Creating `initialize()`
    - Removed `generate_recommended_dosing_schedule()` from constructor
- `pkpd_lum.cpp` :
    - Changed constructor to `pkpd_lum::pkpd_lum(patient age, patient weight)`
    - Changed `median_weight` from 54.0 to 42.0
    - Adding `set_age_and_weight()`
- `main.cpp` :
    - Added `pkpd_artemether::stochastic = true`
    - Edited the constructor call for `therapy_lumefantrine`  and `therapy_AL`
    - Set `G_CLO_PMF` from 1.0 to 10.0

- [ ]  Checking AL efficacy after the above changes were made:
    - Results
        - Weight = 54 Kg
            - `Today's date:  2025-04-13
            Length of simulation:  0:03:51.179762
            The failed treatment count for weight 54kg is: 3631
            The percentage efficacy for weight 54kg is: 63.690000000000005`
- [ ]  Reverted `median_weight` in `pkpd_lum.cpp` from 42.0 to 54.0
    - Results
        - [ ]  Weight = 54 Kg
            - `Today's date:  2025-04-14
            Length of simulation:  0:03:51.406475
            The failed treatment count for weight 54kg is: 3709
            The percentage efficacy for weight 54kg is: 62.91`

### Observation 1: Changing `median_weight` in `pkpd_lum.cpp` from 54.0 to 42.0 did not cause the drop in efficacy

- Will keep `median_weight` as 54.0 as its the same value in `pkpd_artemether.cpp`

### Observation 2: Addition of `set_age_and_weight()` , i.e., scaling the `patient_blood_volume` by `patient_weight` seems to have caused the massive drop in efficacy

- Before checking if scaling by weight did indeed cause the drop in efficacy, let’s see if modifying the dosage schedule of lumefantrine affects the efficacy in any way
    - Results
        - Weight = 54 Kg
            - `Today's date:  2025-04-14 
            Length of simulation:  0:03:54.339982 
            The failed treatment count for weight 54kg is: 3582
            The percentage efficacy for weight 54kg is: 64.18`

### Tuesday, 04/15 to Wednesday, 04/16

### Re-running the simulations but only reverting one change at a time and noting it down in a master results table.

[04/13 - 04/19: Master Results Table ](04%2013%20-%2004%2019%20Master%20Results%20Table%201d6125b3f5d3801daecdc13daf56a8ba.md)

- [ ]  Make sure the simulation gives you the old efficacies for weight = 54 Kg
    - [ ]  `--AL`, `--pmax_lum = 0`
        
        Today's date:  2025-04-15
        Length of simulation:  0:00:28.841660
        The failed treatment count for weight 54kg is: 358
        The percentage efficacy for weight 54kg is: 64.2
        
    - [ ]  `--AL`
        
        Today's date:  2025-04-15
        Length of simulation:  0:00:31.304186
        The failed treatment count for weight 54kg is: 31
        The percentage efficacy for weight 54kg is: 96.9
        
- [ ]  Made the following changes:
    - [ ]  Edited `output_results`  to `output_art_monotherapy` and `output_combination_AL`
    - [ ]  Changed `G_OUTPUT_TYPE` to `G_CLO_OUTPUT_TYPE`
    - [ ]  Edited the output function calls in `--art` and `--AL` modules
    - [ ]  Checking efficacies for weight = 54, n = 1k
        - [ ]  `--art` 
        
        Today's date:  2025-04-15
        Length of simulation:  0:00:25.246640
        The failed treatment count for weight 54kg is: 986
        The percentage efficacy for weight 54kg is: 1.4000000000000057
        - [ ]  `--AL`, `--pmax_lum = 0`
            
            Today's date:  2025-04-15
            Length of simulation:  0:00:32.419496
            The failed treatment count for weight 54kg is: 326
            The percentage efficacy for weight 54kg is: 67.4
            
        - [ ]  `--AL`
            
            Today's date:  2025-04-15
            Length of simulation:  0:00:31.689173
            The failed treatment count for weight 54kg is: 20
            The percentage efficacy for weight 54kg is: 98.0
            
    - [ ]  Changed `--pmax_art` = 0.9997 only for `--art`
    
               Today's date:  2025-04-15 
                Length of simulation:  0:00:25.881040 
               The failed treatment count for weight 54kg is: 0
                The percentage efficacy for weight 54kg is: 100.0
    
- [ ]  Changing the formula to calculate the PMF in `main.cpp`
    - [ ]  from 
    Line 337 `double stepsize_PMF = pow( G_CLO_PMF, 1.0 / (48.0/maximum_enforced_stepsize) );`
        
        
        379`double dd_PMF = stepsize_PMF * ( 1.0 / ( 1.0 + ( dyn1->y0[ dyn1->dim - 1 ] / G_DENSITY_50 ) ) );`
        380`dyn1->y0[ dyn1->dim - 1 ] *= dd_PMF;`
        381 `dyn2->y0[ dyn2->dim - 1 ] *= dd_PMF;` 
        
        to 
        
        `if (true)`
        
        `{`
        
        `double dd_factor = 1.0 / ( 1.0 + ( dyn1->y0[ dyn1->dim - 1 ] / G_DENSITY_50 ) );`
        
        `double dd_48hr_PMF = G_CLO_PMF * dd_factor;`
        
        `double dd_stepsize_PMF = pow( dd_48hr_PMF, 1.0 / (48.0/maximum_enforced_stepsize) );`
        
        `dyn1->y0[ dyn1->dim - 1 ] *= dd_stepsize_PMF;`
        
        `dyn2->y0[ dyn2->dim - 1 ] *= dd_stepsize_PMF;`
        
        `}`
        
        - Results
            - `--AL --pmax_lum 0`
                
                Today's date:  2025-04-15
                Length of simulation:  0:00:30.086937
                The failed treatment count for weight 54kg is: 981
                The percentage efficacy for weight 54kg is: 1.9000000000000057
                
            - `--AL`
                
                Today's date:  2025-04-15 
                Length of simulation:  0:00:28.494704 
                The failed treatment count for weight 54kg is: 77
                The percentage efficacy for weight 54kg is: 92.3
                
    
    - [ ]  Editing artemether schedule dosing schedule in  `pkpd_dha.cpp`
        - From:
            
            298 `void pkpd_dha::generate_recommended_dosing_schedule() {`
            
            `// TODO NEEDS TO BE DONE BY AGE AND WEIGHT`
            
            `// one tablet is 40mg of dihydroartemisinin`
            
            `// TODO: need to get tablet schedule by weight, age, pregnancy status`
            
            `double num_tablets_per_dose;`
            
            `if( weight < 5.0 )`
            
            `{ num_tablets_per_dose = 0.0;`
            
            `} else if( weight < 8.0 )`
            
            `{ num_tablets_per_dose = 0.5;`
            
            `} else if( weight < 11.0 ) {`
            
            `num_tablets_per_dose = 0.75; }`
            
            `else if( weight < 17.0 ) {`
            
            `num_tablets_per_dose = 1.0;` 
            
            `} else if( weight < 25.0 ) {`
            
            `num_tablets_per_dose = 1.5;`
            
            `} else if( weight < 36.0 ) {` 
            
            `num_tablets_per_dose = 2.0;` 
            
            `} else if( weight < 60.0 ) {` 
            
            `num_tablets_per_dose = 3.0; }` 
            
            `else if( weight < 80.0 ) {` 
            
            `num_tablets_per_dose = 4.0; }` 
            
            `else { num_tablets_per_dose = 5.0; }` 
            
            `// TODO REMOVE PLACEHOLDER BELOW`
            
            `num_tablets_per_dose = 1.0;`
            
            `double total_mg_dose = num_tablets_per_dose * 40.0; // one tablet is 40 mg of dihydroartemisinin`
            
            `v_dosing_times.insert( v_dosing_times.begin(), 3, 0.0 );`
            
            `v_dosing_times[1] = 24.0;`
            
            `v_dosing_times[2] = 48.0;`
            
            `v_dosing_amounts.insert( v_dosing_amounts.begin(), 3, total_mg_dose );`
            
        
        355 `}`
        
        To: 
        
        `void pkpd_dha::generate_recommended_dosing_schedule()`
        
        `{  // DOSING GUIDELINES SAY    0,  8, 24, 36, 48, 60`
        
           `// BUT WE CAN JUST DO       0, 12, 24, 36, 48, 60`
        
        `double num_tablets_per_dose;`
        
        `if( weight < 15.0 )`
        
        `{  num_tablets_per_dose = 0.50;`
        
        `} else if( weight < 25.0 )`
        
        `{ num_tablets_per_dose = 1.0;`
        
        `} else if( weight < 35.0 )`
        
        `{ num_tablets_per_dose = 1.50;`
        
        `} else`
        
        `{ num_tablets_per_dose = 2.0;`
        
        `}`
        
        `// Artemether given by weight, twice daily, for a total of three days - WHO guidelines, 2024`
        
        `double total_mg_dose_per_occassion = num_tablets_per_dose * 40.0;`
        
        `v_dosing_times.insert( v_dosing_times.begin(), 6, 0.0 );`
        
        `v_dosing_times[0] = 0.0;`
        
        `v_dosing_times[1] = 12.0;`
        
        `v_dosing_times[2] = 24.0;`
        
        `v_dosing_times[3] = 36.0;`
        
        `v_dosing_times[4] = 48.0;`
        
        `v_dosing_times[5] = 60.0;`
        
        `v_dosing_amounts.insert( v_dosing_amounts.begin(), 6, total_mg_dose_per_occassion );`
        
        `}`
        
    
- [ ]  Editing lumefantrine dosing schedule in pkpd_lum.cpp
    - From
    
       320  `void pkpd_lum::generate_recommended_dosing_schedule()`
    
         `{ // DOSING GUIDELINES SAY    0,  8, 24, 36, 48, 60`
    
    `// BUT WE CAN JUST DO       0, 12, 24, 36, 48, 60`
    
    `double num_tablets_per_dose;`
    
    `if( weight < 5.0 )`
    
    `{ num_tablets_per_dose = 0.0;`
    
    `} else if( weight < 14.0 )` 
    
    `{ num_tablets_per_dose = 1.0;`
    
    `} else if( weight < 24.0 )`
    
    `{ num_tablets_per_dose = 2.0;`
    
    `} else if( weight < 34.0 )`
    
    `{ num_tablets_per_dose = 3.0;`
    
    `} else`
    
    `{ num_tablets_per_dose = 4.0;`
    
    `} // NOTE - do not confuse this with daily dosing - AL is taken twice daily, two occassions per day`
    
         `total_mg_dose_per_occassion = num_tablets_per_dose * 120.0;`
    
    `v_dosing_times.insert( v_dosing_times.begin(), 6, 0.0 );`
    
    `v_dosing_times[0] = 0.0;`
    
    `v_dosing_times[1] = 12.0;`
    
    `v_dosing_times[2] = 24.0;`
    
    `v_dosing_times[3] = 36.0;`
    
    `v_dosing_times[4] = 48.0;`
    
    `v_dosing_times[5] = 60.0;`
    
    `v_dosing_amounts.insert( v_dosing_amounts.begin(), 6, total_mg_dose_per_occassion );`  
    
        363 `}`
    
    - To
        
        `void pkpd_lum::generate_recommended_dosing_schedule()`
        
        `{ // DOSING GUIDELINES SAY    0,  8, 24, 36, 48, 60`
        
        `// BUT WE CAN JUST DO       0, 12, 24, 36, 48, 60`
        
        `double num_tablets_per_dose;` 
        
        `if( weight < 15.0 )`
        
        `{ num_tablets_per_dose = 1.0;`
        
        `} else if( weight < 25.0 )`
        
        `{ num_tablets_per_dose = 2.0;`
        
        `} else if( weight < 35.0 )`
        
        `{ num_tablets_per_dose = 3.0;`
        
        `} else`
        
        `{ num_tablets_per_dose = 4.0;`
        
        `} // NOTE - do not confuse this with daily dosing - AL is taken twice daily, two occassions per day`
        
        `total_mg_dose_per_occassion = num_tablets_per_dose * 120.0;
        v_dosing_times.insert( v_dosing_times.begin(), 6, 0.0 );
        v_dosing_times[0] = 0.0;
        v_dosing_times[1] = 12.0;
        v_dosing_times[2] = 24.0;
        v_dosing_times[3] = 36.0;
        v_dosing_times[4] = 48.0;
        v_dosing_times[5] = 60.0;
        v_dosing_amounts.insert( v_dosing_amounts.begin(), 6, total_mg_dose_per_occassion );`
        
        `}`
        

- [ ]  Implementing age and weight as:
    
    `dyn1->age = G_CLO_AGE;
    dyn1->weight = G_CLO_WEIGHT;
    dyn1-> patient_blood_volume = 5500000.0 * (dyn1-> weight/dyn1-> median_weight);` 
    

`dyn2->age = G_CLO_AGE;
dyn2->weight = G_CLO_WEIGHT;
dyn2-> patient_blood_volume = 5500000.0 * (dyn2-> weight/dyn2-> median_weight);       // 5.5L of blood for an adult individual`

### Observation: The median_weight for artemether is 48.5 and for lumefantrine is 54 (42.0 in initialize params)

- [ ]  `--AL --pmax_lum 0`

`Today's date:  2025-04-16
Length of simulation:  0:00:23.715065
The failed treatment count for weight 54kg is: 1000
The percentage efficacy for weight 54kg is: 0.0`

- [ ]  `--AL`

`Today's date:  2025-04-16
Length of simulation:  0:00:24.085654
The failed treatment count for weight 54kg is: 383
The percentage efficacy for weight 54kg is: 61.7`

- [ ]  Scaling Hill equation by `patient_blood_volume` 
Note: EC50 is also scaled by `patient_blood_volume`
    - [ ]  In `pkpd_dha.cpp`
        
        From 
        
        106 `double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow(y[8],p->pdparam_n) / (pow(y[8],p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)) );`
        
        To 
        
        `double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow((y[8]/p -> patient_blood_volume),p->pdparam_n) / (pow((y[8]/p -> patient_blood_volume),p->pdparam_n) + pow((p->pdparam_EC50/p -> patient_blood_volume),p->pdparam_n)));`
        
    - [ ]  In `pkpd_lum.cpp` 
    From
    
    104 `double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow(y[1],p->pdparam_n) / (pow(y[1],p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)) );`
    
    To 
    
    `double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow((y[1]/p -> patient_blood_volume),p->pdparam_n) / (pow((y[1]/p -> patient_blood_volume),p->pdparam_n) + pow((p->pdparam_EC50/p -> patient_blood_volume),p->pdparam_n)) );`
    
- [ ]  `--AL --pmax_lum 0`
    
    `Today's date:  2025-04-16
    Length of simulation:  0:00:24.992387
    The failed treatment count for weight 54kg is: 1000
    The percentage efficacy for weight 54kg is: 0.0`
    
- [ ]  `--AL`
    
    `Today's date:  2025-04-16
    Length of simulation:  0:00:24.718199
    The failed treatment count for weight 54kg is: 434
    The percentage efficacy for weight 54kg is: 56.6`
    
- [ ]  Making median_weight the same in `pkpd_artemether.cpp` (48.5) and `pkpd_lum.cpp` (54) to 54 Kg
    - [ ]  Editing age and weight implementation
    
    `dyn1->age = G_CLO_AGE;
    dyn1->patient_weight = G_CLO_WEIGHT;
    dyn1->weight = dyn1->patient_weight;
    dyn1-> patient_blood_volume = 5500000.0 * (dyn1-> weight/dyn1-> median_weight);`
    
    Weight not getting set properly as dose is generated within constructor, implementing `intialize()` which calls `generate_recommended_dosing_schedule()` and `initialize_params()`
    
    `dyn2->age = G_CLO_AGE;
    dyn2->patient_weight = G_CLO_WEIGHT;
    dyn2->weight = dyn2->patient_weight;
    dyn2-> patient_blood_volume = 5500000.0 * (dyn2-> weight/dyn2-> median_weight);`
    
    - [ ]  `--AL --pmax_lum 0`
        
        `Today's date:  2025-04-16
        Length of simulation:  0:00:25.210110
        The failed treatment count for weight 54kg is: 1000
        The percentage efficacy for weight 54kg is: 0.0`
        
    - [ ]  `--AL`
        
        `Today's date:  2025-04-16
        Length of simulation:  0:00:24.529019
        The failed treatment count for weight 54kg is: 424
        The percentage efficacy for weight 54kg is: 57.6`
        
    
- [ ]  Moving `intial_log10_parasitemia` from constructor to `initialize_params()` in `pkpd_dha.cpp`
    - [ ]  `--AL --pmax_lum 0`
        
        `Today's date:  2025-04-16
        Length of simulation:  0:00:24.713192
        The failed treatment count for weight 54kg is: 1000
        The percentage efficacy for weight 54kg is: 0.0`
        
    - [ ]  `--AL`
        
        `Today's date:  2025-04-16
        Length of simulation:  0:00:24.204257
        The failed treatment count for weight 54kg is: 407
        The percentage efficacy for weight 54kg is: 59.300000000000004`
        
    
    - [ ]  Scaling Hill equation by `patient_blood_volume` without scaling EC50 is  by `patient_blood_volume`
    
    - [ ]  `pkpd_dha.cpp`
        
        From 
        
        106 `double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow((y[8]/p -> patient_blood_volume),p->pdparam_n) / (pow((y[8]/p -> patient_blood_volume),p->pdparam_n) + pow((p->pdparam_EC50/p -> patient_blood_volume),p->pdparam_n)));`
        
        To
        
        `double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow((y[8]/p -> patient_blood_volume),p->pdparam_n) / (pow((y[8]/p -> patient_blood_volume),p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)));`
        
    - [ ]  `pkpd_lum.cpp`
        
        From 
        
        `double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow((y[1]/p -> patient_blood_volume),p->pdparam_n) / (pow((y[1]/p -> patient_blood_volume),p->pdparam_n) + pow((p->pdparam_EC50/p -> patient_blood_volume),p->pdparam_n)) );` 
        
        To 
        
        `double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow((y[1]/p -> patient_blood_volume),p->pdparam_n) / (pow((y[1]/p -> patient_blood_volume),p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)) );` 
        
    
    - [ ]  `--AL --pmax_lum 0`
        
        
        `Today's date:  2025-04-16
        Length of simulation:  0:00:25.027019
        The failed treatment count for weight 54kg is: 1000
        The percentage efficacy for weight 54kg is: 0.0`
        
    - [ ]  `--AL`
        
        `Today's date:  2025-04-16
        Length of simulation:  0:00:23.549212
        The failed treatment count for weight 54kg is: 1000
        The percentage efficacy for weight 54kg is: 0.0`
        
    
    ### Observation: 
    a) EC50/patient_blood_volume seems to be the right call
    
    ### b) For some reason the parasitaemia, regardless of monotherapy/combination or weight is at 450093.809 at hour 671
    → Reverting to previous Hill Equation formula
    

- [ ]  Copying class and functions in `pkpd_dha.cpp` to `pkpd_artemether.cpp`
- Renamed `pkpd_dha.cpp` to `pkpd_artemether.cpp`
- Renamed pkpd_dha.h to pkpd_artemether.h
- Renamed pkpd_dha_og.cpp to pkpd_dha.cpp
- Renamed pkpd_dha_og.cpp to pkpd_dha.h
    - Also replaced  all occurences of ‘dha’  in `pkpd_artemether.cpp` to ‘artemether’
- Changed command-line flags
    - `--pmax_art` to `--pmax_dha` which is stored in the variable `G_CLO_PMAX_DHA`
    - `--hill_art` to `--hill_dha`  which is stored in the variable `G_CLO_HILL_COEFF_DHA`
    - `--ec50_art` to `--ec50_dha` which is stored in the variable `G_CLO_EC50_DHA`
- Changed function `output_results_monotherapy_art` to `output_results_monotherapy_dha`
- Created `G_CLO_HILL_COEFF_ARTEMETHER` , `G_CLO_EC50_ARTEMETHER` and `G_CLO_PMAX_ARTEMETHER`
- Changed the following in `main.cpp`
    
    `dyn1->pdparam_n = G_CLO_HILL_COEFF_DHA;
    dyn1->pdparam_EC50 = G_CLO_EC50_DHA;
    dyn1->pdparam_Pmax = G_CLO_PMAX_DHA;`
    
    to 
    
    `dyn1->pdparam_n = G_CLO_HILL_COEFF_ARTEMETHER;
    dyn1->pdparam_EC50 = G_CLO_EC50_ARTEMETHER;
    dyn1->pdparam_Pmax = G_CLO_PMAX_ARTEMETHER;`
    
- Added `pkpd_rtemether.cpp` and `pkpd_artemether.h` to CMakeLists.txt
    
    `set(SRC_FILES main.cpp pkpd_adq.cpp pkpd_adq.h pkpd_dha.cpp pkpd_dha.h pkpd_lum.cpp pkpd_lum.h pkpd_ppq.cpp pkpd_ppq.h pkpd_artemether.cpp pkpd_artemether.h)`
    

- [ ]  Checking efficacy
    - [ ]  `--AL --pmax_lum`
        
        `Today's date:  2025-04-16
        Length of simulation:  0:00:23.717361
        The failed treatment count for weight 54kg is: 1000
        The percentage efficacy for weight 54kg is: 0.0`
        
    - [ ]  `--AL`
        
        `Today's date:  2025-04-16 
        Length of simulation:  0:00:24.317273 
        The failed treatment count for weight 54kg is: 424
        The percentage efficacy for weight 54kg is: 57.6`
        
        `Today's date:  2025-04-16
        Length of simulation:  0:00:24.554705
        The failed treatment count for weight 54kg is: 395
        The percentage efficacy for weight 54kg is: 60.5`
        
        `Today's date:  2025-04-16
        Length of simulation:  0:00:25.475640
        The failed treatment count for weight 54kg is: 393
        The percentage efficacy for weight 54kg is: 60.699999999999996`
        

### Wednesday 04/17 and Saturday 04/19: Calibrating AL model to weight = 54 Kg (Wednesday),  re-calibrated model after fixing bug (Saturday)

- Forgot if I calibrated the model weight = 54 using concentration or total units in blood
- Re-calibrated the model after editing

### Results:

![pkpd_AL_pmax_art_test_01_efficacy.png](pkpd_AL_pmax_art_test_01_efficacy.png)

![pkpd_AL_ec50_lum_test_10k_efficacy.png](pkpd_AL_ec50_lum_test_10k_efficacy.png)

![pkpd_AL_weight_calibration_01.png](pkpd_AL_weight_calibration_01.png)