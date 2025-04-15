# pkpd_artemether:

## Required files and corresponding descriptions:

### main.cpp: 
	Contains the main() function, where the program runs and interacts with the class.

### pkpd_artemether.cpp: 
	Contains the implementation of the functions declared in the class (pkpd_artemether).

### pkpd_artemether.h: 
	Declares the class interface (functions and variables) without implementation details.

<br>

## pkpd_artemether.h structure/control of flow:


1. Include guards 

	```` #ifndef PKPD_ARTEMETHER ```` : 
	
	```` #ifndef ```` means "if not defined." It checks whether ```` PKPD_ARTEMETHER ```` has already been defined. If it has NOT been defined, the compiler    proceeds to the next line.

	```` #define PKPD_ARTEMETHER ```` : 
	
	If ```` PKPD_ARTEMETHER ```` was not defined before, this line defines it so that the next time this file is included, the compiler knows it has already been processed.

2. Load required libraries

````
	#include <vector> 
	#include <math.h>
	#include <gsl/gsl_errno.h>
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_odeiv.h>
	#include <gsl/gsl_vector.h>
	#include <gsl/gsl_rng.h> // random number generators from Gnu Scientific Library
	#include <gsl/gsl_cdf.h>
	#include <gsl/gsl_randist.h>  // this includes the gsl_ran_gaussian function
````

3. Include global variables that can be used outside the class: ```` enum parameter_index_artemether ````
	- ```` enum parameter_index_artemether ```` contains named indicies (0 to length of enum)
	- Names of the indicies are: ````i_artemether_KTR, i_artemether_k20, i_artemether_F1_indiv, i_artemether_F1_thisdose, artemether_num_params````

## ```` class pkpd_artemether ```` structure:

### Declare public methods:

1. Class Constructor: ```` pkpd_artemether() ````
2. Class Destructor: ```` ~pkpd_artemether() ````
3. Boolean variable to turn stochasticity on or off:  ```` stochastic ````
4. Vector  of parameters: ```` vprms ````

### Declare ODE params, functions and methods:

2. Variable representing the dimensionality of the model: ```` dim ````
3. GSL (Gnu Scientific Library) function calls to 	````step````, ````control```` and ````evolve````
	
	Note: The GSL ODE functions (```` gsl_odeiv_step, gsl_odeiv_control, and gsl_odeiv_evolve````) are used to control how the integration happens: step size, accuracy, and the evolution of the system over time.

4. Function to integrate the ODEs: ```` predict() ````
5. Variable to store and update system conditions at t=0: ```` y0 ````
6. Jacobian function: ```` jac() ````
7.  Function to calculate the RHS of the ODE: ```` rhs_ode() ````

	Note: During the integration process (in your predict function), the solver takes the current state of the system (the values of the variables), calculates the right-hand side of the ODE (using ```` rhs_ode ````), and updates the system’s state over time.

### Declare PK parameters and functions:

8. Function to initalize parameters and generate dosing schedule: ```` initialize() ````
9. Function to initialize parameters: ```` initialize_params() ````
10. Function to redraw new parameters before administering new dose: ```` redraw_params_before_newdose() ````

### Declare PD parameters and functions:

11. Function to set the parasitemia at t =0: ```` set_parasitemia() ````
12. Variables representing the Hill co-efficient, EC50 for artemether and PMax for artemether: ```` pdparam_n, pdparam_EC50 ```` and ```` pdparam_Pmax ```` respectively.

### Declare parameters and functions related to dosing schedule:

13. Function to generate dosage times, frequency and concentration of dose: ```` generate_recommended_dosing_schedule() ````
14. Vector containing dosing times for artemether: ```` v_dosing_times ````
15. Vector containing the doses for artemether in mg: ```` v_dosing_amounts ````
16. Variable to represent the number of doses already given: ```` num_doses_given ````
17. Boolean variable to represent if additional doses still remain to be given: ```` doses_still_remain_to_be_taken ````
18. Function that checks whether its time to give a dose or not and returns a Boolean: ```` we_are_past_a_dosing_time() ````
19. Function to administer the next dose to patient: ```` give_next_dose_to_patient() ````

### Declare variables and functions related to Patient Characteristics:

20. ```` patient_weight ```` - in kg.
21. ```` median_weight ```` - Median weight for which estimates are calculated for
22. ```` weight ```` - Just another variable to store the actual weight that is gonna be used
23. ```` age ````
24. ```` patient_blood_volume```` : In microliters, so 250,000 mul for infants and ~6,000,000 for large adults
25. boolean variable to represent preganancy status: ```` pregnant ````
26. boolean variable to represent sex: ```` is_male ````

### Storage variables for the dynamics of the PK/PD curves:

27. Vector of the hourly drug concentrations in the blood: ```` v_concentration_in_blood ````
    - Should be mg/ml (check)
28. Vector of hours over 28 days (672 hours): ```` v_concentration_in_blood_hourtimes ````
29. Vector storing parasite density in the blood after drug action: ```` v_parasitedensity_in_blood ````
30. Variable to keep track of the time (?): ```` num_hours_logged ````
31. Pointer (rng) to a GSL random number generator (```` gsl_rng ````), which you'll later initialize, use to generate random numbers, and eventually free once you're done: ```` gsl_rng *rng ````
32. Close class

33. Close include guards

	```` #endif // PKPD_ARTEMETHER ````: 
	
	This closes the conditional ```` #ifndef ```` block. If the header file gets included multiple times, the first inclusion defines ````PKPD_ARTEMETHER```` , and subsequent inclusions skip everything inside the guard.
