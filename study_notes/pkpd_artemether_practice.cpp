// Step 0: Include required headers and libraries

#include "assert.h" 

/* Includes the assert.h header, which provides the assert macro in C++.
What assert.h Does: 
    The assert macro is used for debugging by checking assumptions in your code.
    If the condition inside assert(condition) evaluates to false, the program prints an error message and terminates. 
    
When to Use assert:
    Use it to catch programming errors early in development.
    Typically used for checking preconditions in functions.
    It is removed in release builds if NDEBUG is defined (#define NDEBUG).*/


#include "pkpd_artemether_practice.h" // include/connect to the header file for pkpd_artemether class

#include <vector>


/* ************************************************************************* */

bool pkpd_artemether::stochastic = true; 

// setting variable stochastic belonging to class pkpd_artemether to true
// to control for stochasticity in the model

/* ************************************************************************* */

// Step 1: Define the constructor for object of class pkpd_artemether
// time to define the constructor
// which means, we are now going to specify what happens when a pkpd_artemether object is created

// pkpd_artemether::pkpd_artemether() // Old constructor, error-prone

pkpd_artemether::pkpd_artemether(double patient_age, double patient_weight) {

// Step 2a: Create a vector of doubles to hold the parameters
    // vprms is a vector of doubles that we declared in the header
    // the code below creates a vector of doubles (0.0) of size artemether_num_params
        // artemether_num_params is always the last parameter in the enum parameter_index_artemether
            // which always gives the total number of parameters in the index
            // right now, it is 5 (remember that the last one is not a parameter, but just the number of parameters)

    vprms.resize(artemether_num_params, 0.0); // Resizing the vector to hold artemether_num_params elements

    // assert checks if the size of vprms is equal to artemether_num_params
    // if it is not, the program will terminate and print an error message
    assert( vprms.size()== artemether_num_params ); 

// Step 2b: Create a dynamically allocated array of doubles, y0, to hold the state variables
    // Number of compartments = 9, 10th compartment is for the parasite density
    dim = 10;

    // Creating a pointer to the first element of a dynamically allocated array of doubles

    y0 = new double[dim]();
    //for(int i=0; i<dim; i++) y0[i]=0.0; // initializing the array to zero, not required as the () above automatically initializes it to zero

     // above y0 is initialized to zero, because it's the dosing schedule that creates
     // non-zero initial conditions for the ODEs

    /*  
        Every time new is invoked:
        - A block of memory is allocated on the heap.
        - The variable (y0) receives the address of the first element of the allocated block, and we access the rest of the elements using pointer arithmetic.
        - C++ does not store the array size along with the pointer!
        - The pointer y0 only knows where the first element is, not how many elements follow it.
        - The responsibility of keeping track of dim (the size) falls on the programmer.
        - This is different from a statically allocated array, where the compiler knows the size.
        - If you didn’t initialize it: You might get garbage values (random stuff in memory)
            Or worse: subtle bugs where certain compartments start with huge or NaN values, and your model spirals into chaos */

    // NOTE: REPLACING MANUAL ALLOCATION WITH std::vector<double> y0(dim) AS IT IS MORE EFFICIENT
    // std::vector is a dynamic array that can grow and shrink in size
    // It automatically manages memory for you, so you don't have to worry about allocating and deallocating memory manually
    // It also provides bounds checking, so you can't access elements outside the array size
    // Also comes with a lot of useful member functions like push_back, pop_back, etc.

    //std::vector<double> y0(dim); 
        // No need to initialize to zero, as std::vector does that automatically for type double

    // Set default initial parasitaemia (can be overridden later using set_parasitaemia())
    // the last differential equation is for the parasitaemia; it is a per/ul measure
    y0[dim-1] = 10000.0; // BUT we do need to set the initial parasitaemia at time zero to something positive
                         // this should be obtained from the person class

// Step 3: Initialize the ODE solver

    /* GSL ODE Solver Setup
        - GSL (GNU Scientific Library) is a C library for numerical computing.
        - The ODE solver is used to solve ordinary differential equations (ODEs).
        - The GSL ODE solver requires 
            - a function to define the system of ODEs, 
            - a step size
            - and an initial condition.
    
    1. Method selection: 
        We first set the method we want to use to solve the ODE — in this case, the Runge-Kutta-Fehlberg 4(5) method. 
        This is stored in GSL as the built-in constant: gsl_odeiv_step_rkf45.
    */

    const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45; 

    /* 2. os — Stepper object: We create a stepper object (os) using:
        os = gsl_odeiv_step_alloc(method, dimension) 
        This tells GSL: Use the RKF45 method to solve a system with dim differential equations.
        The stepper handles taking individual numerical steps in the solution process.*/

    os 	= gsl_odeiv_step_alloc(T, dim);

    /* 3. oc — Controller object: We create a controller object to monitor and limit the error for each step. 
        oc = gsl_odeiv_control_y_new(abs_err, rel_err)
        This controls the absolute and relative error tolerance — i.e., how far off we’re allowed to be at each integration step. */
    
    oc 	= gsl_odeiv_control_y_new (1e-6, 0.0);

    /* 4. oe — Evolver object/Evolution system
        - The evolver object is responsible for the actual evolution of the system over time. 
        - It uses the stepper and controller to do this. 
        - This manages the overall integration process — it:
            * Calls the stepper and controller
            * Tracks time
            * Adjusts step sizes
            * Stores internal state */

    oe 	= gsl_odeiv_evolve_alloc(dim);
    
// Step 4: Default patient characteristics 

    weight_initialized = false; // To check if the member variable 'weight' has been initialized
    age_initialized = false;    // To check if the member variable 'age' has been initialized
    patient_blood_volume = 5500000.0; // 5.5L of blood for an adult individual // 5500000 ul gives 5.5 L
    median_weight  = 48.5; // In KG, check Hoglund, Kloprogge papers
    set_age_and_weight(patient_age, patient_weight);    

    num_doses_given = 0;
    num_hours_logged = 0;

    is_male = false;
    is_pregnant = false; // To check if female patient is pregnant

    doses_still_remain_to_be_taken = true;

// Step 5: PD Parameters

    pdparam_n = 20.0; // Hill coefficient (steepness of the dose-response curve)
    pdparam_EC50 = 0.1; //  Drug concentration at which 50% of max effect is seen
    pdparam_Pmax = 0.99997; // max daily killing rate; converted to hourly later
    

    /* rng likely refers to a random number generator pointer
    Setting it to NULL here means:
        “No RNG assigned yet — we’ll probably initialize it later when needed”
    This is common when:
        Some parts of the model (e.g., stochastic dosing, variability in response) use randomness
        But you don’t want to initialize the RNG until it’s actually needed, or until the seed is specified*/
        
    rng = NULL;

} // end of constructor

/* **************************************************************************************************************** */

// Step 6: Define the destructor for object of class pkpd_artemether

/* - It is defined using the same name as the class, preceded by the tilde (~) symbol
   - While defining a destructor, we specify what is deleted/destroyed when an object is deleted or destroyed
   - The destructor is called automatically when the object goes out of scope or is deleted
   - It cleans up the resources allocated in the constructor
   - It frees the memory allocated for the ODE solver objects

    - It is important to free the memory allocated for the ODE solver objects to avoid memory leaks
        - Memory leaks occur when memory is allocated but not freed, leading to increased memory usage over time
        - In this case, we are using the GSL library to solve ODEs, and we need to free the memory allocated for the ODE solver objects

    - It also frees the memory allocated for the y0 array.
        - In this case, we are using a vector, so we don't need to delete it manually
        - The vector will automatically free its memory when it goes out of scope
    */

pkpd_artemether::~pkpd_artemether()
{
    if (y0 != nullptr) {
        delete[] y0;
        y0 = nullptr;
    }

    if(oe) gsl_odeiv_evolve_free(oe);

/*  This is freeing/deallocating memory associated with the evolution object (oe), which was created using:
        - oe = gsl_odeiv_evolve_alloc(dim)
        - The evolution object is part of the GSL ODE solver framework. 
        - It manages the integration process — that is, it keeps track of the 
            - current step size
            - the current time t, 
            - and internally stores helper data needed to carry out the stepping method chosen (like RKF45).
    
    - Why free it?
        - Because it was created using a GSL allocation function, you must manually free it when you're done, to avoid a memory leak.*/

    if (oc) gsl_odeiv_control_free(oc);

/*  - This line deallocates memory used by the controller object oc, which you created earlier like this:
        oc = gsl_odeiv_control_y_new(...);
    - The controller (oc) is what keeps the error of the numerical integration within bounds. 
    - It checks the difference between expected and calculated values and tells the stepper whether to shrink or grow the step size.
    - By calling gsl_odeiv_control_free(oc);, you're telling GSL:
        - “We’re done controlling errors for this system — release the memory!” */
    
    if (os) gsl_odeiv_step_free(os);

/* This frees the memory used by the stepper object os, which you originally created with:
        - os = gsl_odeiv_step_alloc(...); 
    The stepper is the numerical method (e.g., RKF45) that actually computes the solution for each integration step. 
    Once the object is no longer needed — for instance, when the simulation ends — this call safely deallocates the associated memory. */
}

// Summary: evolve, control and step objects have an allocation function and a free function
    // The free function is used to deallocate the memory allocated by the allocation function
    // This is done by the destructor

/* ************************************************************************************************************** */

// Step 7: Time to define some functions!

// Function 1: set_age_and_weight

void pkpd_artemether::set_age_and_weight( double a, double w )
{
    age = a;
    weight = w;
    patient_blood_volume = 5500000.0 * (w/median_weight);
    weight_initialized = true; // this is true if the weight has been set
    age_initialized = true; // this is true if the age has been set
    
}

// Function 2: set_parasitaemia

void pkpd_artemether::set_parasitaemia( double parasites_per_ul )
{
    y0[dim-1] = parasites_per_ul; //  the final ODE equation is always the Pf asexual parasitaemia
}

// Function 3: generate_recommended_dosing_schedule

void pkpd_artemether::generate_recommended_dosing_schedule()
{
    // TODO NEEDS TO BE DONE BY AGE AND WEIGHT
    
    // TODO: need to get tablet schedule by weight, age, pregnancy status

    if (!weight_initialized)
    {
        fprintf(stderr, "\n\tERROR: Weight not initialized. Exiting program.\n");
        return;
    }

    num_tablets_per_dose = 1.0; // one tablet is 40mg of artemether
    
    if( weight < 15.0 )
    {
        num_tablets_per_dose = 0.50; // 0.50 * 40.0 = 20 mg
    }
    else if( weight < 25.0 )
    {
        num_tablets_per_dose = 1.00; // 1.00 * 40.0 = 40 mg
    }
    else if( weight < 35.0 )
    {
        num_tablets_per_dose = 1.50; // 1.50 * 40.0 = 60 mg
    }
    else
    {
        num_tablets_per_dose = 2.00; // 2.00 * 40.0 = 80 mg
    }
    
    total_mg_dose = num_tablets_per_dose * 40.0; // one tablet is 40 mg of artemether

    v_dosing_times.insert( v_dosing_times.begin(), 6, 0.0 );
    v_dosing_times[0] = 0.0;
    v_dosing_times[1] = 12.0;
    v_dosing_times[2] = 24.0; 
    v_dosing_times[3] = 36.0;
    v_dosing_times[4] = 48.0;
    v_dosing_times[5] = 60.0;
    
    v_dosing_amounts.insert( v_dosing_amounts.begin(), 6, total_mg_dose );
}

// Function 4: we_are_past_a_dosing_time 

bool pkpd_artemether::we_are_past_a_dosing_time( double current_time )
{
    // check if there are still any doses left to give
    if( num_doses_given < v_dosing_times.size() )                   // num_doses_given is set to 0 in the constructor
                                                                    // v.dosing_times.size() is the length of the vector; here it is 6
    {
        if( current_time >= v_dosing_times[num_doses_given] )
        {
            return true;
        }
    }

    return false;
}

// Function 5: redraw_params_before_newdose

void pkpd_artemether::redraw_params_before_newdose()
{
    // TODO the function below is a copy-and-paste from the PPQ function; must be modified

    /*  Called before each dose
            Applies inter-dose variability (e.g. by redrawing ETA for KTR)
            Can adjust parameters like KTR or F1 between doses, for the same patient */

    double OCC = 1.0 + static_cast<double>(num_doses_given); // NOTE the RHS here is a class member
                                                             // static_cast<double>(num_doses_given) is used to convert the integer num_doses_given to a double
                                                             // changed (double)num_doses_given from OG code to static_cast<double>(num_doses_given)
    
    if(pkpd_artemether::stochastic) 
    {
        double ETA_rv = gsl_ran_gaussian( rng, sqrt(0.23000) );     // this is ETA6, ETA7, and ETA8
        vprms[i_artemether_KTR] *= exp(-ETA_rv);                    // WARNING this behavior is strange ... check if this is the right way to do it
                                                                    // Check if KTR = 8/MTR is set anywhere else
    }
    
}


// Function 6: give_next_dose_to_patient

void pkpd_artemether::give_next_dose_to_patient( double fractional_dose_taken ) // Called in main.cpp, fractional_dose_taken is set to 1.0 in main.cpp
{
    if( doses_still_remain_to_be_taken ) // doses_still_remain_to_be_taken is set to true in the constructor
    {
        redraw_params_before_newdose(); // these are the dose-specific parameters that you're drawing here
        
        y0[0] +=  v_dosing_amounts[num_doses_given] * fractional_dose_taken;
        
        num_doses_given++;

        if( num_doses_given >= v_dosing_amounts.size() ) doses_still_remain_to_be_taken=false;

    }

}


// Function 7: intialize_parameters

void pkpd_artemether::initialize_parameters() {
    //WARNING - THE AGE MEMBER VARIABLE MUST BE SET BEFORE YOU CALL THIS FUNCTION

    /* Summary:
    1. KTR (the absorption rate constant) varies between patients — meaning each simulated patient can have a different value.
        But it does not vary between doses for the same patient — so once it’s set here, it stays constant across all 6 doses.
    2. You apply ETA-based variability to KTR only in redraw_params_before_newdose() — that’s where inter-dose variability (IOV) happens.
        But here in initialize_params(), you’re just setting a fixed baseline KTR using the MT value (mean transit time), 
        which stays the same across all doses for a given patient.
    3. Called once per patient, during setup.
    4. Sets baseline PK values (like KTR, CL, V2)
    5. Applies inter-patient variability only → i.e., every patient gets their own fixed parameters*/

    double mw=48.5; // median weight used for CL/V2 scaling below

    vprms[i_artemether_F1_thisdose] = 1.0;


    // The log10 total parasite load in the body at t = 0
    /*  - This is the initial parasite density in the blood
        - The log10 of the parasite density is used in the Hill function to calculate the drug effect
        - y0[dim-1] is the last element in the state vector y0 
            → it represents parasite density (per μL)
        - patient_blood_volume is in μL
        - So multiplying them gives total number of parasites
        Then takes the log10 of that 
            → common in PKPD modeling because parasitaemia can span many orders of magnitude*/

    initial_log10_totalparasitaemia = log10( y0[dim-1]*patient_blood_volume ); // Gives the log10 total parasite load in the body at t = 0
                                                                               // If intial_parasitaemia = 20000 parasiteper microliter, and patient_blood_volume = 5500000 microliters
                                                                                 // intial_total parasitaemia = 20000*5500000 = 1.1 * 10^11
                                                                                 // initial_log10_totalparasitaemia = log10(20000*5500000) = 11.041
    double THETA7_pe = 0.278;
    double THETA6_pe = -0.375;

    double TVF1 = 1.0 + THETA7_pe*(initial_log10_totalparasitaemia-3.98); // TVF1 = 1.0 + (0.278) * (11.041-3.98)
                                                                          // TVF1 = 1 + (0.278 * 7.061)
                                                                          // TVF1 = 2.96 if initial_parasitaemia = 20000 and patient_blood_volume = 5500000
    if(is_pregnant) TVF1 *= (1.0+THETA6_pe); 
    
    double F1=TVF1; // F1 = 2.96 if initial_parasitaemia = 20000 and patient_blood_volume = 5500000

    if(pkpd_artemether::stochastic) // Adjusting the between-patient variability in bioavailability (need to check if F1 is indeed bioavailability) if stochastic is true
    {
        double ETA4_rv = gsl_ran_gaussian( rng, sqrt(0.08800) );
        F1 *= ETA4_rv;
    }
 
    vprms[i_artemether_F1_indiv] = F1; 

    double TVMT_pe = 0.982;
    double MT = TVMT_pe; 
    vprms[i_artemether_KTR] = 8.0/MT; // KTR = 8/0.982 = 8.15

    double THETA1_pe = 78.0;
    double THETA2_pe = 129.0;
    double TVCL = THETA1_pe * pow( weight/mw, 0.75 );  // TVCL = 78.0 * ((weight/48.5)^0.75))
                                                       // for weight = 54, TVCL = 78.0 * ((54/48.5)^0.75) = 84.54
    double CL = TVCL;   // CL = 84.54

    double TVV2 = THETA2_pe * (weight/mw);  // TVV2 = 129.0 * (weight/48.5)
                                            // for weight = 54, TVV2 = 129.0 * (54/48.5) = 143.63
    double V2 = TVV2;   // V2 = 143.63
    if(pkpd_artemether::stochastic) 
    {
        double ETA2_rv = gsl_ran_gaussian( rng, sqrt(0.0162) ); 
        V2 *= exp(ETA2_rv); // Meaning V2 = V2 * exp(ETA2_rv)
                            // Example ETA2_rv = 0.0824, exp(0.0824) = 1.09
                            // V2 = 143.63 * 1.09 = 156.56
    }
    
    vprms[i_artemether_k20] = CL/V2; // For weight = 54. CL = 84.54, V2 = 156.56 (if stochastic)
                                    // k20 = 84.54/156.56 = 0.539

}

// Function 7: initialize_simulation

void pkpd_artemether::initialize_simulation() {
    //-- WARNING -- the age member variable must be set before you call this function -- add this check
      // Check added 04/10/25

    if (!age_initialized) {
        fprintf(stderr, "\n\tERROR: Age not initialized. Exiting program.\n");
        return;
    }

    // NOTE must call the two functions below in this order -- dosing schedule needs to be set first
    generate_recommended_dosing_schedule();
    initialize_parameters();
    
}

// Function 8: rhs_ode

//
// "rhs_ode" must be defined with these four arguments for it to work in the GSL ODE routines
// Meaning GSL expects a very specific function signature, if you don’t match it, GSL won't compile or run correctly

/* For GSL to work, your rhs_ode must:
    - Take time as a double
        → GSL gives you the current time t (in hours)

    - Take an array of state values: const double y[]
        → These are the current values of your compartments (e.g., drug amounts)

    - Take an array to write the derivatives/slopes: double f[]
        → You fill this in with how each state variable is changing at time t
        → i.e. f[i] = dy[i]/dt

    - Take a void* pointer to the object (your pkpd_artemether instance)
        → So you can access vprms, patient_blood_volume, etc.
        → We type cast this to pkpd_artemether* 
    
    - Output of this function is an int (0 for success, 1 for failure) 
        → GSL uses this to check if the function executed correctly
        → If it returns 1, GSL will stop the integration and report an error
        → If it returns 0, GSL continues with the next step */ 
        

int pkpd_artemether::rhs_ode(double t, const double y[], double f[], void *pkd_object )
{

    /* What it does:
        - pkd_object is passed in as a void*
        - GSL doesn’t know what kind of pointer it is
        - You cast it to a pointer of the correct type: pkpd_artemether* p
        - Now p can access all the class members like: p->vprms[...], p->patient_blood_volume etc.
    
    In short: This gives us a way to simulate this inside a GSL-compatible function. */

    pkpd_artemether * p = (pkpd_artemether*) pkd_object;

    // this is compartment 1, the GI tract or fixed dose compartment, i.e. the hypothetical compartment where the drug goes in first
    //
    // the value that is tracked here is the total or absolute mg amount of artemether in the patient's blood; 
    // it is NOT a drug concentration value of mg or mol artemether per unit blood volume 
    
    // f[0] tracks dy[0]/dt

    f[0] =  - p->vprms[i_artemether_KTR] * y[0]; // y[0] amount of drug in dose compartment at t = 0

     // these are the seven transit compartments
     f[1] = y[0]*p->vprms[i_artemether_KTR] - y[1]*p->vprms[i_artemether_KTR];
     f[2] = y[1]*p->vprms[i_artemether_KTR] - y[2]*p->vprms[i_artemether_KTR];
     f[3] = y[2]*p->vprms[i_artemether_KTR] - y[3]*p->vprms[i_artemether_KTR];
     f[4] = y[3]*p->vprms[i_artemether_KTR] - y[4]*p->vprms[i_artemether_KTR];
     f[5] = y[4]*p->vprms[i_artemether_KTR] - y[5]*p->vprms[i_artemether_KTR];
     f[6] = y[5]*p->vprms[i_artemether_KTR] - y[6]*p->vprms[i_artemether_KTR];
     f[7] = y[6]*p->vprms[i_artemether_KTR] - y[7]*p->vprms[i_artemether_KTR];
     
     // this is the central compartment (the blood)
     //
     // and the current units here (Aug 7 2024) are simply the total mg of artemether in the blood
     // NOTE it looks like all of our blood "concentrations" in these PK classes simply track "total mg of molecule in blood"
     f[8] = y[7]*p->vprms[i_artemether_KTR] - y[8]*p->vprms[i_artemether_k20];
     
     // this is the per/ul parasite population size
     // Testing: adjusting the concentration in the central compartment/EC50 by the PATIENT BLOOD VOLUME
     double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow((y[8]/p -> patient_blood_volume),p->pdparam_n) / (pow((y[8]/p -> patient_blood_volume),p->pdparam_n) + pow((p->pdparam_EC50/p -> patient_blood_volume),p->pdparam_n)));
 
     f[9] = -a * y[9];
     
     
     return GSL_SUCCESS;
 }


// Function 9: predict 
// 

void pkpd_artemether::predict( double t0, double t1 )
{
    gsl_odeiv_system sys = {pkpd_artemether::rhs_ode, pkpd_artemether::jac, dim, this};   // the fourth argument is a void pointer that you 
                                                                            // are supossed to use freely; you normally
                                                                            // use it to access the paramters
    double t = t0;    
    double h = 1e-6;
    
    
    while (t < t1)
    {
        // check if there are still doses to give
        if( num_doses_given < v_dosing_times.size() )
        {
            // check if time t is equal to or larger than the next scheduled dose
            if( t >= v_dosing_times[num_doses_given]  )
            {
                redraw_params_before_newdose();
                
                // add the new dose amount to the "dose compartment", i.e. the first compartment
                y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_artemether_F1_thisdose]; // Currently F1_thisdose = 1.0 in intialize_params(). Stochasticity not implemented yet.
                
                num_doses_given++;
            }
        }

        // check if time t is equal to or larger than the next scheduled hour to log
        if( t >= ((double)num_hours_logged)  )
        {
            v_concentration_in_blood.push_back( y0[8]/patient_blood_volume); // this is in mg/microliter
            v_parasitedensity_in_blood.push_back( y0[dim-1] );
            v_concentration_in_blood_hourtimes.push_back( t );
            
            num_hours_logged++;
        }

        // carry out the Runge-Kutta integration
        int status = gsl_odeiv_evolve_apply(oe, oc, os, &sys, &t, t1, &h, y0);

        if (status != GSL_SUCCESS)
            break;
    }
    
}

