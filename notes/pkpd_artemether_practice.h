
// Step 1: Include guards,  required libraries and global variables

#ifndef PKPD_ARTEMETHER
#define PKPD_ARTEMETHER

#include <vector>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h> // random number generators from Gnu Scientific Library
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>  // this includes the gsl_ran_gaussian function

/*What is an enum?
  enum (short for enumeration) is a user-defined data type that assigns names to integer values.
  It makes code more readable and avoids using "magic numbers" (hardcoded numbers with no clear meaning).*/
  
  enum parameter_index_artemether { i_artemether_KTR, i_artemether_k20, i_artemether_F1_indiv, i_artemether_F1_thisdose, artemether_num_params };

/*Description of parameter_index_artemether:
  This defines a list of named constants for indexing parameters related to artemether.
  By default, the first value (i_artemether_KTR) is assigned 0, the second (i_artemether_k20) is 1, and so on.
  The last one, artemether_num_params, will automatically be equal to the total number of parameters (in this case, 4).

 Why use enum? 
  Improves readability: Instead of writing parameters[2], you can write parameters[i_artemether_F1_indiv], which makes it clear what the value represents.
  Prevents errors: If you hardcode numbers, it's easy to make mistakes. Using named indices makes your code more robust.
  Easy to update: If more parameters are added, artemether_num_params automatically updates, so loops or arrays using it won’t break.  
  
Example:
  double parameters[artemether_num_params];  
  parameters[i_artemether_KTR] = 0.5;  // Assign value to KTR

  This makes your code much clearer than:
  parameters[0] = 0.5;  // Hard to tell what "0" represents*/


// Step 2: Define class

class pkpd_artemether {

  // Step 3: Declare public methods

public:

    // Step 3.a: Constructor
    // pkpd_artemether(); // Default

    pkpd_artemether(double age, double weight);  // constructor with inputs

    /*A constructor is a special function that runs automatically when an object of the class is created.
      It is used to initialize variables or set up the object properly.
      It has the same name as the class and no return type (not even void).*/
    
    // Step 3.b: Destructor
    ~pkpd_artemether();
    /*A destructor is a special function that runs automatically when an object is destroyed (e.g., when it goes out of scope).
      It is used to clean up resources (e.g., freeing memory if needed).
      It has the same name as the class but with a ~ (tilde) in front and no return type.*/

    // Step 4: Declare public variables
    static bool stochastic;
    /*In your case (static bool stochastic):
      Stochastic is a shared variable for all objects of pkpd_artemether.
      If you have multiple objects of the class, they will all access the same value for stochastic.
      You can modify it from anywhere, even without an object (directly using the class name).*/

    std::vector <double> vprms;

    //-----------------------------------------------------------------------------------------------------

    // Step 5: Declare methods to solve ODE

    // Step 5.a: Size or dimensionality of the ODE system
    size_t dim;
    /*Why size_t is used: 
      It is unsigned, meaning it can only hold non-negative values (positive integers or zero).
      It is typically large enough to hold any memory size, which is why it’s preferred when 
      dealing with the dimensions of large arrays, the number of elements, or the size of a data structure.*/

    gsl_odeiv_step*     os;
    gsl_odeiv_control*  oc;
    gsl_odeiv_evolve*   oe;

    // Step 5.b: Function to integrate the ODEs
    void predict(double t0, double t1);


    // Step 5.c: Add conditions for system at t = 0
    // y0 contains conditions at t = 0
    double* y0;

    // Step 5.d: Define Jacobian and right-hand side of ODEs  
    // It will compute the Jacobian matrix of your system of ODEs
    // It is typically used by numerical solvers like GSL to improve the accuracy and speed of the solution.
    static int jac(double a1, const double* a2, double* a3, double* a4, void* a5) { 
      return 0; 
    };
    /* Why we declare the jacobian in the header:
      It’s static, so it’s only used within the current translation unit.
      It may need to be inlined for performance reasons.
      It’s required to be passed as a function pointer by the ODE solver.
      It might be part of a header-only design for this part of the code.*/

    static int rhs_ode(double t, const double y[], double f[], void* params);
    /* This function computes the right-hand side of the system of ODEs at a given time t, 
       using the current state y[] and returning the derivatives in f[].*/
    
    //---------------------------------------------------------------------------------------------------------------
    
    //Step 6: Declare PK parameters via function calls

    void set_age_and_weight(double a, double w);

    void initialize();

    void initialize_parameters();

    void redraw_params_before_newdose();

    //---------------------------------------------------------------------------------------------------------------                                 

    // Step 7: PD variables and methods

    void set_parasitaemia(double parasites_per_ul);

    double pdparam_n;
    double pdparam_EC50;
    double pdparam_Pmax;

    // TODO: check this (from 2019)  Ricardo says it should be the log-natural of the per/ul parasitaemia
    //          UPDATE - 2024 - this can be deprecated
    double initial_log10_totalparasitaemia;

    //-----------------------------------------------------------------------------------------------------------------------

    // Step 8: Declare dosing schedule functions

    double num_tablets_per_dose; 
    double total_mg_dose;
    void generate_recommended_dosing_schedule();
    std::vector <double> v_dosing_times;
    std::vector <double> v_dosing_amounts;

    int num_doses_given;
    bool doses_still_remain_to_be_taken;
    bool we_are_past_a_dosing_time(double current_time);
    void give_next_dose_to_patient(double fractional_dose_taken);

    //-----------------------------------------------------------------------------------------------------------------------

    // Step 9: Set patient characteristics

    double patient_weight;          // this is the kg weight of the current patient
    double median_weight;           // this is the median weight of a patient that these estimates were calibrated for
    double weight;                  // this is the weight that is actually used in the calculations (it's one of the two above)
    double patient_age;             // parameter
    double age;                     // member/model input   
                                    // this is the age of the patient in years
    double patient_blood_volume;    // in microliters, so should be between 250,000 (infant) to 6,000,000 (large adult)
    bool weight_initialized = false;// this is true if the weight has been set
    void set_age_and_weight(double a, double w);


    //bool pregnant;                // usually means just 2nd or 3rd trimester -- TODO: have this replace the "is_pregnant" bool
    bool is_pregnant;               // TODO: deprecate
    bool is_male;




    //-----------------------------------------------------------------------------------------------------------------------

    // Step 10: Declare variables to log drug concentrations and parasite density


    // an hourly time series of drug concentrations in the blood compartment only
    std::vector <double> v_concentration_in_blood;                // an hourly time series of drug concentrations in the blood compartment only
                                                            // should be in nanograms per milliliter (ng/ml), probably, TODO: Venitha to check
                                                            // specifically for DHA is this is the case
    //vector<double> v_concentration_in_blood_metabolite;   // same as above, but this allows you to keep track of a particular metabolite concentration                                                    
                                                            // no metabolite information is used for lumefantrine
    std::vector <double> v_concentration_in_blood_hourtimes;
    std::vector <double> v_parasitedensity_in_blood;
    int num_hours_logged;
    
    gsl_rng *rng;	
    
  };

  #endif // PKPD_ARTEMETHER