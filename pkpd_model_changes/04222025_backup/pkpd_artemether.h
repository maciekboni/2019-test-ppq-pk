#ifndef PKPD_artemether
#define PKPD_artemether

#include <vector>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h> // random number generators from Gnu Scientific Library
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>  // this includes the gsl_ran_gaussian function

using namespace std;


// these are the parameters
enum parameter_index_artemether { i_artemether_KTR, i_artemether_k20, i_artemether_F1_indiv, i_artemether_F1_thisdose, i_artemether_typical_CL, i_artemether_CL_indiv, i_artemether_typical_V2, i_artemether_V2_indiv, artemether_num_params }; 

//typedef enum parameter_index_artemether i_artemether;

//extern gsl_rng *G_RNG;	

class pkpd_artemether
{   
public:    
    explicit pkpd_artemether();    // constructor
    ~pkpd_artemether();         	  // destructor

                                  
    static bool stochastic;       // just set this to false and the class will run a deterministic model with
                                  // all the population means; age and weight dependence will still be in the model

    vector<double> vprms;         // this holds all the parameters; they are indexed by the enums above
    
    //
    // ----  1  ----  ODE STRUCTURE
    //

    size_t dim;		          // the dimensionality of the ODE system; this has to be of type size_t
                              // if you want to avoid warnings; it's just a "long unsigned int"

    
    // these are pointers to GSL ODE structures that integrate the ODE system
    gsl_odeiv_step* 	os;
    gsl_odeiv_control* 	oc;
    gsl_odeiv_evolve*	oe;


    // integrate the ODEs from t0 to t1; the conditions at t0 are stored in the class as member y0;
    // after predict is done, it updates y0 automatically, so you can continue integrating
    void predict( double t0, double t1 );
    double* y0;
    
    static int jac(double a1, const double* a2, double* a3, double* a4, void* a5) { return 0; };
    static int rhs_ode(double t, const double y[], double f[], void *params);



    //
    // ----  2  ----  PK PARAMETERS
    //

    void initialize();  
    //void initialize_params_w_population_means();
    void initialize_params();

    // TODO: implement this function; we should remove the dosing from the predict function
    void redraw_params_before_newdose();

    

    //
    // ----  3  ----  PD PARAMETERS that appear in the Hill function
    //

    void set_parasitaemia( double parasites_per_ul ); 

    double pdparam_n;
    double pdparam_EC50;
    double pdparam_Pmax;
    
    // TODO: check this (from 2019)  Ricardo says it should be the log-natural of the per/ul parasitaemia
    //          UPDATE - 2024 - this can be deprecated
    double initial_log10_totalparasitaemia;


    //
    // ----  4  ----  DOSING SCHEDULE
    //

    // the members below are used to create and manage the dose schedule; i.e. the course of treatment
    void generate_recommended_dosing_schedule();
    vector<double> v_dosing_times;
    vector<double> v_dosing_amounts;
    int num_doses_given;
    bool doses_still_remain_to_be_taken;
    bool we_are_past_a_dosing_time( double current_time );    
    void give_next_dose_to_patient( double fraction_of_dose_taken ); // TODO: this needs to be implemented



    //
    // ----  5  ----  PATIENT CHARACTERISTICS
    //

    double patient_weight;  // this is the kg weight of the current patient
    double median_weight;   // this is the median weight of a patient that these estimates were calibrated for
    double weight;         // this is the weight that is actually used in the calculations (it's one of the two above)
    double age;
    double patient_blood_volume;    // in microliters, so should be between 250,000 (infant) to 6,000,000 (large adult)
    //bool pregnant;                  // usually means just 2nd or 3rd trimester -- TODO: have this replace the "is_pregnant" bool
    bool is_pregnant; // TODO: deprecate
    bool is_male;



    //
    // ----  6  ----  STORAGE VARIABLES FOR DYNAMICS OF PK AND PD CURVES
    //


    // vector<double> v_dosing_compartment;                // an hourly time series of drug concentrations in the dosing compartment
    // vector<double> v_transit_compartment1;              // an hourly time series of drug concentrations in the transit compartment 1
    // vector<double> v_transit_compartment2;              // an hourly time series of drug concentrations in the transit compartment 2
    // vector<double> v_transit_compartment3;              // an hourly time series of drug concentrations in the transit compartment 3
    // vector<double> v_transit_compartment4;              // an hourly time series of drug concentrations in the transit compartment 4
    // vector<double> v_transit_compartment5;              // an hourly time series of drug concentrations in the transit compartment 5
    // vector<double> v_transit_compartment6;              // an hourly time series of drug concentrations in the transit compartment 6
    // vector<double> v_transit_compartment7;              // an hourly time series of drug concentrations in the transit compartment 7

    // vector<double> v_killing_rate;                            // an hourly time series of the killing rate
    
    // an hourly time series of drug concentrations in the blood compartment only
    vector<double> v_concentration_in_blood;                // an hourly time series of drug concentrations in the blood compartment only
                                                            // should be in nanograms per milliliter (ng/ml), probably, TODO: Venitha to check
                                                            // specifically for artemether is this is the case
    
                                                            //vector<double> v_concentration_in_blood_metabolite;   // same as above, but this allows you to keep track of a particular metabolite concentration                                                    
                                                            // no metabolite information is used for lumefantrine
    vector<double> v_concentration_in_blood_hourtimes;
    vector<double> v_parasitedensity_in_blood;
    int num_hours_logged;
    double total_mg_dose_per_occassion;
    
    gsl_rng *rng;		
    
};

#endif // PKPD_artemether
