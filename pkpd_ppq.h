#ifndef PKPD_PPQ
#define PKPD_PPQ

// #include "all_gp_data.h"
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


// these are the parameters -- the n intercept parameters (for the n slopes for the n patients) are not included in this list
enum parameter_index_ppq { i_ppq_k15, i_ppq_k62, i_ppq_k20, i_ppq_k23, i_ppq_k32, i_ppq_k24, i_ppq_k42, i_ppq_k56, i_ppq_F1_indiv, i_ppq_F1_thisdose, ppq_num_params }; 

// NOTE -- there is no need to declare the typedef below; in practice, you always use integers, and you never cycle through the parameters
// typedef enum parameter_index_ppq i_ppq;

//extern gsl_rng *G_RNG;	

class pkpd_ppq
{   
public:    
    explicit pkpd_ppq();    // constructor
    ~pkpd_ppq();         	  // destructor

    size_t dim;		          // the dimensionality of the ODE system; this has to be of type size_t
                              // if you want to avoid warnings; it's just a "long unsigned int"
                                  
    static bool stochastic;       // just set this to false and the class will run a deterministic model with
                                  // all the population means; age and weight dependence will still be in the model

    vector<double> vprms;         // this holds all the parameters; they are indexed by the enums above
    
    
    // integrate the ODEs from t0 to t1; the conditions at t0 are stored in the class as member y0;
    // after predict is done, it updates y0 automatically, so you can continue integrating
    void predict( double t0, double t1 );
    double* y0;
    
    static int jac(double a1, const double* a2, double* a3, double* a4, void* a5) { return 0; };
    static int rhs_ode(double t, const double y[], double f[], void *params);

    void initialize();
    void initialize_params_w_population_means();
    void initialize_params();

    void set_parasitaemia( double parasites_per_ul );
    
    void redraw_params_before_newdose();
    bool we_are_past_a_dosing_time( double current_time );    
    void give_next_dose_to_patient( double fractional_dose_taken );

    // these are the pharmacodynamic parameters that appear in this hill function
    double pdparam_n;
    double pdparam_EC50;
    double pdparam_Pmax;
    
    
    // the members below are used to create and manage the dose schedule; i.e. the course of treatment
    void generate_recommended_dosing_schedule();
    vector<double> v_dosing_times; // in hours
    vector<double> v_dosing_amounts;
    int num_doses_given;

    // an hourly time series of drug concentrations in the blood compartment only
    vector<double> v_concentration_in_blood;
    vector<double> v_concentration_in_blood_hourtimes;
    vector<double> v_parasitedensity_in_blood;
    int num_hours_logged;
        
    // these are pointers to GSL ODE structures that integrate the ODE system
    gsl_odeiv_step* 	os;
    gsl_odeiv_control* 	oc;
    gsl_odeiv_evolve*	oe;
    
    double patient_weight;  // this is the kg weight of the current patient
    double median_weight;   // this is the median weight of a patient that these estimates were calibrated for
    double weight;         // this is the weight that is actually used in the calculations (it's one of the two above)
    double age;

    double patient_blood_volume; // in microliters
    
    gsl_rng *rng;		
    
};

#endif // PKPD_PPQ
