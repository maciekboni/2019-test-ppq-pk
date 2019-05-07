#ifndef PKPD_DHA
#define PKPD_DHA

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
enum parameter_index_dha { i_dha_KTR, i_dha_k20, i_dha_F1_indiv, i_dha_F1_thisdose, dha_num_params }; 

//typedef enum parameter_index_dha i_dha;

//extern gsl_rng *G_RNG;	

class pkpd_dha
{   
public:    
    explicit pkpd_dha();    // constructor
    ~pkpd_dha();         	  // destructor

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

    //void initialize();
    //void initialize_params_w_population_means();
    void initialize_params();

    // TODO implement this function; we should remove the dosing from the predict function
    void give_patient_new_dose( double fraction_of_dose_taken );
    void redraw_params_before_newdose();
    
    double pdparam_n;
    double pdparam_EC50;
    double pdparam_Pmax;
    
    //TODO check this.  Ricardo says it should be the log-natural of the per/ul parasitaemia
    double initial_log10_totalparasitaemia;
    
    
    // the members below are used to create and manage the dose schedule; i.e. the course of treatment
    void generate_recommended_dosing_schedule();
    vector<double> v_dosing_times;
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
    bool is_male;
    bool is_pregnant;

    double patient_blood_volume; // in microliters
    
    gsl_rng *rng;		
    
};

#endif // PKPD_DHA
