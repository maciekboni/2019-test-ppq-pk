#ifndef PKPD_LUM
#define PKPD_LUM

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


// these are the flow or rate parameters in the ODEs, and some additional dose and scaling parameters
enum parameter_index_lum { i_lum_k12, i_lum_k23, i_lum_k32,  i_lum_k20, i_lum_F1_indiv, i_lum_F1_thisdose, i_lum_pmf, lum_num_params }; 


class pkpd_lum
{   
public:    
    explicit pkpd_lum();    // constructor
    ~pkpd_lum();            // destructor

    gsl_rng *rng;		

    static bool stochastic;       // just set this to false and the class will run a deterministic model with
                                  // all the population means; age and weight dependence will still be in the model



    //
    // ----  1  ----  ODE STRUCTURE
    //

    size_t dim;		          // the dimensionality of the ODE system; this has to be of type size_t
                              // if you want to avoid warnings; it's just a "long unsigned int"
                                  
    // these are pointers to GSL ODE structures that integrate the ODE system
    gsl_odeiv_step* os;
    gsl_odeiv_control* oc;
    gsl_odeiv_evolve* oe;

    // integrate the ODEs from t0 to t1; the conditions at t0 are stored in the class as member y0;
    // after predict is done, it updates y0 automatically, so you can continue integrating
    void predict( double t0, double t1 );
    double* y0;
    
    static int jac(double a1, const double* a2, double* a3, double* a4, void* a5) { return 0; };
    static int rhs_ode(double t, const double y[], double f[], void *params);



    //
    // ----  2  ----  PK PARAMETERS
    //

    vector<double> vprms;         // this holds all the parameters; they are indexed by the enums above
    
    void initialize();
    void initialize_params();
    void redraw_params_before_newdose();    // you may need this if there is inter-occassion variability
                                            // in the PK parameters 

    double central_volume_of_distribution;  // for a particular patient of a particular weight, this is the volume of the blood 
                                            // plus the volume of everything else that is in instantaneous equilibrium with the blood
    


    //
    // ----  3  ----  PD PARAMETERS that appear in the Hill function
    //
    void set_parasitaemia( double parasites_per_ul );
    double parasites_per_ul_at_first_lum_dose;  // this is needed explicitly for lumefantrine because it is a covariate that affects absorption

    double pdparam_n;
    double pdparam_EC50;
    double pdparam_Pmax;



    //
    // ----  4  ----  DOSING SCHEDULE
    //

    void generate_recommended_dosing_schedule();

    vector<double> v_dosing_times;      // in hours
    vector<double> v_dosing_amounts;    // matched to the dosing times above

    int num_doses_given;
    double total_mg_dose_per_occassion; // this is the same as daily dose for most therapies
                                        // but NOT for artemether-lumefantrine or ALAQ
    bool doses_still_remain_to_be_taken;
    bool we_are_past_a_dosing_time( double current_time );    
    void give_next_dose_to_patient( double fractional_dose_taken );



    //
    // ----  5  ----  PATIENT CHARACTERISTICS
    //

    double patient_weight;          // this is the kg weight of the current patient
    double median_weight;           // this is the median weight of a patient that these estimates were calibrated for
    double weight;                  // this is the weight that is actually used in the calculations (it's one of the two above)
    double age;                     // in years
    double patient_blood_volume;    // in microliters, so should be between 250,000 (infant) to 6,000,000 (large adult)
    bool pregnant;                  // usually means just 2nd or 3rd trimester

    void set_age_and_weight(double a, double w);

    //
    // ----  6  ----  STORAGE VARIABLES FOR DYNAMICS OF PK AND PD CURVES
    //

    
    vector<double> v_dosing_compartment;
    vector<double> v_peripheral_concentration;
    vector<double> v_killing_rate;

    vector<double> v_concentration_in_blood;                // an hourly time series of drug concentrations in the blood compartment only
                                                            // should be in nanograms per milliliter (ng/ml)
    //vector<double> v_concentration_in_blood_metabolite;   // same as above, but this allows you to keep track of a particular metabolite concentration                                                    
                                                            // no metabolite information is used for lumefantrine
    vector<double> v_concentration_in_blood_hourtimes;      // this just gives you an x-axis for plotting
    vector<double> v_parasitedensity_in_blood;
    int num_hours_logged;
    
};

#endif // PKPD_LUM
