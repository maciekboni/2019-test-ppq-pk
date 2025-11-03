#include "globals.h"
#include "general_functions.h"

// GLOBAL RANDOM NUMBER GENERATOR
gsl_rng *G_RNG;	
// gsl_rng* G_RNG = nullptr;
std::mt19937 G_RNG_CPP;	// Change to pointer if required

// Parasite dynamics
double G_CLO_PMF = 8.0;
double G_DENSITY_50 = 50011.087;   // from correspondence with Aubrey Cunnington, the parasite density level at which growth is inhibited 
                            // to 50% of its max value occurs at ln(10.82) (10.49), 11.54) parasites per microliter for 
                            // estimated from n=64 Gambian children with uncomplicated malaria (Giorgiadou et al, Nat Microbiol 2019)

/*************************************************************************************************************************************/
// Patient characteristics
double G_CLO_N = 1;
double G_CLO_AGE = 25.0;
double G_CLO_WEIGHT = 54.0;

/*************************************************************************************************************************************/

// PK/PD parameters
    // Pmax values
double G_CLO_PMAX_DHA = 0.983; //pmax_dha = 0.983 gives ~68.9% efficacy for DHA monotherapy, calibrated by Venitha in Dec 2024
                        //Original value = 0.99997
double G_CLO_PMAX_ARTEMETHER = 0.983;
double G_CLO_PMAX_ARTESUNATE = 0.983;
/*--------------------------------*/
double G_CLO_PMAX_AQ = 0.95;
double G_CLO_PMAX_LUM = 0.97;
double G_CLO_PMAX_PPQ = 0.95;

// EC50 values
// Assuming this to be ng/ml
double G_CLO_EC50_DHA = 0.1;
double G_CLO_EC50_ARTEMETHER = 0.1;                                               
double G_CLO_EC50_ARTESUNATE = 0.1;
/*--------------------------------*/
double G_CLO_EC50_AQ = 5.0;
double G_CLO_EC50_LUM = 175.0;
double G_CLO_EC50_PPQ = 5.0;

// Hill Co-efficient
double G_CLO_HILL_COEFF_DHA = 20.0;
double G_CLO_HILL_COEFF_ARTEMETHER = 20.0;
double G_CLO_HILL_COEFF_ARTESUNATE = 20.0;
/*--------------------------------*/
double G_CLO_HILL_COEFF_AQ = 5.0;
double G_CLO_HILL_COEFF_LUM = 6.0;
double G_CLO_HILL_COEFF_PPQ = 5.0;

/*************************************************************************************************************************************/

// Input/Output parameters
double G_CLO_OUTPUT_TYPE = 0;
therapy_type G_CLO_THERAPY = therapy_none;

// Define simulation control variables
double maximum_enforced_stepsize = 0.5;
double t0 = 0.0;
double t1 = maximum_enforced_stepsize;