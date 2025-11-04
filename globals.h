#ifndef GLOBALS_H
#define GLOBALS_H

#include <random>
#include <gsl/gsl_rng.h>

enum therapy_type { 
    therapy_none, 
    therapy_amodiaquine, 
    therapy_artemether, 
    therapy_artesunate, 
    therapy_dha, 
    therapy_lumefantrine, 
    therapy_AL, 
    therapy_DHA_PPQ, 
    therapy_AS_AQ}; 

// Global random number generators
extern std::mt19937 G_RNG_CPP;      // C++ standard library RNG
extern gsl_rng* G_RNG;              // GSL RNG

// Global simulation control variables
extern double maximum_enforced_stepsize;
extern double t0, t1;

// Member variables corresponding to command-line options
extern double G_CLO_PMF;
extern double G_DENSITY_50;
extern int G_CLO_N;
extern double G_CLO_AGE;
extern double G_CLO_WEIGHT;

// Pmax values
extern double G_CLO_PMAX_DHA;
extern double G_CLO_PMAX_ARTEMETHER;
extern double G_CLO_PMAX_ARTESUNATE;
/*--------------------------------*/
extern double G_CLO_PMAX_AQ;
extern double G_CLO_PMAX_LUM;
extern double G_CLO_PMAX_PPQ;

// EC50 values
// Assuming this to be ng/ml
extern double G_CLO_EC50_DHA;
extern double G_CLO_EC50_ARTEMETHER;                                               
extern double G_CLO_EC50_ARTESUNATE; 
/*--------------------------------*/
extern double G_CLO_EC50_AQ;
extern double G_CLO_EC50_LUM;
extern double G_CLO_EC50_PPQ;

// Hill Co-efficient
extern double G_CLO_HILL_COEFF_DHA;
extern double G_CLO_HILL_COEFF_ARTEMETHER;
extern double G_CLO_HILL_COEFF_ARTESUNATE;
/*--------------------------------*/
extern double G_CLO_HILL_COEFF_AQ;
extern double G_CLO_HILL_COEFF_LUM;
extern double G_CLO_HILL_COEFF_PPQ;

extern int G_CLO_OUTPUT_TYPE;
extern therapy_type G_CLO_THERAPY;

#endif // GLOBALS_H