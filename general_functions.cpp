#include "pkpd_ppq.h" 
#include "pkpd_dha.h"
#include "pkpd_lum.h"
#include "pkpd_adq.h"
#include "pkpd_artemether.h"
#include "general_functions.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_rng.h> // random number generators from Gnu Scientific Library
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <random>
#include <string>
#include <string>

double G_CLO_PMF = 8.0;
double G_DENSITY_50 = 50011.087;

int G_CLO_N = 1;
double G_CLO_AGE = 25.0;
double G_CLO_WEIGHT = 54.0;

double G_CLO_HILL_COEFF_DHA = 20.0;
double G_CLO_HILL_COEFF_ARTEMETHER = 20.0;
double G_CLO_HILL_COEFF_LUM = 15.0;

double G_CLO_EC50_DHA = 0.1;
double G_CLO_EC50_ARTEMETHER = 0.1;                                               // Assuming this to be ng/ml
double G_CLO_EC50_LUM = exp( 0.525 * log (2700)); // use natural log, 63.30907617 // Assuming this to be ng/ml

double G_CLO_PMAX_DHA = 0.983; //pmax_dha = 0.983 gives ~68.9% efficacy for DHA monotherapy, calibrated by Venitha in Dec 2024 
                               //Original value = 0.99997

double G_CLO_PMAX_ARTEMETHER = 0.983; 
double G_CLO_PMAX_LUM = 0.9995;

int G_CLO_OUTPUT_TYPE = 0;

// Trying to study why AL efficacy changes with weight
double G_CLO_BLOOD_VOLUME_EXPONENT = 1;
double G_CLO_CENTRAL_VOLUME_EXPONENT = 1;

enum therapy_type { therapy_none , therapy_amodiaquine, therapy_artemether, therapy_artesunate , therapy_dha, therapy_lumefantrine , therapy_AL, therapy_DHA_PPQ, therapy_AS_AQ}; 
enum therapy_type G_CLO_THERAPY = therapy_none;

void output_results_monotherapy_dha(int pi, pkpd_dha *dyn)
{
    if (G_CLO_OUTPUT_TYPE == 1) {
        int j = dyn->v_concentration_in_blood.size()-1;
        fprintf(stdout, "%d , %10.3f , %10.3f , %10.3f \n", pi, dyn->v_concentration_in_blood_hourtimes[j], dyn->v_concentration_in_blood[j], dyn->v_parasitedensity_in_blood[j] );
    }
    else {
        for(int j=0; j<dyn->v_concentration_in_blood.size(); j++ )
        {
            fprintf(stdout, "%d , %10.3f , %10.3f , %10.3f \n", pi, dyn->v_concentration_in_blood_hourtimes[j], dyn->v_concentration_in_blood[j], dyn->v_parasitedensity_in_blood[j] );
        }
    }
}

void output_results_combination_AL(int pi, pkpd_artemether *dyn1, pkpd_lum *dyn2)
{
    if (G_CLO_OUTPUT_TYPE == 1) {
        int j = dyn1->v_concentration_in_blood.size()-1;
        fprintf(stdout, "%d %10.3f %10.3f %10.3f %10.3f \n", pi, dyn1->v_concentration_in_blood_hourtimes[j], dyn1->v_concentration_in_blood[j], dyn2->v_concentration_in_blood[j], dyn1->v_parasitedensity_in_blood[j] );        
    }
    else {
        for(int j=0; j<dyn1->v_concentration_in_blood.size(); j++ )
        {
            fprintf(stdout, "%d %10.3f %10.3f %10.3f %10.3f \n", pi, dyn1->v_concentration_in_blood_hourtimes[j], dyn1->v_concentration_in_blood[j], dyn2->v_concentration_in_blood[j], dyn1->v_parasitedensity_in_blood[j] );
        }
    }
}

void output_results_combination_DHA_PPQ(int pi, pkpd_dha *dyn1, pkpd_ppq *dyn2)
{
    if (G_CLO_OUTPUT_TYPE == 1) {
        int j = dyn1->v_concentration_in_blood.size()-1;
        fprintf(stdout, "%d %10.3f %10.3f %10.3f %10.3f \n", pi, dyn1->v_concentration_in_blood_hourtimes[j], dyn1->v_concentration_in_blood[j], dyn2->v_concentration_in_blood[j], dyn1->v_parasitedensity_in_blood[j] );        
    }
    else {
        for(int j=0; j<dyn1->v_concentration_in_blood.size(); j++ )
        {
            fprintf(stdout, "%d %10.3f %10.3f %10.3f %10.3f \n", pi, dyn1->v_concentration_in_blood_hourtimes[j], dyn1->v_concentration_in_blood[j], dyn2->v_concentration_in_blood[j], dyn1->v_parasitedensity_in_blood[j] );
        }
    }
}

void ParseArgs(int argc, char **argv)
{
    int i, start;
    start=1;
    
    string str;
    for(i=start; i<argc; i++)
    {
        str = argv[i];
        

        if( str == "--adq" ) 		                    G_CLO_THERAPY  		                = therapy_amodiaquine;
        else if( str == "--artemether" )		        G_CLO_THERAPY  		                = therapy_artemether;
        else if( str == "--artesunate" )		        G_CLO_THERAPY  		                = therapy_artesunate;
        else if( str == "--dha" ) 		                G_CLO_THERAPY  		                = therapy_dha;
        else if( str == "--lum" )		                G_CLO_THERAPY  		                = therapy_lumefantrine;

        else if( str == "--AL" ) 		                G_CLO_THERAPY  		                = therapy_AL;
        else if( str == "--DHA-PPQ" ) 		            G_CLO_THERAPY  		                = therapy_DHA_PPQ;
        else if( str == "--AS-AQ" ) 		            G_CLO_THERAPY  		                = therapy_AS_AQ;

        else if( str == "--age" ) 		                G_CLO_AGE  		                    = atof( argv[++i] );
        else if( str == "--weight" ) 	                G_CLO_WEIGHT	                    = atof( argv[++i] );
        else if( str == "-n" ) 	                        G_CLO_N	                            = atoi( argv[++i] );
        else if( str == "--pmf" ) 	                    G_CLO_PMF	                        = atof( argv[++i] );
        
        else if( str == "--hill_dha" ) 		            G_CLO_HILL_COEFF_DHA	            = atof( argv[++i] );
        else if( str == "--hill_artemether" ) 		    G_CLO_HILL_COEFF_ARTEMETHER	        = atof( argv[++i] );
        else if( str == "--hill_lum" ) 		            G_CLO_HILL_COEFF_LUM	            = atof( argv[++i] );
        
        else if( str == "--ec50_dha" ) 		            G_CLO_EC50_DHA	                    = atof( argv[++i] );
        else if( str == "--ec50_artemether" ) 		    G_CLO_EC50_ARTEMETHER	            = atof( argv[++i] );
        else if( str == "--ec50_lum" ) 		            G_CLO_EC50_LUM	                    = atof( argv[++i] );

        else if( str == "--pmax_dha" ) 		            G_CLO_PMAX_DHA	                    = atof( argv[++i] );
        else if (str == "--pmax_artemether" ) 		    G_CLO_PMAX_ARTEMETHER               = atof( argv[++i] );
        else if (str == "--pmax_lum" ) 		            G_CLO_PMAX_LUM	                    = atof( argv[++i] );

        // BV exponent can be removed
        else if( str == "--blood_volume_exponent" )     G_CLO_BLOOD_VOLUME_EXPONENT	        = atof( argv[++i] );
        else if (str == "--central_volume_exponent" )   G_CLO_CENTRAL_VOLUME_EXPONENT       = atof( argv[++i] );

        else if (str == "-o" ) 		                    G_CLO_OUTPUT_TYPE	                = atoi( argv[++i] );
        
        /*else if( str == "--endttr" ) 		            G_CLO_END_TITER	                    = atof( argv[++i] );
        else if( str == "--chainlength" )               G_CLO_CHAIN_LENGTH	                = atoi( argv[++i] );
        else if( str == "--showevery" ) 		        G_CLO_SHOW_EVERY	                = atoi( argv[++i] );
        else if( str == "--burnin" ) 		            G_CLO_BURNIN	                    = atoi( argv[++i] );
        else if( str == "--iltau" ) 		            G_CLO_INTEGRATION_LIMIT_FOR_TAU	    = atof( argv[++i] );
        else if( str == "--profile" ) 		            G_CLO_PROFILE  		                = true;
        else if( str == "--censored" ) 		            G_CLO_CENSORED 		                = true;*/
	else
        {
            fprintf(stderr, "\n\tUnknown option [%s] on command line.\n\n", argv[i]);
            exit(-1);
        }
    }

    return;
}