#include "general_functions.h"
#include "globals.h"

#include "pkpd_artemether.h"
#include "pkpd_artesunate.h"
#include "pkpd_dha.h"

#include "pkpd_adq.h"
#include "pkpd_lum.h"
#include "pkpd_ppq.h" 

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


general_functions::general_functions()
{
      // ✅ Now this will work
}

general_functions::~general_functions()
{
      // ✅ Now this will work
}

void general_functions::output_results_monotherapy_dha(int pi, pkpd_dha *dyn)
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

void general_functions::output_results_monotherapy_lum(int pi, pkpd_lum *dyn)
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

void general_functions::output_results_combination_AL(int pi, pkpd_artemether *dyn1, pkpd_lum *dyn2)
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

void general_functions::output_results_combination_DHA_PPQ(int pi, pkpd_dha *dyn1, pkpd_ppq *dyn2)
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

void general_functions::output_results_combination_AS_AQ(int pi, pkpd_artesunate *dyn1, pkpd_adq*dyn2)
{
    if (G_CLO_OUTPUT_TYPE == 1) {
        int j = dyn1->v_concentration_in_blood.size()-1;
        fprintf(stdout, "%d %10.3f %10.3f %10.3f %10.3f %10.3f \n", pi, dyn1->v_concentration_in_blood_hourtimes[j], dyn1->v_concentration_in_blood[j], dyn2->v_concentration_in_blood[j], dyn2->v_concentration_in_blood_metabolite[j], dyn1->v_parasitedensity_in_blood[j] );        
    }
    else {
        for(int j=0; j<dyn1->v_concentration_in_blood.size(); j++ )
        {
            fprintf(stdout, "%d %10.3f %10.3f %10.3f %10.3f %10.3f \n", pi, dyn1->v_concentration_in_blood_hourtimes[j], dyn1->v_concentration_in_blood[j], dyn2->v_concentration_in_blood[j], dyn2->v_concentration_in_blood_metabolite[j], dyn1->v_parasitedensity_in_blood[j] );
        }
    }
}

void general_functions::ParseArgs(int argc, char **argv)
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