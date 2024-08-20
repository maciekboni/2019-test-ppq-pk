#include "pkpd_ppq.h" 
#include "pkpd_dha.h"
#include "pkpd_lum.h"
#include "pkpd_adq.h"


#include <gsl/gsl_rng.h> // random number generators from Gnu Scientific Library
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <chrono>
#include <string>

// GLOBAL RANDOM NUMBER GENERATOR
gsl_rng *G_RNG;		

// GLOBAL VARIABLES
bool G_CLO_LUM = false;
bool G_CLO_ADQ = false;

enum therapy_type { none , therapy_artemisinin , therapy_lumefantrine , therapy_amodiaquine , therapy_AL }; 
enum therapy_type G_CLO_THERAPY = none;

double G_CLO_AGE = 25.0;
double G_CLO_WEIGHT = 54.0;

// PARASITE MULTIPLICATION FACTOR (PMF) - this is the number of parasites that are produced by each parasite in the blood per 48h life cycle
// The default value is 10.0 for a 48h cycle. Change this in the differential equations to make sure it is scaled correctly
double G_CLO_PMF = 1.0;

int G_CLO_N = 1; // this is the number of patients

// from correspondence with Aubrey Cunnington, the parasite density level at which growth is inhibited 
// to 50% of its max value occurs at ln(10.82) (10.49), 11.54) parasites per microliter 
// estimated from n=64 Gambian children with uncomplicated malaria (Giorgiadou et al, Nat Microbiol 2019)
double G_DENSITY_50 = 50011.087; // calculated as (e^10.82)

// Adding the following parameters for customizing the hill coefficient, EC50 and Pmax for artemisinin and lumefantrine

// enum hill_coefficient { none , hill_art, hill_lum };
// enum hill_coefficient G_CLO_HILL_COEFF = none;

// enum ec50 { none , ec50_art, ec50_lum };
// enum ec50 G_CLO_EC50 = none;

// enum pmax { none, pmax_art, pmax_lum };
// enum pmax G_CLO_PMAX = none;

// FUNCTION DECLARATIONS
void ParseArgs(int argc, char **argv);

int main(int argc, char* argv[])
{
    ParseArgs(argc,argv);

    // get random number seed from current time
    auto now = std::chrono::high_resolution_clock::now();
    auto milliseconds =  std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch());
    auto seed = static_cast<unsigned long>(milliseconds.count());
    //printf("\n seed : %lu\n", seed);
    
    // make random number generator (RNG) the Mersenne Twister which has period 2^19937 - 1
    const gsl_rng_type *TT_RAND = gsl_rng_mt19937;
    G_RNG = gsl_rng_alloc(TT_RAND);
    gsl_rng_set( G_RNG, seed ); // seed the RNG    

    // these will always be stochastic unless we are debugging something    
    pkpd_dha::stochastic = true;
    pkpd_ppq::stochastic = true;
    pkpd_lum::stochastic = true;
    pkpd_adq::stochastic = true;
    
    // declaring and initializing the standard time-keeping variables for the main ODE loop
    double maximum_enforced_stepsize = 0.5; // in hours
    double t0=0.0;
    double t1=maximum_enforced_stepsize;
    
    /*while( t0 < 52.0 )
    {
        dyn_dha->predict(t0, t1);
	t0+=maximum_enforced_stepsize; t1+=maximum_enforced_stepsize;
    }
    
        
    for(int j=0; j<dyn_dha->v_concentration_in_blood.size(); j++ )
    {
        if(j%1==0) printf("%10.3f \t %10.3f \n", dyn_dha->v_concentration_in_blood_hourtimes[j], dyn_dha->v_concentration_in_blood[j] );
    }*/

    
    
    
    //FILE *fp, *fp2, *fp3, *fp4;


    //fp = fopen("out.ppq.allpatients.v20240318.csv","w");
    //fp2 = fopen("out.dha.allpatients.v20240318.csv","w");
    
    //fp4 = fopen("out.adq.allpatients.v20240328.csv","w");
    //fprintf(fp, "PID,HOUR,COMP2CONC,PARASITEDENSITY\n" );
    //fprintf(fp2, "PID,HOUR,COMP2CONC,PARASITEDENSITY\n" );
    
    //fprintf(fp4, "PID,HOUR,COMP2CONC,PARASITEDENSITY\n" );
    
    // pi is patient index
    /*for(int pi=0; pi < 100; pi++)
    {
        auto dyn = new pkpd_ppq();
        dyn->rng = G_RNG;    
        dyn->age = 25.0;
        dyn->initialize();
        dyn->set_parasitaemia( 1000.0 );
        t0=0.0;
        t1=maximum_enforced_stepsize;
        

        //BEGIN - INTEGRATION
        while( t0 < 168.0*4.0 )
        {
            if( dyn->we_are_past_a_dosing_time(t0) )
            {
                dyn->give_next_dose_to_patient(1.0);    // 1.0 means the full dose is given
                                                        // if no dose remains to be given, function does nothing
            }

            dyn->predict(t0, t1);
            t0+=maximum_enforced_stepsize; t1+=maximum_enforced_stepsize;
        }
        //END - INTEGRATION 

        for(int j=0; j<dyn->v_concentration_in_blood.size(); j++ )
        {
            fprintf(fp, "%d , %10.3f , %10.3f , %10.3f \n", pi, dyn->v_concentration_in_blood_hourtimes[j], dyn->v_concentration_in_blood[j], dyn->v_parasitedensity_in_blood[j] );
        }
        
        delete dyn;
    }*/ 


    if( G_CLO_THERAPY == therapy_lumefantrine )
    {
        //fp3 = fopen("out.lum.allpatients.20240401.csv","w");
        fprintf(stdout, "PID,HOUR,COMP2CONC,PARASITEDENSITY\n" );

        fprintf(stderr, "\n");
        // pi is patient index
        for(int pi=0; pi < G_CLO_N; pi++)
        {
            auto dyn = new pkpd_lum();
            //fprintf(stderr, "\tlum object created pi = %d \r", pi); fflush(stderr);
            dyn->set_parasitaemia( 20000.0 );    
            dyn->parasites_per_ul_at_first_lum_dose = 1000.0;   // YOU MUST DO THIS SEPARATELY because the parasitaemia level "at first
                                                            // lum dose" is a special quantity that affects the lum absorption

            dyn->rng = G_RNG;    
            dyn->age = G_CLO_AGE;
            dyn->weight = G_CLO_WEIGHT;
            dyn->initialize();                              // NB: parasitaemia must be set before initializing parameters
            //dyn->vprms[i_lum_pmf] = G_CLO_PMF; // TODO: deprecate this pmf enum
            //fprintf(stderr, "\n\tlum object initialized pi=%d", pi); fflush(stderr);
        
            t0=0.0;
            t1=maximum_enforced_stepsize;
            double stepsize_PMF = pow( G_CLO_PMF, 1.0 / (48.0/maximum_enforced_stepsize) );

            //BEGIN - INTEGRATION
            while( t0 < 168.0*4.0 )
            {
                if( dyn->doses_still_remain_to_be_taken )
                {
                    if( dyn->we_are_past_a_dosing_time(t0) )   
                    {
                        dyn->give_next_dose_to_patient(1.0);    // 1.0 means the full dose is given
                                                                // if no dose remains to be given, function does nothing
                    }
                }
                dyn->predict(t0, t1);

                // after integrating the differential equations in the predict function above,
                // we need to ---- GROW THE PARASITES ---- for half-an-hour (i.e the maximum_enforced_stepsize)   
                dyn->y0[ dyn->dim - 1 ] *= stepsize_PMF; 

                t0 += maximum_enforced_stepsize; t1 += maximum_enforced_stepsize;
            }
            //END - INTEGRATION 

            //for(int j=0; j<dyn->v_concentration_in_blood.size(); j++ )
            for(int j=0; j<96; j++ )
            {
                fprintf(stdout, "%d , %10.3f , %10.3f , %10.3f \n", pi, dyn->v_concentration_in_blood_hourtimes[j], dyn->v_concentration_in_blood[j], dyn->v_parasitedensity_in_blood[j] );
            }
        
            delete dyn;
        }

        fprintf(stderr, "\n");
        //fclose(fp3);
    }


    if( G_CLO_THERAPY == therapy_artemisinin )
    {
        //fp3 = fopen("out.lum.allpatients.20240401.csv","w");
        fprintf(stdout, "PID,HOUR,COMP2CONC,PARASITEDENSITY\n" );

        fprintf(stderr, "\n");
        // pi is patient index
        for(int pi=0; pi < G_CLO_N; pi++)
        {
            auto dyn = new pkpd_dha;
            //fprintf(stderr, "\tlum object created pi = %d \r", pi); fflush(stderr);
            dyn->set_parasitaemia( 20000.0 );    

            dyn->rng = G_RNG;    
            dyn->age = G_CLO_AGE;
            dyn->weight = G_CLO_WEIGHT;
            dyn->initialize_params();                              // NOTE parasitaemia must be set before initializing parameters
            // dyn->vprms[i_lum_pmf] = G_CLO_PMF;
            //fprintf(stderr, "\n\tlum object initialized pi=%d", pi); fflush(stderr);
        
            t0=0.0;
            t1=maximum_enforced_stepsize;
            double stepsize_PMF = pow( G_CLO_PMF, 1.0 / (48.0/maximum_enforced_stepsize) );
        

            //BEGIN - INTEGRATION
            while( t0 < 168.0*4.0 )
            {
                if( dyn->doses_still_remain_to_be_taken )
                {
                    if( dyn->we_are_past_a_dosing_time(t0) )   
                    {
                        dyn->give_next_dose_to_patient(1.0);    // 1.0 means the full dose is given
                                                                // if no dose remains to be given, function does nothing
                    }
                }
                dyn->predict(t0, t1);

                // after integrating the differential equations in the predict function above,
                // we need to ---- GROW THE PARASITES ---- for half-an-hour (i.e the maximum_enforced_stepsize)   
                dyn->y0[ dyn->dim - 1 ] *= stepsize_PMF; 

                t0 += maximum_enforced_stepsize; t1 += maximum_enforced_stepsize;
            }
            //END - INTEGRATION 

            for(int j=0; j<dyn->v_concentration_in_blood.size(); j++ )
            //for(int j=0; j<60; j++ )
            {
                fprintf(stdout, "%d , %10.3f , %10.3f , %10.3f \n", pi, dyn->v_concentration_in_blood_hourtimes[j], dyn->v_concentration_in_blood[j], dyn->v_parasitedensity_in_blood[j] );
            }
        
            delete dyn;
        }

        fprintf(stderr, "\n");
        //fclose(fp3);
    }




    if( G_CLO_THERAPY == therapy_AL )
    {
        //fp3 = fopen("out.lum.allpatients.20240401.csv","w");
        fprintf(stdout, "PID,HOUR,COMP2CONC,PARASITEDENSITY\n" );

        fprintf(stderr, "\n");
        // pi is patient index
        for(int pi=0; pi < G_CLO_N; pi++)
        {
            auto dyn1 = new pkpd_dha();
            auto dyn2 = new pkpd_lum();
            //fprintf(stderr, "\tlum object created pi = %d \r", pi); fflush(stderr);
            dyn1->set_parasitaemia( 20000.0 ); // NOTE: you must set both of these to the same thing    
            dyn2->set_parasitaemia( 20000.0 );    
            
            dyn2->parasites_per_ul_at_first_lum_dose = 20000.0;      // NOTE YOU MUST DO THIS SEPARATELY because the parasitaemia level "at first
                                                                    // lum dose" is a special quantity that affects the lum absorption

            dyn1->rng = G_RNG;    
            dyn1->age = G_CLO_AGE;
            dyn1->weight = G_CLO_WEIGHT;
            dyn1->initialize_params();                      // TODO: change the name of this function to initialize
            dyn2->rng = G_RNG;    
            dyn2->age = G_CLO_AGE;
            dyn2->weight = G_CLO_WEIGHT;
            dyn2->initialize();                             // NB: parasitaemia must be set before initializing parameters
            
            
        
            t0=0.0;
            t1=maximum_enforced_stepsize;           // normally set to 0.5 hours

            // this is the PMF adjusted to the 30-minute stepsize 
            double stepsize_PMF = pow( G_CLO_PMF, 1.0 / (48.0/maximum_enforced_stepsize) );
        

            //BEGIN - INTEGRATION
            while( t0 < 168.0*4.0 )
            {

                // ---- first, calculate artemisinin clearance and killing over a 30-minute period (maximum_enforced_stepsize)
                if( dyn1->doses_still_remain_to_be_taken )
                {
                    if( dyn1->we_are_past_a_dosing_time(t0) )   
                    {
                        dyn1->give_next_dose_to_patient(1.0);    // 1.0 means the full dose is given
                                                                // if no dose remains to be given, function does nothing
                    }
                }
                dyn1->predict(t0, t1);

                // now that we have killed some parasites with the art component of the therapy, we need to
                // adjust the parasite density in the lum object (dyn2) so that it matches the parasite density
                // in the art object (dyn1)
                dyn2->y0[ dyn2->dim - 1 ] = dyn1->y0[ dyn1->dim - 1 ];


                // ---- then, calculate lumefantrine clearance and killing over a 30-minute period (maximum_enforced_stepsize)
                if( dyn2->doses_still_remain_to_be_taken )
                {
                    if( dyn2->we_are_past_a_dosing_time(t0) )   
                    {
                        dyn2->give_next_dose_to_patient(1.0);    // 1.0 means the full dose is given
                                                                // if no dose remains to be given, function does nothing
                    }
                }
                dyn2->predict(t0, t1);

                // now that we have killed some parasites with the lum component of the therapy, we need to
                // adjust the parasite density in the art object (dyn1) so that it matches the parasite density
                // in the lum object (dyn2)
                dyn1->y0[ dyn1->dim - 1 ] = dyn2->y0[ dyn2->dim - 1 ];

                // after integrating the differential equations in the predict functions above,
                // we need to ---- GROW THE PARASITES ---- for half-an-hour (i.e the maximum_enforced_stepsize)   
                double dd_PMF = stepsize_PMF * ( 1.0 / ( 1.0 + ( dyn1->y0[ dyn1->dim - 1 ] / G_DENSITY_50 ) ) );
                dyn1->y0[ dyn1->dim - 1 ] *= dd_PMF; 
                dyn2->y0[ dyn2->dim - 1 ] *= dd_PMF; 


                t0 += maximum_enforced_stepsize; t1 += maximum_enforced_stepsize;



            }
            //END - INTEGRATION 

            for(int j=0; j<dyn1->v_concentration_in_blood.size(); j++ )
            {
                fprintf(stdout, "%d , %10.3f , %10.3f , %10.3f , %10.3f \n", pi, dyn1->v_concentration_in_blood_hourtimes[j], dyn1->v_concentration_in_blood[j], dyn2->v_concentration_in_blood[j], dyn1->v_parasitedensity_in_blood[j] );
            }
        
            delete dyn1;
            delete dyn2;
        }

        fprintf(stderr, "\n");
        //fclose(fp3);
    }


    /* fprintf(stderr, "\n\n");
    for(int pi=0; pi < 10; pi++)
    {
        auto dyn = new pkpd_adq();
        //fprintf(stderr, "\tlum object created pi = %d \r", pi); fflush(stderr);
        dyn->set_parasitaemia( 1000.0 );    

        dyn->rng = G_RNG;    
        dyn->age = 25.0;
        dyn->initialize();                              // NB: parasitaemia must be set before initializing parameters
        fprintf(stderr, "\tADQ object initialized for patiend id, pi=%4d   \r", pi); fflush(stderr);
        
        t0=0.0;
        t1=maximum_enforced_stepsize;
        

        //BEGIN - INTEGRATION
        while( t0 < 168.0*4.0 )
        {
            if( dyn->we_are_past_a_dosing_time(t0) ) // not very efficient ,  need to change this up a bit
            {
                dyn->give_next_dose_to_patient(1.0);    // 1.0 means the full dose is given
                                                        // if no dose remains to be given, function does nothing
            }

            dyn->predict(t0, t1);
            t0 += maximum_enforced_stepsize; t1 += maximum_enforced_stepsize;
        }
        //END - INTEGRATION 

        for(int j=0; j<dyn->v_concentration_in_blood.size(); j++ )
        {
            fprintf(fp4, "%d , %10.3f , %10.3f , %10.3f \n", pi, dyn->v_concentration_in_blood_hourtimes[j], dyn->v_concentration_in_blood[j], dyn->v_parasitedensity_in_blood[j] );
        }
        
        delete dyn;
    }
    fprintf(stderr, "\n\n"); */



    
    /*for(int pi=0; pi<100; pi++)
    {
        auto dyn = new pkpd_dha();
        dyn->rng = G_RNG;    
        dyn->age = 25.0;
        dyn->initialize_params();
        t0=0.0;
        t1=maximum_enforced_stepsize;
        
        while( t0 < 168.0 )
        {
            dyn->predict(t0, t1);
            t0+=maximum_enforced_stepsize; t1+=maximum_enforced_stepsize;
        }

        for(int j=0; j<dyn->v_concentration_in_blood.size(); j++ )
        {
            fprintf(fp2, "%d , %10.3f , %10.3f , %10.3f \n", pi, dyn->v_concentration_in_blood_hourtimes[j], dyn->v_concentration_in_blood[j], dyn->v_parasitedensity_in_blood[j] );
        }
        
        delete dyn;
    }*/

    //for(int k=0; k<num_params; k++) printf("\n\t %3.5f",dyn_ppq->v_individual_pk_prms[k]);
    
    /*fclose(fp);
    fclose(fp2);
    
    fclose(fp4);*/
    
    

    // free memory
    // delete dyn_ppq;
    // delete dyn_dha;
    

    gsl_rng_free( G_RNG );
        
    return 0;
}


// parses command line arguments
void ParseArgs(int argc, char **argv)
{
    int i, start;
    start=1;
    
    string str;
    for(i=start; i<argc; i++)
    {
        str = argv[i];
             if( str == "--lum" )		                G_CLO_THERAPY  		                = therapy_lumefantrine;
	    else if( str == "--adq" ) 		                G_CLO_THERAPY  		                = therapy_amodiaquine;
        else if( str == "--AL" ) 		                G_CLO_THERAPY  		                = therapy_AL;
        else if( str == "--art" ) 		                G_CLO_THERAPY  		                = therapy_artemisinin;
        else if( str == "--age" ) 		                G_CLO_AGE  		                    = atof( argv[++i] );
        else if( str == "--weight" ) 	                G_CLO_WEIGHT	                    = atof( argv[++i] );
        else if( str == "-n" ) 	                        G_CLO_N	                            = atoi( argv[++i] );
        else if( str == "--pmf" ) 	                    G_CLO_PMF	                        = atof( argv[++i] );
        
        //else if( str == "--hill_art" ) 		        G_CLO_HILL_COEFF	                = atoi( argv[++i] );
        //else if( str == "--hill_lum" ) 		        G_CLO_HILL_COEFF	                = atoi( argv[++i] );
        //else if( str == "--ec50_art" ) 		        G_CLO_EC50	                        = atoi( argv[++i] );
        //else if( str == "--ec50_lum" ) 		        G_CLO_EC50	                        = atoi( argv[++i] );
        //else if( str == "--pmax_art" ) 		        G_CLO_PMAX	                        = atoi( argv[++i] );
        //else if (str == "--pmax_lum" ) 		        G_CLO_PMAX	                        = atoi( argv[++i] );
        
        /*else if( str == "--endttr" ) 		            G_CLO_END_TITER	                    = atof( argv[++i] );
        else if( str == "--chainlength" )                G_CLO_CHAIN_LENGTH	                = atoi( argv[++i] );
        else if( str == "--showevery" ) 		            G_CLO_SHOW_EVERY	                = atoi( argv[++i] );
        else if( str == "--burnin" ) 		            G_CLO_BURNIN	                    = atoi( argv[++i] );
        else if( str == "--iltau" ) 		                G_CLO_INTEGRATION_LIMIT_FOR_TAU	    = atof( argv[++i] );
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



