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

double G_CLO_AGE = 25.0;
double G_CLO_WEIGHT = 54.0;

int G_CLO_N = 1; // this is the number of patients


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


    
    auto dyn_ppq = new pkpd_ppq();
    auto dyn_dha = new pkpd_dha();
    //pkpd_ppq* dyn_ppq = new pkpd_ppq();   // the above syntax is the more modern way 
    
    pkpd_dha::stochastic = true;
    pkpd_ppq::stochastic = true;
    pkpd_lum::stochastic = true;
    pkpd_adq::stochastic = true;
    
    dyn_ppq->rng = G_RNG;    
    dyn_ppq->age = 25.0;
    dyn_ppq->initialize();
    dyn_dha->rng = G_RNG;    
    dyn_dha->age = 25.0;
    dyn_dha->initialize_params();

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


    if( G_CLO_LUM )
    {
        //fp3 = fopen("out.lum.allpatients.20240401.csv","w");
        fprintf(stdout, "PID,HOUR,COMP2CONC,PARASITEDENSITY\n" );

        fprintf(stderr, "\n");
        // pi is patient index
        for(int pi=0; pi < G_CLO_N; pi++)
        {
            auto dyn = new pkpd_lum();
            //fprintf(stderr, "\tlum object created pi = %d \r", pi); fflush(stderr);
            dyn->set_parasitaemia( 1000.0 );    
            dyn->parasites_per_ul_at_first_lum_dose = 1000.0;   // YOU MUST DO THIS SEPARATELY because the parasitaemia level "at first
                                                            // lum dose" is a special quantity that affects the lum absorption

            dyn->rng = G_RNG;    
            dyn->age = G_CLO_AGE;
            dyn->weight = G_CLO_WEIGHT;
            dyn->initialize();                              // NB: parasitaemia must be set before initializing parameters
            //fprintf(stderr, "\n\tlum object initialized pi=%d", pi); fflush(stderr);
        
            t0=0.0;
            t1=maximum_enforced_stepsize;
        

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
                t0 += maximum_enforced_stepsize; t1 += maximum_enforced_stepsize;
            }
            //END - INTEGRATION 

            for(int j=0; j<dyn->v_concentration_in_blood.size(); j++ )
            {
                fprintf(stdout, "%d , %10.3f , %10.3f , %10.3f \n", pi, dyn->v_concentration_in_blood_hourtimes[j], dyn->v_concentration_in_blood[j], dyn->v_parasitedensity_in_blood[j] );
            }
        
            delete dyn;
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
    delete dyn_ppq;
    delete dyn_dha;
    

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
             if( str == "-lum" )		G_CLO_LUM  		= true;
	    else if( str == "-adq" ) 		G_CLO_ADQ  		= true;
        else if( str == "-age" ) 		G_CLO_AGE  		= atof( argv[++i] );
        else if( str == "-weight" ) 	G_CLO_WEIGHT	= atof( argv[++i] );
        else if( str == "-n" ) 	        G_CLO_N	        = atoi( argv[++i] );
        /*else if( str == "-endttr" ) 		G_CLO_END_TITER	        = atof( argv[++i] );
        else if( str == "-chainlength" )        G_CLO_CHAIN_LENGTH	= atoi( argv[++i] );
        else if( str == "-showevery" ) 		G_CLO_SHOW_EVERY	= atoi( argv[++i] );
        else if( str == "-burnin" ) 		G_CLO_BURNIN	        = atoi( argv[++i] );
        else if( str == "-iltau" ) 		G_CLO_INTEGRATION_LIMIT_FOR_TAU	= atof( argv[++i] );
        else if( str == "-profile" ) 		G_CLO_PROFILE  		= true;
        else if( str == "-censored" ) 		G_CLO_CENSORED 		= true;*/
	else
        {
            fprintf(stderr, "\n\tUnknown option [%s] on command line.\n\n", argv[i]);
            exit(-1);
        }
    }

    return;
}



