#include "pkpd_ppq.h"
#include "pkpd_dha.h"

#include <gsl/gsl_rng.h> // random number generators from Gnu Scientific Library
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <chrono>

// GLOBAL RANDOM NUMBER GENERATOR
gsl_rng *G_RNG;		


int main(int argc, char* argv[])
{
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
    
    pkpd_dha::stochastic = false;
    pkpd_ppq::stochastic = true;
    
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

    
    
    
    FILE *fp, *fp2;
    fp = fopen("out.ppq.allpatients.v20191109.csv","w");
    fp2 = fopen("out.dha.allpatients.csv","w");
    fprintf(fp, "PID , HOUR, COMP2CONC , PARASITEDENSITY\n" );
    fprintf(fp2, "PID , HOUR, COMP2CONC , PARASITEDENSITY\n" );
    
    for(int pi=0; pi < 5; pi++)
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
    }

    
    for(int pi=0; pi<10; pi++)
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
    }

    //for(int k=0; k<num_params; k++) printf("\n\t %3.5f",dyn_ppq->v_individual_pk_prms[k]);
    
    fclose(fp);
    fclose(fp2);
    
    

    // free memory
    delete dyn_ppq;
    delete dyn_dha;
    gsl_rng_free( G_RNG );
    


    
    return 0;
}



