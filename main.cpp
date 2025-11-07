#include "globals.h"
#include "general_functions.h"
#include "pkpd_simulation_functions.h"

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


int main(int argc, char* argv[])
{
    // Create an instance of general_functions to load general functions
    general_functions general_func;

    // Parse command line arguments
    general_func.ParseArgs(argc,argv);

    // Create an instance of pkpd_simulation to load pkpd_simulation functions
    pkpd_simulation pkpd_sim;

    // get random number seed from current time
    auto now = std::chrono::high_resolution_clock::now();
    auto milliseconds =  std::chrono::duration_cast<std::chrono::microseconds>(now.time_since_epoch());
    auto seed = static_cast<unsigned long>(milliseconds.count());
    
    // make random number generator (RNG) the Mersenne Twister which has period 2^19937 - 1
    const gsl_rng_type *TT_RAND = gsl_rng_mt19937;
    G_RNG = gsl_rng_alloc(TT_RAND);
    gsl_rng_set( G_RNG, seed ); // seed the RNG   
    
    // Seeding the GLOBAL RNG using a C++ friendly version from <random>
    G_RNG_CPP.seed(seed);

    // these will always be stochastic unless we are debugging something    
    pkpd_artemether::stochastic = true;
    pkpd_artesunate::stochastic = true;
    pkpd_dha::stochastic = true;

    pkpd_adq::stochastic = true;
    pkpd_lum::stochastic = true;
    pkpd_ppq::stochastic = true;
    
    if(G_CLO_THERAPY == therapy_AL){
        pkpd_sim.simulate_ACT_PKPD_AL();
    }
    else if(G_CLO_THERAPY == therapy_DHA_PPQ){
        pkpd_sim.simulate_ACT_PKPD_DHA_PPQ();
    }
    else if(G_CLO_THERAPY == therapy_AS_AQ)
    {
        pkpd_sim.simulate_ACT_PKPD_AS_AQ();
    }
    else if(G_CLO_THERAPY == therapy_dha ){
        pkpd_sim.simulate_ACT_PKPD_DHA();
    }
    else if(G_CLO_THERAPY == therapy_lumefantrine )
    {
        pkpd_sim.simulate_ACT_PKPD_LUM();
    }
    
    gsl_rng_free( G_RNG );
    return 0;
}






