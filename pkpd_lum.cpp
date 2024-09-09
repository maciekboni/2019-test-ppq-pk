//#include <iostream>
//#include <string>
//#include <cstdlib>

#include "assert.h"
#include "pkpd_lum.h"

bool pkpd_lum::stochastic = true;

// constructor
pkpd_lum::pkpd_lum(  )
{
    
    vprms.insert( vprms.begin(), lum_num_params, 0.0 );
    assert( vprms.size()==lum_num_params );
    
    // this is the dimensionality of the ODE system
    dim = 4; // three pk equations and one pd clearance equation

    // this is the main vector of state variables
    y0 = new double[dim];
    for(int i=0; i<dim; i++) y0[i]=0.0;
    // above y0 is initialized to zero, because it's the dosing schedule that creates
    // non-zero initial conditions for the ODEs
    
    // BUT we do need to set the initial parasitaemia at time zero to something positive
    // this should be obtained from the person class
    // the last differential equation is for the parasitaemia; it is a per/ul measure
    // y0[dim-1] = 10000.0;
    set_parasitaemia( 10000.0 ); // simply a wrapper for the statement above

    parasites_per_ul_at_first_lum_dose = 10000.0;   // YOU MUST DO THIS SEPARATELY because the parasitaemia level "at first
                                                    // lum dose" is a special quantity that affects the lum absorption


    const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45;
    os 	= gsl_odeiv_step_alloc(T, dim);
    oc 	= gsl_odeiv_control_y_new (1e-6, 0.0);
    oe 	= gsl_odeiv_evolve_alloc(dim);
    
    patient_weight = -1.0;
    median_weight  =  54.0;     // in kilograms 
    weight = median_weight;     // this is the weight that is actually used in the calculations
    pregnant = false;

    num_doses_given = 0;
    doses_still_remain_to_be_taken = true; 
    num_hours_logged = 0;    
    total_mg_dose_per_occassion = -99.0;    // meaning it is not set yet
    
    age = 25.0;
    patient_blood_volume = 5500000.0;       // 5.5L of blood for an adult individual
    central_volume_of_distribution = -99.0; // meaning it is not set yet

    // the parameters 15, exp( 0.525 * log(2700)), and 0.9 give about a 90% drug efficacy for an initial parasitaemia of 10,000/ul (25yo patient, 54kg)
    pdparam_n = 15.0; // default parameter if CLO is not specified
    pdparam_EC50 = exp( 0.525 * log(2700)); // this is about 63.3, TODO: determine where this came from
                                            // default parameter if CLO is not specified
    pdparam_Pmax = 0.9995; // here you want to enter the max daily killing rate; it will be converted to hourly later
                           // default parameter if CLO is not specified

    rng=NULL;

    //printf("\n\nFinished constructing pkpd_lum, the lumefantrine PKPD object that manages drug clearance and parasite clearance.\n");
}

// destructor
pkpd_lum::~pkpd_lum()
{
    delete[] y0;
    gsl_odeiv_evolve_free(oe);
    gsl_odeiv_control_free(oc);
    gsl_odeiv_step_free(os);
}

void pkpd_lum::set_parasitaemia( double parasites_per_ul )
{
    y0[dim-1] = parasites_per_ul; //  the final ODE equation is always the Pf asexual parasitaemia
}


// MUST REMEMBER THAT THE FUNCTION BELOW IS A STATIC MEMBER FUNCTIONS OF THIS CLASS
// 
// i think "rhs_ode" must be defined with these four arguments for it to work in
// in the GSL ODE routines (i.e. in practice, this is what works)
//
int pkpd_lum::rhs_ode(double t, const double y[], double f[], void *pkd_object )
{
    pkpd_lum* p = (pkpd_lum*) pkd_object;

    // these are the right-hand sides of the derivatives of the six compartments
    
    // this is compartment 1, the fixed dose compartment, i.e. the hypothetical compartment
    // where the drug goes in first, ka is the transition from this fixed dose compartment to the central compartment (blood)
    f[0] =  - p->vprms[i_lum_k12] * y[0];

    // this is compartment 2, the central compartment, i.e. the blood
    f[1] =  y[0]*p->vprms[i_lum_k12]  +  y[2]*p->vprms[i_lum_k32]  -  y[1]*p->vprms[i_lum_k23]  -  y[1]*p->vprms[i_lum_k20];
    
    // this is compartment 3, the only peripheral compartment in this model
    f[2] = y[1]*p->vprms[i_lum_k23]  -  y[2]*p->vprms[i_lum_k32];
    
    // this is the per/ul parasite population size
    double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow(y[1],p->pdparam_n) / (pow(y[1],p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)) );
    f[3] = - a * y[3];          // NOTE there is no parasite growth here because the PMF factor for parasite growth is done
                                // manually in the main diff-eq loop
    //f[3] = 0.02 * y[3]  -  a * y[3];
    
    
    return GSL_SUCCESS;
}


void pkpd_lum::give_next_dose_to_patient( double fractional_dose_taken )
{
    if( doses_still_remain_to_be_taken )
    {
        // redraw_params_before_newdose(); // these are the dose-specific parameters that you're drawing here
        
        // basically, for LUM, we do not redraw, bc there is no inter-occassion variability
        // so this whole function just adds a dose or fractional dose

        // for LUM there is no need to get a new F1 param for each dose
        // y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_lum_F1_thisdose] * fractional_dose_taken;

        y0[0] +=  v_dosing_amounts[num_doses_given] * fractional_dose_taken;
        
        num_doses_given++;

        if( num_doses_given >= v_dosing_amounts.size() ) doses_still_remain_to_be_taken=false;
    }

}


void pkpd_lum::predict( double t0, double t1 )
{
    //int i;
    gsl_odeiv_system sys = {pkpd_lum::rhs_ode, pkpd_lum::jac, dim, this};   // the fourth argument is a void pointer that you are supposed 
                                                                            // to use freely; you normally
                                                                            // use it to access the parameters
    double t = t0;    
    double h = 1e-6;
    //int counter=0;
    
    
    while (t < t1)
    {
        // check if time t is equal to or larger than the next scheduled hour to log
        if( t >= ((double)num_hours_logged)  )
        {
            // this is in ng/ml 
            v_concentration_in_blood.push_back( y0[1] * 1000.0 / central_volume_of_distribution  );                                                      
            v_parasitedensity_in_blood.push_back( y0[dim-1] );
            v_concentration_in_blood_hourtimes.push_back( t );
            
            num_hours_logged++;
        }

        // carry out the Runge-Kutta integration
        int status = gsl_odeiv_evolve_apply(oe, oc, os, &sys, &t, t1, &h, y0);

        if (status != GSL_SUCCESS)
            break;
    }
    
}


void pkpd_lum::initialize( void )
{
    
    //-- WARNING -- the age member variable must be set before you call this function -- add this check
    

    // NOTE must call the two functions below in this order -- dosing schedule needs to be set first
    generate_recommended_dosing_schedule();
    initialize_params();
    
    
}




void pkpd_lum::initialize_params( void )
{
    
     // all 10 below are point estimates
    double THETA1 = 1.35;
    double THETA2 = 11.2;           // NOTE: this is the "central volume of distribution" in liters in a 42kg adult which in non-PKPD language
                                    // is the volume of the blood plus the volume of everything else that is in 
                                    // instantaneous equilibrium with the blood ... you will need this when reporting
                                    // concentration, bc concentration = mg / this volume
    double THETA3 = 0.344;
    double THETA4 = 59.0;
    double THETA5 = 0.0386;
    double THETA6 = 1.0;
    double THETA7 = 3.86;
    double THETA8 = 0.352;
    double THETA9  =  -0.643;
    double THETA10 =  -0.343;


    //
    // the variables below are random draws from normal distributions; they will be exponentiated to have the final 
    // random variate '_rv' be log-normally distributed
    //
    double ETA1_rv = 0.0;   // this is fixed bc there was not enough data to identify/infer variability around this process
    double ETA2_rv = 0.0;
    double ETA3_rv = 0.0;   // this is fixed bc there was not enough data to identify/infer variability around this process
    double ETA4_rv = 0.0;   // this is fixed bc there was not enough data to identify/infer variability around this process
    double ETA5_rv = 0.0;   // this is fixed bc there was not enough data to identify/infer variability around this process
    double ETA6_rv = 0.0;

    if( pkpd_lum::stochastic )
    {
        ETA2_rv = gsl_ran_gaussian( rng, sqrt(1.12) );  
        ETA6_rv = gsl_ran_gaussian( rng, sqrt(0.402) );     //  "_rv" means random variate this represents
                                                            //  between-patient variability in bioavailability
    }






    // in the parameter calculations below
    // "TV" means typical value or population mean for some parameter



    // F1_indiv is the relative absorbtion level for this individual
    
    // before we take into account the effects of dose, the value of TVF1 is one; this is the bioavailability parameter
    //
    // we will need to multiply this by a DS parameter (something to do w mg/kg scaling) and a "PARASITE" parameter that tells
    // us what the reduced bioavailability of lumefantrine is when parasitaemia is high
    

    double BXPAR = -0.343; // parameter from box-cox transformation
    double PHI = exp( ETA6_rv );
    double ETATR = ( pow(PHI, BXPAR) - 1.0  ) / BXPAR ; 

    double D50 = THETA7; // 

    //fprintf(stderr, "\n\tinside initialize_params function -- about to set DS param"); fflush(stderr);
    

    // this you keep fixed, and you use the total mg dose per occassion, and NOT any randomly drawn number
    double DS = 1.0 - ( total_mg_dose_per_occassion/weight ) / ( ( total_mg_dose_per_occassion/weight ) + D50  );


    
    
    // PARASITE = ((LNPC /4.20)**THETA(9)) -- TODO: -- check the log type on the parasitaemia (CONFIRMED on 3/31/2024 that it is log-10)
    //      check if it's parasites/microliter (CONFIRMED also on 3/31/2024)
    //      
    //fprintf(stderr, "\n\tinside initialize_params function -- about to set PARASITE param"); fflush(stderr);
    double PARASITE = pow( log10( parasites_per_ul_at_first_lum_dose ) / 4.20 , THETA9 );
    //fprintf(stderr, "\n\tinside initialize_params function -- PARASITE param set"); fflush(stderr);

    double TVF1 = THETA6 * DS * PARASITE;
    double F1 = TVF1 * exp(ETATR);

    double TVQ = THETA3 * pow( weight/42.0 , 0.75 );  // allometric scaling for weight on the Q parameter
    double Q = TVQ * exp( ETA3_rv );

    double TVV = THETA2 * pow( weight/42.0 , 1.0 ); 
    double V = TVV * exp( ETA2_rv );  
    central_volume_of_distribution = V;

    double TVCL = THETA1 * pow( weight/42.0 , 0.75 );  // allometric scaling for weight on the clearance parameter
    double CL = TVCL * exp( ETA1_rv );

    double TVVP = THETA4 * pow( weight/42.0 , 1.0 );
    double VP = TVVP * exp( ETA4_rv );

    double PREGNANCY = pregnant ? (1.0 + THETA8) : 1.0;
    double TVKA = THETA5 * PREGNANCY;
    double KA = TVKA * exp( ETA5_rv );

    //fprintf(stderr, "\n\tinside initialize_params function -- calculations finished"); fflush(stderr);

    
    vprms[i_lum_k12] = KA;
    vprms[i_lum_k23] = Q/V;
    vprms[i_lum_k32] = Q/VP;
    vprms[i_lum_k20] = CL/V;
    vprms[i_lum_F1_indiv] = F1;
    //vprms[i_lum_pmf] = 10.0; // this will get assigned in the main file; the reason is that we don't know what it is and we need to experiment with it a bit
}


void pkpd_lum::redraw_params_before_newdose()
{
     
    // ### first, you don't receive the full dose.  You may receive 80% or 110% of the dose depending
    // on whether you're sitting or standing, whether you've recently had a big meal, whether some gets stuck
    // between your teeth; below, we set the parameter F1 (with some random draws) to adjust this initial dose
    
    // NO INTER-OCCASSION VARIABILITY SO NOTHING FOR THIS FUNCTION TO DO
    
}

bool pkpd_lum::we_are_past_a_dosing_time( double current_time )
{
    // check if there are still any doses left to give
    if( num_doses_given < v_dosing_times.size() )
    {
        if( current_time >= v_dosing_times[num_doses_given] )
        {
            return true;
        }
    }

    return false;
}

void pkpd_lum::generate_recommended_dosing_schedule()
{

    // DOSING GUIDELINES SAY    0,  8, 24, 36, 48, 60
    // BUT WE CAN JUST DO       0, 12, 24, 36, 48, 60

    
    double num_tablets_per_dose;
    
    if( weight < 5.0 )
    {
        num_tablets_per_dose = 0.0;
    }
    else if( weight < 14.0 )
    {
        num_tablets_per_dose = 1.0;
    }
    else if( weight < 24.0 )
    {
        num_tablets_per_dose = 2.0;
    }
    else if( weight < 34.0 )
    {
        num_tablets_per_dose = 3.0;
    }
    else
    {
        num_tablets_per_dose = 4.0;
    }
   
    // NOTE - do not confuse this with daily dosing - AL is taken twice daily, two occassions per day
    total_mg_dose_per_occassion = num_tablets_per_dose * 120.0;
    
    v_dosing_times.insert( v_dosing_times.begin(), 6, 0.0 );
    v_dosing_times[0] = 0.0;
    v_dosing_times[1] = 12.0;
    v_dosing_times[2] = 24.0; 
    v_dosing_times[3] = 36.0;
    v_dosing_times[4] = 48.0;
    v_dosing_times[5] = 60.0;

    v_dosing_amounts.insert( v_dosing_amounts.begin(), 6, total_mg_dose_per_occassion );
    
}




