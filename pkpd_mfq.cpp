// Model and parameter values from Guidi et al., 2019

//#include <iostream>
//#include <string>
//#include <cstdlib>

#include "assert.h"
#include "pkpd_mfq.h"

bool pkpd_mfq::stochastic = true;

// constructor
pkpd_mfq::pkpd_mfq(  )
{
    
    vprms.insert( vprms.begin(), mfq_num_params, 0.0 );
    assert( vprms.size()==mfq_num_params );
    
    // this is the dimensionality of the ODE system
    dim = 4; // three pk equations and one pd clearance equation

    // this is the main vector of state variables
    y0 = new double[dim];
    for(int i=0; i<dim; i++) y0[i]=0.0;
    // above y0 is initialized to zero, because it's the dosing schedule that creates
    // non-zero initial conditions for the ODEs
    
    // BUT we do need to set the initial parasitaemia at time zero to something positive
    // this should be obtained from the person class
    set_parasitaemia( 10000.0 ); // simply a wrapper for the statement above

    const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45;
    os 	= gsl_odeiv_step_alloc(T, dim);
    oc 	= gsl_odeiv_control_y_new (1e-6, 0.0);
    oe 	= gsl_odeiv_evolve_alloc(dim);
    
    patient_weight = -1.0;
    median_weight  =  54.0;     // in kilograms 
    weight = median_weight;     // this is the weight that is actually used in the calculations

    num_doses_given = 0;
    num_hours_logged = 0;    
    total_mg_dose_per_occassion = -99.0; // meaning it is not set yet
    doses_still_remain_to_be_taken = true;
    
    age = 25.0;
    patient_blood_volume = 5500000.0; // 5.5L of blood for an adult individual
    
    // 
    pdparam_n = 15.0;
    pdparam_EC50 = exp( 0.525 * log(2700));
    pdparam_Pmax = 0.9; // here you want to enter the max daily killing rate; it will be converted to hourly later


    rng=NULL;
}

// destructor
pkpd_mfq::~pkpd_mfq()
{
    delete[] y0;
    gsl_odeiv_evolve_free(oe);
    gsl_odeiv_control_free(oc);
    gsl_odeiv_step_free(os);
}

void pkpd_mfq::set_parasitaemia( double parasites_per_ul )
{
    y0[dim-1] = parasites_per_ul; //  the final ODE equation is always the Pf asexual parasitaemia
}


// NOTE :: MUST REMEMBER THAT THE FUNCTION BELOW IS A STATIC MEMBER FUNCTION OF THIS CLASS
//
// i think "rhs_ode" must be defined with these four arguments for it to work in
// in the GSL ODE routines (i.e. in practice, this is what works)
//
int pkpd_mfq::rhs_ode(double t, const double y[], double f[], void *pkd_object )
{
    pkpd_mfq* p = (pkpd_mfq*) pkd_object;

    // these are the right-hand sides of the derivatives of the six compartments
    
    // this is compartment 1, the fixed dose compartment, i.e. the hypothetical compartment
    // where the drug goes in first, ka is the transition from this fixed dose compartment to the transition compartments
    f[0] =  - p->vprms[i_mfq_k12] * y[0];

    // this is compartment 2, the central compartment, i.e. the blood
    f[1] =  y[0]*p->vprms[i_mfq_k12]  +  y[2]*p->vprms[i_mfq_k32]  -  y[1]*p->vprms[i_mfq_k23]  -  y[1]*p->vprms[i_mfq_k20];
    
    // this is compartment 3, the peripheral compartment for AQ in the blood
    f[2] = y[1]*p->vprms[i_mfq_k23]  -  y[2]*p->vprms[i_mfq_k32];
    
    // this is the per/ul parasite population size
    double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow(y[1],p->pdparam_n) / (pow(y[1],p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)) );
    f[3] = -a * y[3];
    
    
    return GSL_SUCCESS;
}


void pkpd_mfq::give_next_dose_to_patient( double fractional_dose_taken )
{
    if( num_doses_given >= v_dosing_amounts.size() )
    {
        return;
    }
    else
    {
        // if there were any substantial inter-occassion variability, you would call this function to redraw params
        // redraw_params_before_newdose(); 

        // there is just one minor change here : you need to reset the aborption rate at the second dose

        if( num_doses_given == 1 )
        {
            vprms[i_mfq_k12] = vprms[i_mfq_k12_late];   // TODO - this happens regardless of whether this dose is taken or not
                                                        // in other words, you need to change this *at* the 24-hour mark not at 
                                                        // the second dose, which could be skipped
                                                        //
                                                        // 
        } 
                
        num_doses_given++;
    }

}


void pkpd_mfq::predict( double t0, double t1 )
{
    //int i;
    gsl_odeiv_system sys = {pkpd_mfq::rhs_ode, pkpd_mfq::jac, dim, this};   // the fourth argument is a void pointer that you are supposed 
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
            v_concentration_in_blood.push_back( y0[1] );
            v_parasitedensity_in_blood.push_back( y0[dim-1] );
            v_concentration_in_blood_hourtimes.push_back( t );
            
            num_hours_logged++;
            if( num_hours_logged >= 24 ) vprms[i_mfq_k12] = vprms[i_mfq_k12_late];
        }


        // carry out the Runge-Kutta integration
        int status = gsl_odeiv_evolve_apply(oe, oc, os, &sys, &t, t1, &h, y0);

        if (status != GSL_SUCCESS)
            break;
    }
    
}


void pkpd_mfq::initialize( void )
{
    
    //-- WARNING - THE AGE MEMBER VARIABLE MUST BE SET BEFORE YOU CALL THIS FUNCTION
    
    generate_recommended_dosing_schedule();
    initialize_params();
    
    
}




void pkpd_mfq::initialize_params( void )
{
     // all 10 below are point estimates .. they should end with a '_pe' for point estimate
    double THETA1 = 0.45;
    double THETA2 = 95.0;           // NB : this is the "central volume of distribution" in liters IN AN ADULT which in non-PKPD language
                                    // is the volume of the blood plus the volume of everything else that is in 
                                    // instantaneous equilibrium with the blood ... you will need this when reporting
                                    // concentration, bc concentration = mg / this "central volume of distribution"
    double THETA3 = 0.35;
    double THETA4 = 60.0;
    double THETA5 = 0.17;           // this is the early-stage absorption rate
    double THETA6 = 1.0;

    double THETA7 = 0.4;            // this is the late-stage absorption rate

    double THETA8 = -0.67;


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

    if( pkpd_mfq::stochastic )
    {
        ETA1_rv = gsl_ran_gaussian( rng, sqrt(0.1521) );      // need to validate this w original authors

        ETA5_rv = gsl_ran_gaussian( rng, sqrt(0.8281) );      // need to validate this w original authors
        ETA6_rv = gsl_ran_gaussian( rng, sqrt(0.0784) );      // need to validate this w original authors                                                  

        //ETA7_rv = gsl_ran_gaussian( rng, sqrt(0.8281) );    // need to validate this w original authors
    }






    // in the parameter calculations below
    // "TV" means typical value or population mean for some parameter



    // F1_indiv is the relative absorbtion level for this individual
    
    // before we take into account the effects of dose, the value of TVF1 is one; this is the bioavailability paramete    


    double TVF1 = THETA6;
    double F1 = TVF1 * exp(ETA6_rv);

    double TVQ = THETA3 * pow( weight/12.2 , 0.75 );  // allometric scaling for weight on the Q parameter
    double Q = TVQ * exp( ETA3_rv );

    double TVV = THETA2 * pow( weight/12.2 , 1.0 ); 
    double V = TVV * exp( ETA2_rv );  // TODO - make sure you use this V variable as the "individual central volume of distribution" when presenting a drug concentration in the blood 
    
    double TVCL = THETA1 * pow( weight/12.2 , 0.75 );  // allometric scaling for weight on the clearance parameter
    double CL = TVCL * exp( ETA1_rv );

    double TVVP = THETA4 * pow( weight/12.2 , 1.0 );
    double VP = TVVP * exp( ETA4_rv );

    double COV = 1.0 + THETA8*( age - 2.6 ); // definitely MUST BE VALIDATED BY THE AUTHORS

    double TVKA = THETA5 * COV;         // WARNING -- THIS CAN GO NEGATIVE FOR AGES 4 AND HIGHER
    double KA = TVKA * exp( ETA5_rv );
    double KA_EARLY = KA;

    double TVKA_LATE = THETA7 * COV;
    //double KA_LATE = TVKA_LATE * exp( ETA7_rv );
    double KA_LATE = TVKA_LATE;


    vprms[i_mfq_k12] = KA;
    vprms[i_mfq_k12_early] = KA_EARLY;
    vprms[i_mfq_k12_late] = KA_LATE;

    vprms[i_mfq_k23] = Q/V;
    vprms[i_mfq_k32] = Q/VP;
    vprms[i_mfq_k20] = CL/V;
    vprms[i_mfq_F1_indiv_first_dose] = F1;
    //vprms[i_mfq_F1_indiv_later_dose] = F1; // Need to figure out if bioavailability changes after the first dose

}


void pkpd_mfq::redraw_params_before_newdose()
{
     
    // first dose of amodiaquine is absorbed more poorly with 22.4% lower absorbtion
    // this is the THETA14 parameter
    
    // NO INTER-OCCASSION VARIABILITY SO NOTHING FOR THIS FUNCTION TO DO

    // NO NEED TO REDRAW, SIMPLY CHANGE K12 TRANSITION PARAMETER


    
}

bool pkpd_mfq::we_are_past_a_dosing_time( double current_time )
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

void pkpd_mfq::generate_recommended_dosing_schedule()
{
    
    double num_tablets_per_dose;
    
    if( weight < 5.0 )
    {
        num_tablets_per_dose = 0.0;
    }
    else if( weight < 9.0 )
    {
        num_tablets_per_dose = 1.0;
    }
    else if( weight < 18.0 )
    {
        num_tablets_per_dose = 2.0;
    }
    else if( weight < 30.0 )
    {
        num_tablets_per_dose = 4.0;
    }
    else
    {
        num_tablets_per_dose = 8.0;
    }
   
    // 
    double total_mg_dose = num_tablets_per_dose * 50.0; // this is 50mg of mefloquine, which is 55mg of MQ-hydrochloride
    
    v_dosing_times.insert( v_dosing_times.begin(), 3, 0.0 );
    v_dosing_times[1] = 24.0;
    v_dosing_times[2] = 48.0; 

    v_dosing_amounts.insert( v_dosing_amounts.begin(), 3, total_mg_dose );
    
}




