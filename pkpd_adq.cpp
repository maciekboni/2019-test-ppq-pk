// Model and parameters adapted from Ali et al., 2018

#include "assert.h"
#include "pkpd_adq.h"

#include <iostream>
#include <string>

bool pkpd_adq::stochastic = true;

// constructor
pkpd_adq::pkpd_adq()
{
    
    vprms.insert( vprms.begin(), adq_num_params, 0.0 );
    assert( vprms.size()==adq_num_params );
    
    // this is the dimensionality of the ODE system
    dim = 9; // eight pk equations and one pd clearance equation

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
    
    patient_id = 0;
    patient_weight = -1.0;
    patient_age = 25.0;
    pregnant = false;

    num_doses_given = 0;
    num_hours_logged = 0;    
    doses_still_remain_to_be_taken = true;
    total_mg_dose_per_occassion = -99.0;    // meaning it is not set yet
    
    patient_age = 25.0;
    patient_blood_volume = 5500000.0; // 5.5L of blood for an adult individual
    
    // 
    pdparam_n = 15.0;
    pdparam_EC50 = exp( 0.525 * log(2700));
    pdparam_Pmax = 0.9; // here you want to enter the max daily killing rate; it will be converted to hourly later


    rng=NULL;
}

// destructor
pkpd_adq::~pkpd_adq()
{
    delete[] y0;
    gsl_odeiv_evolve_free(oe);
    gsl_odeiv_control_free(oc);
    gsl_odeiv_step_free(os);
}

void pkpd_adq::set_parasitaemia( double parasites_per_ul )
{
    y0[dim-1] = parasites_per_ul; //  the final ODE equation is always the Pf asexual parasitaemia
}


// NOTE :: MUST REMEMBER THAT THE FUNCTION BELOW IS A STATIC MEMBER FUNCTION OF THIS CLASS
//
// i think "rhs_ode" must be defined with these four arguments for it to work in
// in the GSL ODE routines (i.e. in practice, this is what works)
//
int pkpd_adq::rhs_ode(double t, const double y[], double f[], void *pkd_object )
{
    pkpd_adq* p = (pkpd_adq*) pkd_object;

    // these are the right-hand sides of the derivatives of the six compartments
    
    // this is compartment 1, the fixed dose compartment, i.e. the hypothetical compartment
    // where the drug goes in first, ka is the transition from this fixed dose compartment to the transition compartments
    f[0] =  - p->vprms[i_adq_k17] * y[0];

    // this is compartment 2, the central compartment, i.e. the blood, for AQ only
    f[1] =  y[7]*p->vprms[i_adq_k82]  +  y[2]*p->vprms[i_adq_k32]  -  y[1]*p->vprms[i_adq_k23]  -  y[1]*p->vprms[i_adq_k24];
    
    // this is compartment 3, the peripheral compartment for AQ in the blood
    f[2] = y[1]*p->vprms[i_adq_k23]  -  y[2]*p->vprms[i_adq_k32];

    // this is compartment 4, the central compartment for DEAQ (a metabolite of AQ) in the blood
    f[3] = p->vprms[i_adq_CF]*y[1]*p->vprms[i_adq_k24]  +  y[4]*p->vprms[i_adq_k54]  +  y[5]*p->vprms[i_adq_k64]  -  y[3]*p->vprms[i_adq_k45]  -  y[3]*p->vprms[i_adq_k46]  -  y[3]*p->vprms[i_adq_k40];

    // this is compartment 5, a peripheral compartment for DEAQ in the blood
    f[4] = y[3]*p->vprms[i_adq_k45]  -  y[4]*p->vprms[i_adq_k54];

    // this is compartment 6, a peripheral compartment for DEAQ in the blood
    f[5] = y[3]*p->vprms[i_adq_k46]  -  y[5]*p->vprms[i_adq_k64];

    // this is compartment 7, a transit compartment from dose to central
    f[6] = y[0]*p->vprms[i_adq_k17]  -  y[6]*p->vprms[i_adq_k78];

    // this is compartment 8, a transit compartment from dose to central
    f[7] = y[6]*p->vprms[i_adq_k78]  -  y[7]*p->vprms[i_adq_k82];
    
    // this is the per/ul parasite population size
    double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow(y[1],p->pdparam_n) / (pow(y[1],p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)) );
    f[8] = -a * y[8];
    
    
    return GSL_SUCCESS;
}


void pkpd_adq::give_next_dose_to_patient( double fractional_dose_taken )
{
    if( num_doses_given >= v_dosing_amounts.size() )
    {
        return;
    }
    else
    {
        // if there were any substantial inter-occassion variability, you would call this function to redraw params
        // redraw_params_before_newdose(); 

        // there is just one minor non-random inter-occassion variability effect that we model here, see below

        if( num_doses_given == 0 )
        {
            y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_adq_F1_indiv_first_dose] * fractional_dose_taken;
        } 
        else
        {
            y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_adq_F1_indiv_later_dose] * fractional_dose_taken;
        }
                
        num_doses_given++;
    }

}


void pkpd_adq::predict( double t0, double t1 )
{
    //int i;
    gsl_odeiv_system sys = {pkpd_adq::rhs_ode, pkpd_adq::jac, dim, this};   // the fourth argument is a void pointer that you are supposed 
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
            v_concentration_in_blood.push_back( y0[1] );                // TODO: you need to push_back 1000 * y0[1] / V
                                                                        // and the result is microgram per liter = ng/mL
            v_concentration_in_blood_metabolite.push_back( y0[3] );     // TODO: you need to push_back 1000 * y0[3] / VC_DEAQ
                                                                        // and the result is microgram per liter = ng/mL
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


void pkpd_adq::initialize_pkpd_object( void )
{
    
    //-- WARNING - THE AGE MEMBER VARIABLE MUST BE SET BEFORE YOU CALL THIS FUNCTION
    
    generate_recommended_dosing_schedule();
    initialize_params();
    
}


void pkpd_adq::initialize_params( void )
{

    // this is the median weight of a patient that these estimates were calibrated for
    double population_median_weight = 42.0;
    
     // all 18 below are point estimates
    double THETA1 = 2960.0;
    double THETA2 = 13500.0;
    double THETA3 = 2310.0;
    double THETA4 = 22700.0;
    double THETA5 = 0.589;
    double THETA6 = 1.0;
    double THETA7 = 258.0;
    double THETA8 = 32.6;
    double THETA9  = 154.0;
    double THETA10 = 2460.0;
    double THETA11 = 5580.0;
    double THETA12 = 31.3;
    double THETA13 = 0.236;
    double THETA14 = 0.224;
    double THETA15 = 11.8;
    double THETA16 = 3.6;
    double THETA17 = 12.9;
    double THETA18 = 3.22;


    //
    // the variables below are random draws from normal distributions; they will be exponentiated to have the final 
    // random variate '_rv' be log-normally distributed
    //
    double ETA1_rv = 0.0;   
    double ETA2_rv = 0.0;
    double ETA3_rv = 0.0;   
    double ETA4_rv = 0.0;   
    double ETA5_rv = 0.0;   
    double ETA6_rv = 0.0;
    double ETA7_rv = 0.0;
    double ETA8_rv = 0.0;
    double ETA9_rv = 0.0;
    double ETA10_rv = 0.0;
    double ETA11_rv = 0.0;
    double ETA12_rv = 0.0;
    double ETA13_rv = 0.0;

    if( pkpd_adq::stochastic )
    {
        ETA1_rv = gsl_ran_gaussian( rng, sqrt(0.103684) );  //  "_rv" means random variate 
        ETA2_rv = gsl_ran_gaussian( rng, sqrt(0.281961) );     
        //ETA3_rv = gsl_ran_gaussian( rng, sqrt(0.0) );       // fixed
        //ETA4_rv = gsl_ran_gaussian( rng, sqrt(0.0) );       // fixed
        ETA5_rv = gsl_ran_gaussian( rng, sqrt(0.616225) );                                                                 
        ETA6_rv = gsl_ran_gaussian( rng, sqrt(0.095481) );
        ETA7_rv = gsl_ran_gaussian( rng, sqrt(0.451584) );
        ETA8_rv = gsl_ran_gaussian( rng, sqrt(0.004) );
        // ETA 9, 10, 11, 12 are fixed
        ETA13_rv = gsl_ran_gaussian( rng, sqrt(0.872356) );
    }


    // in the parameter calculations below
    // "TV" means typical value or population mean for some parameter

    // this block of code if just for slower clearance in neonates, infants, and small toddlers
    // all of this saturates to 1.0 for people 2yo and older
    double PMA = (patient_age * 12.0) + 9.0; // this is some strange age adjustment
    double ME_AQ   = pow( PMA, THETA16 ) / ( pow( PMA, THETA16 ) + pow( THETA15, THETA16 ) );
    double ME_DEAQ = pow( PMA, THETA18 ) / ( pow( PMA, THETA18 ) + pow( THETA17, THETA18 ) );

    double TVCL = THETA1 * pow(patient_weight/50.0 , 0.75 ) * ME_AQ;
    double CL   = TVCL * exp(ETA1_rv);

    double TVV_AQ = THETA2 * pow(patient_weight/50.0 , 1.0 );
    double V_AQ   = TVV_AQ * exp( ETA2_rv );

    double TVQ = THETA3 * pow(patient_weight/50.0 , 0.75 );
    double Q   = TVQ * exp( ETA3_rv );

    double TVVP_AQ = THETA4 * pow(patient_weight/50.0 , 1.0 );
    double VP_AQ   = TVVP_AQ * exp( ETA4_rv );

    double TVKA = THETA5;
    double KA   = TVKA * exp( ETA5_rv );

    // ----BLOCK1 COMMENTED OUT FOR EFFICIENCY----
    // double TCOV          =  1.0;
    // double F1_D          = (1.0 - THETA14*TCOV);
    // double TVF1          = THETA6 * F1_D;
    // double F1            = TVF1 * exp( ETA6_rv );
    // double F1_later_dose = THETA6 * exp( ETA6_rv );
    // ----BLOCK1 COMMENTED OUT FOR EFFICIENCY----
    double F1_later_dose = exp( ETA6_rv );
    double F1            = (1.0 - THETA14) * F1_later_dose;

    double TVVC_DEAQ = THETA7 * pow(patient_weight/50.0 , 1.0 );
    double VC_DEAQ   = TVVC_DEAQ * exp( ETA7_rv );

    double TVCL_DEAQ = THETA8 * pow(patient_weight/50.0 , 0.75 ) * ME_AQ;
    double CL_DEAQ   = TVCL_DEAQ * exp( ETA8_rv );

    double TVQ2 = THETA9 * pow(patient_weight/50.0 , 0.75 );
    double Q2   = TVQ2 * exp( ETA9_rv );

    double TVVP2 = THETA10 * pow(patient_weight/50.0 , 1.0 );
    double VP2   = TVVP2 * exp( ETA10_rv );

    double TVVP3 = THETA11 * pow(patient_weight/50.0 , 1.0 );
    double VP3   = TVVP3 * exp( ETA11_rv );

    double TVQ3 = THETA12 * pow(patient_weight/50.0 , 0.75 );
    double Q3   = TVQ3 * exp( ETA12_rv );

    double TVMT = THETA13;
    double MT   = TVMT * exp( ETA13_rv );

    // simple re-scaling from mg/L to ng/mL (or perhaps the other way around)
    double S2 = V_AQ / 1000.0;
    double S4 = VC_DEAQ / 1000.0;

    double NN = 2.0;
    double KTR = (NN+1.0)/MT; // this is the transition param btw dosing and central compartments

    // ----BLOCK2 COMMENTED OUT FOR EFFICIENCY----
    // double MW_AQ   =  355.87                             ; Molecular weight AQ (g/mol)
    // double MW_DEAQ =  327.813                            ; Molecular weight DEAQ (g/mol)
    // double CF = MW_DEAQ / MW_AQ                          ; Conversion factor AQ -> DEAQ molecular mass
    // double CF = 0.921159412;
    // ----BLOCK2 COMMENTED OUT FOR EFFICIENCY----


    // these are the three transitions params from dose to central; normally they are all the same
    // so this particular structure is a bit strange; but we keep it for consistency with the original 
    // NONMEM model fit
    vprms[i_adq_k17] = KTR;
    vprms[i_adq_k78] = KTR;
    vprms[i_adq_k82] = KA;

    vprms[i_adq_k23] = Q/V_AQ;
    vprms[i_adq_k32] = Q/VP_AQ;

    vprms[i_adq_k24] = CL/V_AQ;

    vprms[i_adq_k45] = Q2/VC_DEAQ;
    vprms[i_adq_k54] = Q2/VP2;
    vprms[i_adq_k46] = Q3/VC_DEAQ;
    vprms[i_adq_k64] = Q3/VP3;

    vprms[i_adq_k40] = CL_DEAQ/VC_DEAQ;

    vprms[i_adq_CF] = 0.921159412; // this is the commented out value from BLOCK2

    vprms[i_adq_F1_indiv_first_dose] = F1;
    vprms[i_adq_F1_indiv_later_dose] = F1_later_dose; // basically, this is just a rv log-normally distributed around 1.0
    
}


void pkpd_adq::redraw_params_before_newdose()
{
     
    // TODO -- first dose of amodiaquine is absorbed more poorly with 22.4% lower absorbtion
    // this is the THETA14 parameter
    
    // NO INTER-OCCASSION VARIABILITY SO NOTHING FOR THIS FUNCTION TO DO
    
}

bool pkpd_adq::we_are_past_a_dosing_time( double current_time )
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

void pkpd_adq::generate_recommended_dosing_schedule()
{
    
    // CHECK -- that this dosing schedule function indeed corresponds to amodiaquine

    double num_tablets_per_dose = -99.0;
    
    // this is in kilograms
    if( patient_weight < 8.0 )
    {
        num_tablets_per_dose = 1.0;
    }
    else if( patient_weight < 17.0 )
    {
        num_tablets_per_dose = 2.0;
    }
    else if( patient_weight < 35.0 )
    {
        num_tablets_per_dose = 4.0;
    }
    else
    {
        num_tablets_per_dose = 8.0;
    }
   
    // 
    total_mg_dose_per_occassion = num_tablets_per_dose * 67.5;  // this is 67.5mg of amodiaquine, which is the pediatric strength
                                                                // PS - in reality, there are different tablet sizes
    
    v_dosing_times.insert( v_dosing_times.begin(), 3, 0.0 );
    v_dosing_times[1] = 24.0;
    v_dosing_times[2] = 48.0; 

    v_dosing_amounts.insert( v_dosing_amounts.begin(), 3, total_mg_dose_per_occassion );
    
}




