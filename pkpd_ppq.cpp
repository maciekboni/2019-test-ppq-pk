//#include <iostream>
//#include <string>
//#include <cstdlib>

#include "assert.h"
#include "pkpd_ppq.h"

bool pkpd_ppq::stochastic = true;

// constructor
pkpd_ppq::pkpd_ppq(  )
{
    
    vprms.insert( vprms.begin(), ppq_num_params, 0.0 );
    assert( vprms.size()==ppq_num_params );
    
    // this is the dimensionality of the ODE system
    dim = 7;

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

    const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45;
    os 	= gsl_odeiv_step_alloc(T, dim);
    oc 	= gsl_odeiv_control_y_new (1e-6, 0.0);
    oe 	= gsl_odeiv_evolve_alloc(dim);
    
    patient_weight = -1.0;
    median_weight  =  54.0;
    weight = median_weight;  // this is the weight that is actually used in the calculations

    num_doses_given = 0;
    num_hours_logged = 0;    
    
    age = 25.0;
    patient_blood_volume = 5500000.0; // 5.5L of blood for an adult individual
    
    // the parameters 15, exp( 0.525 * log(2700)), and 0.9 give about a 90% drug efficacy for an initial parasitaemia of 10,000/ul (25yo patient, 54kg)
    pdparam_n = 15.0;
    pdparam_EC50 = exp( 0.525 * log(2700));
    pdparam_Pmax = 0.9; // here you want to enter the max daily killing rate; it will be converted to hourly later


    rng=NULL;
}

// destructor
pkpd_ppq::~pkpd_ppq()
{
    delete[] y0;
    gsl_odeiv_evolve_free(oe);
    gsl_odeiv_control_free(oc);
    gsl_odeiv_step_free(os);
}

void pkpd_ppq::set_parasitaemia( double parasites_per_ul )
{
    y0[dim-1] = parasites_per_ul;
}


// NOTE :: MUST REMEMBER THAT THE FUNCTION BELOW IS A STATIC MEMBER FUNCTIONS OF THIS CLASS
//
// i think "rhs_ode" must be defined with these four arguments for it to work in
// in the GSL ODE routines
//
int pkpd_ppq::rhs_ode(double t, const double y[], double f[], void *pkd_object )
{
    pkpd_ppq* p = (pkpd_ppq*) pkd_object;

    // these are the right-hand sides of the derivatives of the six compartments
    
    // this is compartment 1, the fixed dose compartment, i.e. the hypothetical compartment
    // where the drug goes in first
    f[0] =  - p->vprms[i_ppq_k15] * y[0];

    // this is compartment 2, the central compartment, i.e. the blood
    f[1] =  y[5]*p->vprms[i_ppq_k62] - y[1]*p->vprms[i_ppq_k20] - y[1]*p->vprms[i_ppq_k23] + y[2]*p->vprms[i_ppq_k32] - y[1]*p->vprms[i_ppq_k24] + y[3]*p->vprms[i_ppq_k42];
    
    // this is compartment 3, the first peripheral compartment; it is V_P1 in the Hoglund paper
    f[2] = y[1]*p->vprms[i_ppq_k23] - y[2]*p->vprms[i_ppq_k32];

    // this is compartment 4, the second peripheral compartment; it is V_P2 in the Hoglund paper
    f[3] = y[1]*p->vprms[i_ppq_k24] - y[3]*p->vprms[i_ppq_k42];
    
    // this is compartment 5, which is the first transit compartment from the dose compartment to the central compartment
    f[4] = y[0]*p->vprms[i_ppq_k15] - y[4]*p->vprms[i_ppq_k56];

    // this is compartment 6, which is the second transit compartment from the dose compartment to the central compartment
    f[5] = y[4]*p->vprms[i_ppq_k56] - y[5]*p->vprms[i_ppq_k62];
    
    // this is the per/ul parasite population size
    double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow(y[1],p->pdparam_n) / (pow(y[1],p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)) );
    f[6] = -a * y[6];
    
    
    return GSL_SUCCESS;
}


void pkpd_ppq::give_next_dose_to_patient( double fractional_dose_taken )
{
    redraw_params_before_newdose(); // these are the dose-specific parameters that you're drawing here
    y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_ppq_F1_thisdose] * fractional_dose_taken;
    num_doses_given++;
}


void pkpd_ppq::predict( double t0, double t1 )
{
    //int i;
    gsl_odeiv_system sys = {pkpd_ppq::rhs_ode, pkpd_ppq::jac, dim, this};   // the fourth argument is a void pointer that you are supossed to use freely; you normally
                                                                                        // use it to access the paramters
    double t = t0;    
    double h = 1e-6;
    //int counter=0;
    
    
    while (t < t1)
    {
        // check if there are still doses to give
        //if( num_doses_given < v_dosing_times.size() )
        //{
            // check if time t is equal to or larger than the next scheduled dose
            //if( t >= v_dosing_times[num_doses_given]  )
            //{
                //redraw_params_before_newdose();
                
                // add the new dose amount to the "dose compartment", i.e. the first compartment
                //y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_ppq_F1_thisdose];
                
                //num_doses_given++;
            //}
        //}

        // check if time t is equal to or larger than the next scheduled hour to log
        if( t >= ((double)num_hours_logged)  )
        {
            v_concentration_in_blood.push_back( y0[1] );
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


void pkpd_ppq::initialize( void )
{
    
    //WARNING - THE AGE MEMBER VARIABLE MUST BE SET BEFORE YOU CALL THIS FUNCTION
    
    if( stochastic ) 
    {
        initialize_params();
    }
    else
    {
        initialize_params_w_population_means();
    }
    
    generate_recommended_dosing_schedule();
}



void pkpd_ppq::initialize_params_w_population_means( void )
{

    // ### ###  as a safety this should default to one; you multiply the dose amount given (in the dose compartment)
    //          by this factor F1; the dose compartment is compartment 0 here and compartment 1 in the PLoS Med paper
    vprms[i_ppq_F1_indiv] = 1.0;

 
    // ### ### KTR is the transition rate among the first three compartments
    double TVMT_pe = 2.11; // this is the point estimate (_pe) for TVMT; there is no need to draw a random variate here
    double MT = TVMT_pe; // * exp( ETA7_rv );
    double KTR = 3.0/MT;
    vprms[i_ppq_k15] = KTR;
    vprms[i_ppq_k56] = KTR;
    vprms[i_ppq_k62] = KTR;

    
    // ### ### this is the transition rate from the central compt to peripheral compt #1
    double ACL = 0.75;  // this is the allometric scaling parameter for weight's influence on the Q1 parameter
    double AV = 1.0;    // this is the allometric scaling parameter for weight's influence on the V2 parameter

    double THETA2_pe = 2910.0;
    double THETA3_pe = 310.0;
    
    double Q1 = THETA3_pe * pow( weight/median_weight, ACL ); 
    double V2 = THETA2_pe * pow( weight/median_weight, AV ); 

//     double ETA2_rv = gsl_ran_gaussian( rng, sqrt(0.371) );  
//     V2 *= exp( ETA2_rv );

    // NOTE ETA3 is fixed at zero; so we do not draw
    // double ETA3_rv = gsl_ran_gaussian( rng, sqrt(0.0) );
    // Q1 *= exp( ETA3 ); // should be zero
    
    vprms[i_ppq_k23] = Q1/V2;
    
    
    // ### ### this is the transition rate from to peripheral compt #1 back to the central compt 
    //         TVV3 = THETA(4)*(WT/M_WE)**AV;
    //         V3 = TVV3*EXP(ETA(4));
    double THETA4_pe = 4910.0;
    //double ETA4_rv = gsl_ran_gaussian( rng, sqrt(0.0558) ); 
    double TVV3 = THETA4_pe * pow( weight/median_weight, AV ); 
    double V3 = TVV3; // * exp(ETA4_rv);
    
    vprms[i_ppq_k32] = Q1/V3;
    
    
    // ### ###  this is the transition rate from the central compt to peripheral compt #2
    //          Q2 = TVQ2*EXP(ETA(5));
    //          TVQ2 = THETA(5)*(WT/M_WE)**ACL;
    double THETA5_pe = 105.0;
    //double ETA5_rv = gsl_ran_gaussian( rng, sqrt(0.0541) ); 
    double TVQ2 = THETA5_pe * pow( weight/median_weight, ACL ); 
    double Q2 = TVQ2; //* exp(ETA5_rv); 
    vprms[i_ppq_k24] = Q2/V2;

    
    // ### ### this is the transition rate from to peripheral compt #2 back to the central compt 
    //         V4 = TVV4*EXP(ETA(6));
    //         TVV4 = THETA(6)*(WT/M_WE)**AV;
    double THETA6_pe = 30900.0;
    //double ETA6_rv = gsl_ran_gaussian( rng, sqrt(0.114) ); 
    double TVV4 = THETA6_pe * pow( weight/median_weight, AV ); 
    double V4 = TVV4; // * exp(ETA6_rv); 
    vprms[i_ppq_k42] = Q2/V4;
    
    
    // ### ### this is the exit rate from the central compartment (the final exit rate in the model
    //         CL = TVCL*EXP(ETA(1));
    //         TVCL = THETA(1)*MF*(WT/M_WE)**ACL;
    double THETA1_pe = 55.4; 
    //double ETA1_rv = gsl_ran_gaussian( rng, sqrt(0.0752) ); 
    double HILL = 5.51;
    double EM50 = 0.575; 
    double MF = pow(age,HILL) / ( pow(age,HILL) + pow(EM50,HILL) );
    
    double TVCL = THETA1_pe * MF * pow( weight/median_weight, ACL ); 
    double CL = TVCL; //* exp(ETA1_rv); 

    vprms[i_ppq_k20] = CL/V2;    

    
}


void pkpd_ppq::initialize_params( void )
{
    //NOTE --- in this function alone, THE KTR PARAM HERE HAS INTER-PATIENT VARIABILITY BUT NO INTER-DOSE VARIABILITY
    //     --- THE F1 PARAMETER HERE IS JUST DEFUALTED TO ONE
    //     --- BOTH OF THESE WILL BE REDRAWN IN THE "redraw_params_before_newdose" FUNCTION IF IT IS CALLED 
    
 

    // ### ### F1_indiv is the relative absorbtion level for this individual
    
    // before we take into account the effects of dose, the value of TVF1 is one
    double TVF1 = 1.0;
    double ETA8_rv = gsl_ran_gaussian( rng, sqrt(0.158) ); // this represents between-patient variability 
    double F1 = TVF1 * exp(ETA8_rv);  
    vprms[i_ppq_F1_indiv] = F1;
    
    

    
    // ### ### KTR is the transition rate among the first three compartments
    double TVMT_pe = 2.11; // this is the point estimate (_pe) for TVMT; there is no need to draw a random variate here
    double ETA7_rv = gsl_ran_gaussian( rng, sqrt(0.135) );  // _rv means random variate
                                                            // NOTE - the second argument to this function call needs to 
                                                            //        be the STANDARD DEVIATION not the variance
    double MT = TVMT_pe * exp( ETA7_rv );
    
    // NOTE at this point you have an MT value without any effect of dose order (i.e. whether it's dose 1, dose 2, etc.
    //      later, you must draw another mean-zero normal rv, and multiply by the value above
    
    double KTR = 3.0/MT;
    vprms[i_ppq_k15] = KTR;
    vprms[i_ppq_k56] = KTR;
    vprms[i_ppq_k62] = KTR;

    
    // ### ### this is the transition rate from the central compt to peripheral compt #1
    double ACL = 0.75;  // this is the allometric scaling parameter for weight's influence on the Q1 parameter
    double AV = 1.0;    // this is the allometric scaling parameter for weight's influence on the V2 parameter

    double THETA2_pe = 2910.0;
    double THETA3_pe = 310.0;
    
    double Q1 = THETA3_pe * pow( weight/median_weight, ACL ); 
    double V2 = THETA2_pe * pow( weight/median_weight, AV ); 

    double ETA2_rv = gsl_ran_gaussian( rng, sqrt(0.371) );  
    V2 *= exp( ETA2_rv );

    // NOTE ETA3 is fixed at zero; so we do not draw
    // double ETA3_rv = gsl_ran_gaussian( rng, sqrt(0.0) );
    // Q1 *= exp( ETA3 ); // should be zero
    
    vprms[i_ppq_k23] = Q1/V2;
    
    
    // ### ### this is the transition rate from to peripheral compt #1 back to the central compt 
    //         TVV3 = THETA(4)*(WT/M_WE)**AV;
    //         V3 = TVV3*EXP(ETA(4));
    double THETA4_pe = 4910.0;
    double ETA4_rv = gsl_ran_gaussian( rng, sqrt(0.0558) ); 
    double TVV3 = THETA4_pe * pow( weight/median_weight, AV ); 
    double V3 = TVV3 * exp(ETA4_rv);
    
    vprms[i_ppq_k32] = Q1/V3;
    
    
    // ### ###  this is the transition rate from the central compt to peripheral compt #2
    //          Q2 = TVQ2*EXP(ETA(5));
    //          TVQ2 = THETA(5)*(WT/M_WE)**ACL;
    double THETA5_pe = 105.0;
    double ETA5_rv = gsl_ran_gaussian( rng, sqrt(0.0541) ); 
    double TVQ2 = THETA5_pe * pow( weight/median_weight, ACL ); 
    double Q2 = TVQ2 * exp(ETA5_rv); 
    vprms[i_ppq_k24] = Q2/V2;

    
    // ### ### this is the transition rate from to peripheral compt #2 back to the central compt 
    //         V4 = TVV4*EXP(ETA(6));
    //         TVV4 = THETA(6)*(WT/M_WE)**AV;
    double THETA6_pe = 30900.0;
    double ETA6_rv = gsl_ran_gaussian( rng, sqrt(0.114) ); 
    double TVV4 = THETA6_pe * pow( weight/median_weight, AV ); 
    double V4 = TVV4 * exp(ETA6_rv); 
    vprms[i_ppq_k42] = Q2/V4;
    
    
    // ### ### this is the exit rate from the central compartment (the final exit rate in the model
    //         CL = TVCL*EXP(ETA(1));
    //         TVCL = THETA(1)*MF*(WT/M_WE)**ACL;
    double THETA1_pe = 55.4; 
    double ETA1_rv = gsl_ran_gaussian( rng, sqrt(0.0752) ); 
    double HILL = 5.51;
    double EM50 = 0.575; 
    double MF = pow(age,HILL) / ( pow(age,HILL) + pow(EM50,HILL) );
    
    double TVCL = THETA1_pe * MF * pow( weight/median_weight, ACL ); 
    double CL = TVCL * exp(ETA1_rv); 

    vprms[i_ppq_k20] = CL/V2;    
    
}


void pkpd_ppq::redraw_params_before_newdose()
{
     
    // ### first, you don't receive the full dose.  You may receive 80% or 110% of the dose depending
    // on whether you're sitting or standing, whether you've recently had a big meal, whether some gets stuck
    // between your teeth; below, we set the parameter F1 (with some random draws) to adjust this initial dose
    
    double IOV_rv=0.0;
    if(pkpd_ppq::stochastic)
    {
        IOV_rv = gsl_ran_gaussian( rng, sqrt(0.252) );  
        //IOV_rv = gsl_ran_gaussian( rng, sqrt(0.000001) );
    }
    
    // this is the mean relative increase in bioavailability(?) from dose to dose
    double F1D_pe = 0.236; 

    // this is the "dose occassion", i.e. the order of the dose (first, second, etc)
    double OCC = 1.0 + (double)num_doses_given; // NOTE the RHS here is a class member
    
    double F1COVD = (1.0 + F1D_pe*(OCC-1.0));   // THE REASON THIS EXISTS IS THAT DOSE ABSORBTION REALLY DOES INCREASE FOR PATIENTS 
                                                // FROM DOSE TO DOSE, **ONLY** IN THE PPQ DATA; THIS MAY NOT OCCUR FOR OTHER DRUGS
    // double THETA8 = 1.0;  // this is just fixed at one
    // double TVF1 = THETA8*F1COVD;
    double TVF1 = F1COVD;
    double F1 =  vprms[i_ppq_F1_indiv] * TVF1 * exp(IOV_rv); // IOV is the between dose variability
    vprms[i_ppq_F1_thisdose] = F1;
    

    
    // ### second, you redraw a specific KTR parameter for this dose, using a draw of the variable IOV2
    //     draw a random variate to get the value of the IOV2 variable
    double IOV2_rv = 0.0;
    if(pkpd_ppq::stochastic)
    {
        IOV2_rv = gsl_ran_gaussian( rng, sqrt(0.195) );
        //IOV2_rv = gsl_ran_gaussian( rng, sqrt(0.00001) );
    
        // MT *= exp( IOV2_rv ); - this is what you want to do, but MT goes into the rate variables below as 1/MT
        // so it's simpler to execute the three lines below
    
        vprms[i_ppq_k15] *= exp( - IOV2_rv );
        vprms[i_ppq_k56] *= exp( - IOV2_rv );
        vprms[i_ppq_k62] *= exp( - IOV2_rv );
    }
    
}

bool pkpd_ppq::we_are_past_a_dosing_time( double current_time )
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

void pkpd_ppq::generate_recommended_dosing_schedule()
{
    // TODO NEEDS TO BE DONE BY AGE AND WEIGHT
	
    // NEED TO MULTIPLY THE DOSING AMOUNT BY 0.57 (THE SCALING PARAM) THAT
    // GIVES YOU THE SCALING IN MOLECULAR WEIGHT FROM PPQ-PHOSPHATE TO PPQ
    
    double num_tablets_per_dose;
    
    if( weight < 5.0 )
    {
        num_tablets_per_dose = 0.0;
    }
    else if( weight < 8.0 )
    {
        num_tablets_per_dose = 0.5;
    }
    else if( weight < 11.0 )
    {
        num_tablets_per_dose = 0.75;
    }
    else if( weight < 17.0 )
    {
        num_tablets_per_dose = 1.0;
    }
    else if( weight < 25.0 )
    {
        num_tablets_per_dose = 1.5;
    }
    else if( weight < 36.0 )
    {
        num_tablets_per_dose = 2.0;
    }
    else if( weight < 60.0 )
    {
        num_tablets_per_dose = 3.0;
    }
    else if( weight < 80.0 )
    {
        num_tablets_per_dose = 4.0;
    }
    else
    {
        num_tablets_per_dose = 5.0;
    }
    
    double total_mg_dose = num_tablets_per_dose * 320.0;
    
    v_dosing_times.insert( v_dosing_times.begin(), 3, 0.0 );
    v_dosing_times[1] = 24.0;
    v_dosing_times[2] = 48.0;
    
    v_dosing_amounts.insert( v_dosing_amounts.begin(), 3, total_mg_dose*0.577 );
    
}




