#include "assert.h"
#include "pkpd_ppq.h"

#include <filesystem>
#include <string>
#include <fstream>
#include <iostream>

bool pkpd_ppq::stochastic = true;

// constructor
pkpd_ppq::pkpd_ppq(  )
{
    
    vprms.insert( vprms.begin(), ppq_num_params, 0.0 );
    assert( vprms.size()==ppq_num_params );
    
    // this is the dimensionality of the ODE system
    // 1 Dose Compartment, 2 Transit Compartments, 1 Central Compartment, 2 Peripheral Compartments, and 1 Parasite 'Compartment'
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
    set_parasitaemia(10000.0); // simply a wrapper for the statement above

    const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45;
    os 	= gsl_odeiv_step_alloc(T, dim);
    oc 	= gsl_odeiv_control_y_new (1e-6, 0.0);
    oe 	= gsl_odeiv_evolve_alloc(dim);

    // Patient Characteristics
    patient_id = 0;                         // Updated in main.cpp
    patient_age = 25.0;
    patient_weight = 54.0;                  // default weight of the patient in kg, can be overwritten via command line input
    patient_blood_volume = 5500000.0;       // 5.5L of blood for an adult individual of weight 54kg. 
                                            // Scaled later according to patient_weight in main function
    
    // Individual Patient Dosing Log
    num_doses_given = 0;
    num_hours_logged = 0;
    last_logged_hour = -1.0; 
    total_mg_dose_per_occasion = -99.0;     // Moved to constructor for uniformity with other classes
    doses_still_remain_to_be_taken = true; 

    // For testing
    central_volume_exponent = 1;
    
    // PD parameters
    pdparam_n = 15.0;                       // default parameter if CLO is not specified
    pdparam_EC50 = exp( 0.525 * log(2700)); // default parameter if CLO is not specified, ng/microliter
    pdparam_Pmax = 0.9;                     // here you want to enter the max daily killing rate; it will be converted to hourly later

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

    y0[dim-1] = parasites_per_ul;           //  the final ODE equation is always the Pf asexual parasitaemia

}


// NOTE :: MUST REMEMBER THAT THE FUNCTION BELOW IS A STATIC MEMBER FUNCTION OF THIS CLASS
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
    double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow((y[1]/p -> vprms[i_ppq_central_volume_of_distribution_indiv]),p->pdparam_n) / (pow((y[1]/p -> vprms[i_ppq_central_volume_of_distribution_indiv]),p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)) );
    
    f[6] = -a * y[6];
    
    return GSL_SUCCESS;
}

void pkpd_ppq::give_next_dose_to_patient( double fractional_dose_taken )
{
    if( doses_still_remain_to_be_taken )
    {
        // The individual F and KA/KTR is modified by IOV by calling this function
        // these are the dose-specific parameters that you're drawing here
        redraw_PK_params_before_newdose(); 
        
        // add the new dose amount to the "dose compartment", i.e. the first compartment
        // Implementing IOV_F on dose; F (i.e., vprms[i_ppq_bioavailability_F_thisdose]) has already been adjusted for IIV
        // Refer to Fig 2 of Hoglund et al., 2017
        y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_ppq_bioavailability_F_thisdose] * fractional_dose_taken;
        
        num_doses_given++;

        if( num_doses_given >= v_dosing_amounts.size() ) {
            doses_still_remain_to_be_taken=false;
        }

    }

}

void pkpd_ppq::predict( double t0, double t1 )
{
    //int i;
    gsl_odeiv_system sys = {pkpd_ppq::rhs_ode, pkpd_ppq::jac, dim, this};   // the fourth argument is a void pointer that you are supossed to use freely; you normally
                                                                            // use it to access the paramters
    double t = t0;    
    double h = 1e-6;
    
    
    while (t < t1)
    {
        // check if time t is equal to or larger than the next scheduled hour to log
        if( t >= ((double)num_hours_logged)  )
        {
            indiv_central_volume_millilitres = vprms[i_ppq_central_volume_of_distribution_indiv] * 1000;          // Converting L to ml as Central Volume is in L
            v_concentration_in_blood.push_back( (y0[1] * pow(10, 6)) / indiv_central_volume_millilitres);         // The concentration in the CC is now converted to ng/ml
                                                                                                                  // This is only the output! 
                                                                                                                  // The actual hill equation uses drug concentration in mg/L 
                                                                                                                  // mg/L == ng/microliter numerically 
                                                                                                                  // This was done as the unit of ec50 ng/microliter
                                                                                                                  // Reporting drug concentration in the blood as ng/ml
 
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

void pkpd_ppq::initialize_PK_params( void )
{
    // this is the median weight of a patient that these estimates were calibrated for
    double population_median_weight  =  54.0;
    
    //NOTE --- in this function alone, THE KTR PARAM HERE HAS INTER-PATIENT VARIABILITY BUT NO INTER-DOSE VARIABILITY
    //     --- THE F1 PARAMETER HERE IS JUST DEFUALTED TO ONE
    //     --- BOTH OF THESE WILL BE REDRAWN IN THE "redraw_params_before_newdose" FUNCTION IF IT IS CALLED 
    
    // indiv_bioavailability_F is the relative absorption level for this individual

    // before we take into account the effects of dose, the value of typical_bioavailability_TVF is one
    double typical_bioavailability_TVF = 1.0;
    double ETA8_rv = 0.0;

    // Used to incorporate IIV in F
    // The variance is from Table 3 of Hoglund 2017, and is calculated as ln((41.4/100)^2+1) = 0.1581962
    if( pkpd_ppq::stochastic )
    {
        ETA8_rv = gsl_ran_gaussian( rng, sqrt(0.1581962) ); // "_rv" means random variate
                                                        // this represents between-patient variability 
    }
    typical_bioavailability_TVF *= exp(ETA8_rv); 
    double indiv_bioavailability_F = typical_bioavailability_TVF;
    
    // KTR is the transition rate among the first three compartments
    // this is the point estimate (_pe) for TVMT; there is no need to draw a random variate here
    double typical_value_TVMT_pe = 2.11;
    double ETA7_rv = 0.0;

    // Used to incorporate IIV in MTT
    // The variance is from Table 3 of Hoglund 2017, and is calculated as ln((38.0/100)^2+1) = 0.1348804822
    if( pkpd_ppq::stochastic )
    {
        ETA7_rv = gsl_ran_gaussian( rng, sqrt(0.1348804822) );         // _rv means random variate
                                                                // NOTE - the second argument to this function call needs to 
    }                                                           //        be the STANDARD DEVIATION not the variance
    double indiv_MT = typical_value_TVMT_pe * exp( ETA7_rv );
    vprms[i_ppq_MT_indiv] = indiv_MT;
    
    // NOTE at this point you have an MT value without any effect of dose order (i.e. whether it's dose 1, dose 2, etc.
    //      later, you must draw another mean-zero normal rv, and multiply by the value above
    
    // this is the transition rate from the central compt to peripheral compt #1
    double THETA2_pe = 2910.0; // TVV of the central compartment
    double THETA3_pe = 310.0;
    
    double indiv_intercompartmental_clearance_Q1 = THETA3_pe * pow(patient_weight/population_median_weight, 0.75 ); 
    
    double typical_volume_TVV = THETA2_pe * pow(patient_weight/population_median_weight, 1.0 ); 
    double ETA2_rv = 0.0;

    // Used to incorporate IIV in MTT
    // The variance is from Table 3 of Hoglund 2017, and is calculated as ln((67.0/100)^2+1) = 0.3708046478
    if( pkpd_ppq::stochastic )
    {
        ETA2_rv = gsl_ran_gaussian( rng, sqrt(0.371) );  
    }

    typical_volume_TVV *= exp( ETA2_rv );
    double indiv_volume_V = typical_volume_TVV;
    
    // this is the transition rate from to peripheral compt #1 back to the central compt 
    double THETA4_pe = 4910.0; // TVVP1 of the first peripheral compartment
    double ETA4_rv = 0.0;

    // Used to incorporate IIV in MTT
    // The variance is from Table 3 of Hoglund 2017, and is calculated as ln((24.0/100)^2+1) = 0.05600219012
    if( pkpd_ppq::stochastic )
    {
        //ETA4_rv = gsl_ran_gaussian( rng, sqrt(0.0558) ); 
        ETA4_rv = gsl_ran_gaussian( rng, sqrt(0.0560) ); 
    }

    double typical_volume_peripheral_TVP1 = THETA4_pe * pow(patient_weight/population_median_weight, 1.0 ); 
    double individual_volume_peripheral_VP1 = typical_volume_peripheral_TVP1 * exp(ETA4_rv);
    
    // this is the transition rate from the central compt to peripheral compt #2
    double THETA5_pe = 105.0; // TVQ2 of the second peripheral compartment
    double ETA5_rv = 0.0;


    // Used to incorporate IIV in MTT
    // The variance is from Table 3 of Hoglund 2017, and is calculated as ln((23.6/100)^2+1) = 0.05420026505
    if( pkpd_ppq::stochastic )
    {
        //ETA5_rv = gsl_ran_gaussian( rng, sqrt(0.0541) ); 
        ETA5_rv = gsl_ran_gaussian( rng, sqrt(0.0542) ); 
    }

    double typical_intercompartmental_clearance_Q2 = THETA5_pe * pow(patient_weight/population_median_weight, 0.75 ); 
    double indiv_intercompartmental_clearance_Q2 = typical_intercompartmental_clearance_Q2 * exp(ETA5_rv); 
    
    // this is the transition rate from to peripheral compt #2 back to the central compt 
    double THETA6_pe = 30900.0;
    double ETA6_rv = 0.0;
    
    // Used to incorporate IIV in TVP2
    // The variance is from Table 3 of Hoglund 2017, and is calculated as ln((34.7/100)^2+1) = 0.1136937972
    if( pkpd_ppq::stochastic )
    {
        ETA6_rv = gsl_ran_gaussian( rng, sqrt(0.1136937972) ); 
    }

    double typical_volume_peripheral_TVP2 = THETA6_pe * pow(patient_weight/population_median_weight, 1.0 ); 
    double indiv_volume_peripheral_VP2 = typical_volume_peripheral_TVP2 * exp(ETA6_rv); 
    
    
    // this is the exit rate from the central compartment (the final exit rate in the model)
    double THETA1_pe = 55.4; 
    double ETA1_rv = 0.0;
    double HILL = 5.51;
    double MF50 = 0.575; 

    // Used to incorporate IIV in CL
    // The variance is from Table 3 of Hoglund 2017, and is calculated as ln((27.9/100)^2+1) = 0.07495996625
    if( pkpd_ppq::stochastic )
    {
        //ETA1_rv = gsl_ran_gaussian( rng, sqrt(0.0752) ); 
        ETA1_rv = gsl_ran_gaussian( rng, sqrt(0.07495996625) ); 
    }

    double maturation_function_MF = pow(patient_age,HILL) / ( pow(patient_age,HILL) + pow(MF50,HILL) );
    double typical_clearance_TVCL = THETA1_pe * maturation_function_MF * pow(patient_weight/population_median_weight, 0.75); 
    double indiv_clearance_CL = typical_clearance_TVCL * exp(ETA1_rv);

    vprms[i_ppq_bioavailability_F_indiv] = indiv_bioavailability_F;
    vprms[i_ppq_bioavailability_F_thisdose] = vprms[i_ppq_bioavailability_F_indiv]; // Only for the first dose
    vprms[i_ppq_CL_indiv] = indiv_clearance_CL;
    vprms[i_ppq_V_indiv] = indiv_volume_V;
    vprms[i_ppq_central_volume_of_distribution_indiv] = vprms[i_ppq_V_indiv]; // in Litres
    vprms[i_ppq_Q1_indiv] = indiv_intercompartmental_clearance_Q1;
    vprms[i_ppq_Q2_indiv] = indiv_intercompartmental_clearance_Q2;
    vprms[i_ppq_VP1_indiv] = individual_volume_peripheral_VP1;
    vprms[i_ppq_VP2_indiv] = indiv_volume_peripheral_VP2;

    vprms[i_ppq_KTR_indiv] = 3.0/vprms[i_ppq_MT_indiv];
    vprms[i_ppq_k15] = vprms[i_ppq_KTR_indiv];
    vprms[i_ppq_k23] = vprms[i_ppq_Q1_indiv]/vprms[i_ppq_V_indiv];
    vprms[i_ppq_k24] = vprms[i_ppq_Q2_indiv]/vprms[i_ppq_V_indiv];
    vprms[i_ppq_k32] = vprms[i_ppq_Q1_indiv]/vprms[i_ppq_VP1_indiv];
    vprms[i_ppq_k42] = vprms[i_ppq_Q2_indiv]/vprms[i_ppq_VP2_indiv];
    vprms[i_ppq_k56] = vprms[i_ppq_KTR_indiv];
    vprms[i_ppq_k62] = vprms[i_ppq_KTR_indiv];
    vprms[i_ppq_k20] = vprms[i_ppq_CL_indiv]/vprms[i_ppq_V_indiv];    
    
}

void pkpd_ppq::initialize_pkpd_ppq_object( void )
{
    //-- WARNING - THE AGE MEMBER VARIABLE MUST BE SET BEFORE YOU CALL THIS FUNCTION

    initialize_PK_params();
    generate_recommended_dosing_schedule();
}

void pkpd_ppq::redraw_PK_params_before_newdose()
{

    double DOSE_F = 0.237;                      //The DOSE_F, 23.7%, increase in relative F between each dosing occasion
    double OCC = 1.0 + (double)num_doses_given; // NOTE the RHS here is a class member

    // Increase in relative bioavailability F with each dose occasion 
    // Note: The first dose has no increase, i.e., F1COVD = 1.0 when OCC = 1
    double FCOVD = (1.0 + DOSE_F*(OCC-1.0));  
   
    // IOV in relative bioavailability F
    // Used to incorporate IOV in F
    // The variance is from Table 3 of Hoglund 2017, and is calculated as ln((53.5/100)^2+1) = 0.2517115716
    if(pkpd_ppq::stochastic)
    {
        double IOV_F_rv = gsl_ran_gaussian( rng, sqrt(0.2517115716) );  

        // Calculating the typical value of F for this dose
        // The original TVF stays at the value initialized before drug administration begins 
        double typical_bioavailability_F_this_dose = 1.0 * FCOVD;                                                                 
        double F_thisdose = vprms[i_ppq_bioavailability_F_indiv] * typical_bioavailability_F_this_dose * exp(IOV_F_rv);
        vprms[i_ppq_bioavailability_F_thisdose] = F_thisdose;
    } else {
        double typical_bioavailability_F_this_dose = 1.0 * FCOVD;                                                                 
        double F_thisdose = vprms[i_ppq_bioavailability_F_indiv] * typical_bioavailability_F_this_dose;
        vprms[i_ppq_bioavailability_F_thisdose] = F_thisdose;
    }

    // Used to incorporate IOV in MTT
    // The variance is from Table 3 of Hoglund 2017, and is calculated as ln((46.4/100)^2+1) = 0.1949876685
    if(pkpd_ppq::stochastic)
    {
        double IOV_MTT_rv = gsl_ran_gaussian( rng, sqrt(0.1949876685) );
        // MT *= exp( IOV_MTT_rv ); - this is what you want to do, but MT goes into the rate variables below as 1/MT
        // so it's simpler to execute the three lines below

        // Note: i_ppq_MT_indiv already has the individual MT value with IIV included
        vprms[i_ppq_KTR_thisdose] = 3.0/vprms[i_ppq_MT_indiv] * exp(-IOV_MTT_rv);
        vprms[i_ppq_KTR_indiv] = vprms[i_ppq_KTR_thisdose];
        vprms[i_ppq_k15] = vprms[i_ppq_KTR_indiv];
        vprms[i_ppq_k56] = vprms[i_ppq_KTR_indiv];
        vprms[i_ppq_k62] = vprms[i_ppq_KTR_indiv];
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
    // TODO NEEDS TO BE DONE BY AGE
	
    // NEED TO MULTIPLY THE DOSING AMOUNT BY 0.57 (THE SCALING PARAM) THAT
    // GIVES YOU THE SCALING IN MOLECULAR WEIGHT FROM PPQ-PHOSPHATE TO PPQ

    // DHA-PPQ comes in two fixed-dose combinations: 20/160 and 40/320
    // Using the pediatric combination of 20/160 and adjusting the dose accordingly
 
    // Updated dosing schedule by weight with accordance to the latest WHO guidelines dated 13 August 2025 (copied below):

    // Target dose and range: 
    // Adults and children weighing ≥ 25 kg: 4 (2–10) mg/kg bw per day DHA and 18 (16–27) mg/kg bw per day PPQ given once a day for 3 days 
    // Children weighing < 25 kg:            4 (2.5–10) mg/kg bw per day DHA and 24 (20–32) mg/kg bw per day PPQ once a day for 3 days
    
    // Revised dose recommendation for DHA + PPQ in young children (2015)
    // Children weighing < 25 kg should receive at least 
    // 2.5 mg/kg bw DHA and 20 mg/kg bw PPQ to achieve the same exposure as children weighing ≥ 25 kg and adults.

    double num_tablets_per_dose;
    
   
    if( patient_weight < 8.0 )
    {
        num_tablets_per_dose = 1;
    }
    else if( patient_weight >= 8.0 && patient_weight < 11.0 )
    {
        num_tablets_per_dose = 1.50;
    }
    else if(patient_weight >= 11.0 && patient_weight < 17.0 )
    {
        num_tablets_per_dose = 2.0;
    }
    else if(patient_weight >= 17.0 && patient_weight < 25.0 )
    {
        num_tablets_per_dose = 3.0;
    }
    else if( patient_weight >= 25.0 && patient_weight < 36.0 )
    {
        num_tablets_per_dose = 4.0;
    }
    else if( patient_weight >= 36.0 && patient_weight < 60.0 )
    {
        num_tablets_per_dose = 6.0;
    }
    else if( patient_weight >= 60.0 && patient_weight < 80.0 )
    {
        num_tablets_per_dose = 8.0;
    }
    else if( patient_weight >= 80.0 )
    {
        num_tablets_per_dose = 10.0;
    } // Adding an error message just in case
    else {
        std::cerr << "Error: Weight not supported." << std::endl;
    }

    // PQP - Piperaquine phosphate, what the patient takes orally
    total_mg_dose_per_occasion_pqp = (num_tablets_per_dose * 160.0);  // A single tablet contains 160 mg of PPQ

    // piperaquine base on which the population PK model is built (Hoglund 2017)
    total_mg_dose_per_occasion = total_mg_dose_per_occasion_pqp * 0.577; // Multiplying by 0.577 to convert from PPQ-phosphate to PPQ base
                                                                         // From the paper: Piperaquine was administered as piperaquine phosphate,
                                                                         // which was converted to piperaquine base with a scale factor of 57.7%   

    v_dosing_times.insert( v_dosing_times.begin(), 3, 0.0 );
    v_dosing_times[0] = 0.0;
    v_dosing_times[1] = 24.0;
    v_dosing_times[2] = 48.0;
    
    v_dosing_amounts.insert( v_dosing_amounts.begin(), 3, total_mg_dose_per_occasion); 
}


