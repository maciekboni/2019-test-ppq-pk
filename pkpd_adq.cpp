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
    // 1 dose compartment + 2 transit compartments + 1 absorption compartment + 2 AQ compartments + 3 DEAQ compartments + 1 PD equation
    // nine pk equations and one pd clearance equation
    dim = 10; 

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
    patient_weight = 50;
    patient_age = 25.0;
    pregnant = false;

    num_doses_given = 0;
    num_hours_logged = 0;    
    doses_still_remain_to_be_taken = true;
    total_mg_dose_per_occassion = -99.0;    // meaning it is not set yet
    patient_blood_volume = patient_weight * 70.0 * 1000.0; // in microliters

    // Patient PD parameters
    pdparam_n = 15.0;
    pdparam_EC50 = 0.1; // ng/microliter
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
    
    // this is the fixed dose compartment
    // IOV in Bioavailability F is implemented in give_next_dose_to_patient
    f[0] =  - p->vprms[i_adq_KTR_thisdose] * y[0]; 

    // the absorption compartment
    f[1] = y[8]*p->vprms[i_adq_KTR_thisdose] - y[1]*p->vprms[i_adq_Ka_thisdose];

    // this is compartment 2, the central compartment for AQ
    f[2] =  y[1]*p->vprms[i_adq_Ka_thisdose]  + y[3]*p->vprms[i_adq_k32] -  y[2]*p->vprms[i_adq_k23]  -  y[2]*p->vprms[i_adq_k24];
    
    // this is compartment 3, the peripheral compartment for AQ in the blood
    f[3] = y[2]*p->vprms[i_adq_k23]  -  y[3]*p->vprms[i_adq_k32];

    // this is compartment 4, the central compartment for DEAQ (a metabolite of AQ) in the blood
    f[4] = p->vprms[i_adq_CF]*y[2]*p->vprms[i_adq_k24]  +  y[5]*p->vprms[i_adq_k54]  +  y[6]*p->vprms[i_adq_k64]  -  y[4]*p->vprms[i_adq_k45]  -  y[4]*p->vprms[i_adq_k46]  -  y[4]*p->vprms[i_adq_k40];

    // this is compartment 5, first peripheral compartment for DEAQ in the blood
    f[5] = y[4]*p->vprms[i_adq_k45]  -  y[5]*p->vprms[i_adq_k54];

    // this is compartment 6, second peripheral compartment for DEAQ in the blood
    f[6] = y[4]*p->vprms[i_adq_k46]  -  y[6]*p->vprms[i_adq_k64];

    // this is compartment 7, the first transit compartment from dose to absorption
    f[7] = y[0]*p->vprms[i_adq_KTR_thisdose]  -  y[7]*p->vprms[i_adq_KTR_thisdose];

    // this is compartment 8, the second transit compartment from dose to absorption
    f[8] = y[7]*p->vprms[i_adq_KTR_thisdose]  -  y[8]*p->vprms[i_adq_KTR_thisdose];

    // this is the per/ul parasite population size
    // Using DEAQ concentration for killing effect, need to revisit later if AQ killing effect needs to be modelled too
    // indiv_central_volume_of_distribution (L) =! patient_blood_volume
    // drug concentration units mg/L, ec50 units ng/microliter, numerically the same
    double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow(y[4]/p->vprms[i_adq_central_volume_of_distribution_DEAQ_indiv],p->pdparam_n) / (pow(y[4]/p->vprms[i_adq_central_volume_of_distribution_DEAQ_indiv],p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)) );
    f[9] = -a * y[9];
    
    
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
        redraw_params_before_newdose(); 

        if( num_doses_given == 0 )
        {
            y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_adq_F_indiv_first_dose] * fractional_dose_taken;
        } 
        else
        {
            y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_adq_F_indiv_later_dose] * fractional_dose_taken;
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
            indiv_central_volume_millilitres = vprms[i_adq_central_volume_of_distribution_DEAQ_indiv] * 1000; 
            v_concentration_in_blood.push_back( (y0[2]*pow(10,6) / indiv_central_volume_millilitres));                   // The concentration in the CC is now converted to ng/ml
                                                                                                                         // This is only the output! 
                                                                                                                         // The actual hill equation uses drug concentration in mg/L 
                                                                                                                         // mg/L == ng/microliter numerically 
                                                                                                                         // This was done as the unit of ec50 ng/microliter
                                                                                                                         // Reporting drug concentration in the blood as ng/ml
            v_concentration_in_blood_metabolite.push_back( (y0[4] * pow(10, 6) / indiv_central_volume_millilitres));     // The concentration in the CC of DEAQ is now converted to ng/ml
                                                                                                                                                                                                                            
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


void pkpd_adq::initialize_pkpd_adq_object( void )
{
    
    //-- WARNING - THE AGE MEMBER VARIABLE MUST BE SET BEFORE YOU CALL THIS FUNCTION
    
    generate_recommended_dosing_schedule();
    initialize_PK_params();
    
}


void pkpd_adq::initialize_PK_params( void )
{

    // this is the median weight of a patient that these estimates were calibrated for
    double population_median_weight = 50.0;
    
     // all 18 below are point estimates
    double THETA1 = 2960.0;                 // CL_AQ
    double THETA2 = 13500.0;                // Vc_AQ
    double THETA3 = 2310.0;                 // Q_AQ
    double THETA4 = 22700.0;                // Vp_AQ    
    double THETA5 = 0.589;                  // Ka_AQ
    double THETA6 = 1.0;                    // F_AQ
    double THETA7 = 258.0;                  // Vc_DEAQ
    double THETA8 = 32.6;                   // CL_DEAQ
    double THETA9  = 154.0;                 // Q1_DEAQ
    double THETA10 = 2460.0;                // Vp1_DEAQ
    double THETA11 = 5580.0;                // Vp2_DEAQ
    double THETA12 = 31.3;                  // Q2_DEAQ    
    double THETA13 = 0.236;                 // MTT

    double THETA14 = 0.224;                 // Effect of first dose on F
                                            // In Ali 2018, given as -22.4% i.e., 0.224 reduction
                                            // Or -0.224
                                            // Missing sign in value here, but implemented correctly as (1-THETA14) * F
                                                                                                    
    double THETA15 = 11.8;                  // PMA50 for AQ, PMA at which maturation is 50%
    double THETA16 = 3.6;                   // Hill factor for PMA50, AQ
    double THETA17 = 12.9;                  // PMA50 for DEAQ
    double THETA18 = 3.22;                  // Hill factor for PMA50, DEAQ


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
        //  "_rv" means random variate 

        // The variances are from Table 3 of Ali 2018
        // The variance is calculated using the exact log-normal formula: ln((CV/100)² + 1)
        ETA1_rv = gsl_ran_gaussian( rng, sqrt(0.098654) );      // IIV in CL_AQ 
        ETA2_rv = gsl_ran_gaussian( rng, sqrt(0.248391) );      // IIV in Vc_AQ

        ETA8_rv = gsl_ran_gaussian( rng, sqrt(0.039221) );      // IIV in CL_DEAQ
        ETA7_rv = gsl_ran_gaussian( rng, sqrt(0.372655) );      // IIV in Vc_DEAQ

        ETA5_rv = gsl_ran_gaussian( rng, sqrt(0.480093) );      // IOV in Ka, also implemented in give_next_dose_to_patient
        //ETA6_rv = gsl_ran_gaussian( rng, sqrt(0.091194) );    // IOV in F, implemented in give_next_dose_to_patient
        //ETA13_rv = gsl_ran_gaussian( rng, sqrt(0.627198) );   // IOV in MTT, implemented in give_next_dose_to_patient
    }


    // in the parameter calculations below
    // "TV" means typical value or population mean for some parameter
    // Postmenstrual age of patient, PMA 
    double patient_PMA = (patient_age * 12.0) + 9.0; // PMA in months, (age in years) * 12 + 9; 9 months represents average gestational age at birth

    // PMA50, PMA at which enzyme maturation is 50%
    double PMA50_AQ = THETA15;
    double PMA50_DEAQ = THETA17;

    double enzyme_maturation_AQ   = pow( patient_PMA, THETA16 ) / ( pow( patient_PMA, THETA16 ) + pow( PMA50_AQ, THETA16 ) );
    double enzyme_maturation_DEAQ = pow( patient_PMA, THETA18 ) / ( pow( patient_PMA, THETA18 ) + pow( PMA50_DEAQ, THETA18 ) );

    // CL_AQ is the clearance of AQ from the central compartment
    double TVCL_AQ = THETA1 * pow(patient_weight/population_median_weight , 0.75 ) * enzyme_maturation_AQ;
    double CL_AQ   = TVCL_AQ * exp(ETA1_rv);

    // Vc_AQ is the volume of distribution of the central compartment for AQ
    double TVVc_AQ = THETA2 * pow(patient_weight/population_median_weight , 1.0 );
    double Vc_AQ   = TVVc_AQ * exp( ETA2_rv );

    // Q_AQ is inter-compartmental clearance between central and peripheral compartments for AQ
    double TVQ_AQ = THETA3 * pow(patient_weight/population_median_weight , 0.75 );
    double Q_AQ   = TVQ_AQ * exp( ETA3_rv );

    // VP_AQ is the volume of distribution of the peripheral compartment for AQ
    double TVVp_AQ = THETA4 * pow(patient_weight/population_median_weight , 1.0 );
    double Vp_AQ   = TVVp_AQ * exp( ETA4_rv );

    // Ka is the absorption rate constant, the rate at which drug moves from the absorption compartment to the central AQ compartment
    double TVKa = THETA5;
    double Ka   = TVKa * exp( ETA5_rv );

    // double F_later_dose               = exp( ETA6_rv );
    // double F_indiv_first_dose         = (1.0 - THETA14) * F_later_dose;

    double F_indiv_first_dose = (1.0 - THETA14) * exp(ETA6_rv);

    // Vc_DEAQ is the volume of distribution of the central compartment for DEAQ
    double TVVc_DEAQ        = THETA7 * pow(patient_weight/population_median_weight , 1.0 );
    double Vc_DEAQ          = TVVc_DEAQ * exp( ETA7_rv );

    // CL_DEAQ is the clearance of DEAQ from the central compartment
    double TVCL_DEAQ        = THETA8 * pow(patient_weight/population_median_weight , 0.75 ) * enzyme_maturation_DEAQ;
    double CL_DEAQ          = TVCL_DEAQ * exp( ETA8_rv );

    // Q1_DEAQ is inter-compartmental clearance between the central and first peripheral compartment for DEAQ
    double TVQ1_DEAQ        = THETA9 * pow(patient_weight/population_median_weight , 0.75 );
    double Q1_DEAQ          = TVQ1_DEAQ * exp( ETA9_rv );

    // Vp1_DEAQ is the volume of distribution of the first peripheral compartment for DEAQ
    double TVVp1_DEAQ       = THETA10 * pow(patient_weight/population_median_weight , 1.0 );
    double Vp1_DEAQ         = TVVp1_DEAQ * exp( ETA10_rv );

    // Vp2_DEAQ is the volume of distribution of the second peripheral compartment for DEAQ
    double TVVp2_DEAQ       = THETA11 * pow(patient_weight/population_median_weight , 1.0 );
    double Vp2_DEAQ         = TVVp2_DEAQ * exp( ETA11_rv );

    // Q2_DEAQ is inter-compartmental clearance between the central and second peripheral compartment for DEAQ
    double TVQ2_DEAQ        = THETA12 * pow(patient_weight/population_median_weight , 0.75 );
    double Q2_DEAQ          = TVQ2_DEAQ * exp( ETA12_rv );

    // MTT is the mean transit time from dose to central compartment
    double TVMT             = THETA13;
    double MT               = TVMT * exp( ETA13_rv );

    // simple re-scaling from mg/L to ng/mL (or perhaps the other way around)
    double S2               = Vc_AQ / 1000.0;
    double S4               = Vc_DEAQ / 1000.0;

    // NN is the number of transit compartments from dose to absorption compartment
    double NN = 2.0;
    double KTR = (NN+1.0)/MT; // this is the transition param btw dosing and absorption compartments

    // ----BLOCK2 COMMENTED OUT FOR EFFICIENCY----
    // double MW_AQ   =  355.87                             ; Molecular weight AQ (g/mol)
    // double MW_DEAQ =  327.813                            ; Molecular weight DEAQ (g/mol)
    // double CF = MW_DEAQ / MW_AQ                          ; Conversion factor AQ -> DEAQ molecular mass
    // double CF = 0.921159412;
    // ----BLOCK2 COMMENTED OUT FOR EFFICIENCY----


    // these are the three transitions params from dose to central; normally they are all the same
    // so this particular structure is a bit strange; but we keep it for consistency with the original 
    // NONMEM model fit

    vprms[i_adq_KTR_thisdose] = KTR;
    vprms[i_adq_Ka_thisdose] = Ka;

    vprms[i_adq_k23] = Q_AQ/Vc_AQ;
    vprms[i_adq_k24] = CL_AQ/Vc_AQ;
    vprms[i_adq_k32] = Q_AQ/Vp_AQ;
    
    vprms[i_adq_k45] = Q1_DEAQ/Vc_DEAQ;
    vprms[i_adq_k54] = Q1_DEAQ/Vp1_DEAQ;
    vprms[i_adq_k46] = Q2_DEAQ/Vc_DEAQ;
    vprms[i_adq_k64] = Q2_DEAQ/Vp2_DEAQ;

    vprms[i_adq_central_volume_of_distribution_AQ_indiv] = Vc_AQ;
    vprms[i_adq_central_volume_of_distribution_DEAQ_indiv] = Vc_DEAQ;
    vprms[i_adq_k40] = CL_DEAQ/Vc_DEAQ;

    vprms[i_adq_CF] = 0.921159412; //Conversion factor AQ -> DEAQ molecular mass
                                   // this is the commented out value from BLOCK2

    vprms[i_adq_F_indiv_first_dose] = F_indiv_first_dose;
    
}


void pkpd_adq::redraw_params_before_newdose()
{
     
    // First dose of amodiaquine is absorbed more poorly with 22.4% lower absorption 
    // this is the THETA14 parameter and has already been implemented in initialize_PK_params()
    
    // Implementing IOV in Ka, F and MTT 

    double THETA5 = 0.589;  // Ka_AQ
    double THETA13 = 0.236; // MTT

    double ETA5_rv = 0.0;   // IOV in Ka
    double ETA6_rv = 0.0;   // IOV in F  
    double ETA13_rv = 0.0;  // IOV in MTT

    if (pkpd_adq::stochastic) {
        // The variance is from Table 3 of Ali 2018,  and is calculated as ln((CV/100)² + 1)
        ETA5_rv = gsl_ran_gaussian( rng, sqrt(0.480093) );    // IOV in Ka
        ETA6_rv = gsl_ran_gaussian( rng, sqrt(0.091194) );    // IOV in F
        ETA13_rv = gsl_ran_gaussian( rng, sqrt(0.627198) );   // IOV in MTT
    }

    // Redraw absorption rate constant Ka
    vprms[i_adq_Ka_thisdose] = THETA5 * exp( ETA5_rv );

    // Redraw bioavailability for subsequent doses
    double F_later_dose = exp(ETA6_rv);
    vprms[i_adq_F_indiv_later_dose] = 1.0 * F_later_dose;

    // Redraw MTT
    double MT = THETA13 * exp( ETA13_rv );
    double NN = 2.0;
    vprms[i_adq_KTR_thisdose] = (NN+1.0)/MT;
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
    
    // Edited schedule based on WHO guidelines for amodiaquine dosing
    /* From WHO:  
    The target dose (and range) are 4 (2–10) mg/kg bw per day artesunate and 10 (7.5–15) mg/kg bw per day amodiaquine once a day for 3 days. 
    A total therapeutic dose range of 6–30 mg/kg bw per day artesunate and 22.5–45 mg/kg bw per dose amodiaquine is recommended. */

    double num_tablets_per_dose;
    
    // this is in kilograms
    if( patient_weight >= 5.0 && patient_weight < 9.0 )
    {
        num_tablets_per_dose = 1.0; // 67.5 mg 
    }
    else if( patient_weight >= 9.0 && patient_weight < 18.0  )
    {
        num_tablets_per_dose = 2.0; // 135 mg
    }
    else if( patient_weight >= 18.0 && patient_weight < 36.0  )
    {
        num_tablets_per_dose = 4.0; // 270 mg
    }
    else if (patient_weight >= 36.0)
    {
        num_tablets_per_dose = 8.0; // 540 mg
    } 

    // Adding an error message saying weight not supported

    else {
        std::cerr << "Error: Weight not supported." << std::endl;
    }
    // 
    total_mg_dose_per_occassion = num_tablets_per_dose * 67.5;  // this is 67.5mg of amodiaquine, which is the pediatric strength
                                                                // PS - in reality, there are different tablet sizes
    
    v_dosing_times.insert( v_dosing_times.begin(), 3, 0.0 );
    v_dosing_times[0] = 0.0;
    v_dosing_times[1] = 24.0;
    v_dosing_times[2] = 48.0; 

    v_dosing_amounts.insert( v_dosing_amounts.begin(), 3, total_mg_dose_per_occassion );
    
}




