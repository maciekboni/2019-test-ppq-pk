// Model and parameter values from Guidi et al., 2018

#include "assert.h"
#include "pkpd_artesunate.h"

#include <filesystem>
#include <string>
#include <fstream>
#include <iostream>

bool pkpd_artesunate::stochastic = true;

// constructor
pkpd_artesunate::pkpd_artesunate( )
{
    vprms.insert( vprms.begin(), artesunate_num_params, 0.0 );
    assert( vprms.size()==artesunate_num_params );
    
    // this is the dimensionality of the ODE system
    // there are 9 PK compartments and 1 variable for the parasite density
    dim = 10;

    // this is the main vector of state variables
    y0 = new double[dim];
    for(int i=0; i<dim; i++) y0[i]=0.0;

    //y0[dim-1] = 10000.0;
    set_parasitaemia(10000.0); // simply a wrapper for the statement above
    
    const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45;
    os 	= gsl_odeiv_step_alloc(T, dim);
    oc 	= gsl_odeiv_control_y_new (1e-6, 0.0);
    oe 	= gsl_odeiv_evolve_alloc(dim);
    
    // Patient Characteristics

    patient_id = 0;                     // Updated in main.cpp
    patient_weight = 54.0;              // default weight of the patient in kg, can be overwritten via command line input
    is_male=false;
    is_pregnant=false;
    patient_age = 25.0;
    patient_blood_volume = 5500000.0;   // 5.5L of blood for an adult individual of weight 54kg. 
                                        // Scaled later according to patient_weight in main function
    
    // Individual Patient Dosing Log
    num_doses_given = 0;
    num_hours_logged = 0;  
    last_logged_hour = -1.0; 
    total_mg_dose_per_occasion = -99.0; // Moved to constructor for uniformity with other classes
    doses_still_remain_to_be_taken = true;

    // Parasite PD Characteristics
    // the parameters 15, exp( 0.525 * log(2700)), and 0.9 give about a 90% drug efficacy for an initial parasitaemia of 10,000/ul (25yo patient, 54kg)
    
    pdparam_n = 20.0; // default parameter if CLO is not specified
    pdparam_EC50 = 0.1; // default parameter if CLO is not specified, ng/microliter
    pdparam_Pmax = 0.99997; // default parameter if CLO is not specified
    //pdparam_Pmax = 0.983; // here you want to enter the max daily killing rate; it will be converted to hourly later
                            
    rng=NULL;

}

// destructor
pkpd_artesunate::~pkpd_artesunate()
{
    delete[] y0;
    gsl_odeiv_evolve_free(oe);
    gsl_odeiv_control_free(oc);
    gsl_odeiv_step_free(os);
}

void pkpd_artesunate::set_parasitaemia( double parasites_per_ul )
{
    y0[dim-1] = parasites_per_ul; //  the final ODE equation is always the Pf asexual parasitaemia
}

// NOTE MUST REMEMBER THAT THE FUNCTION BELOW IS A STATIC MEMBER FUNCTIONS OF THIS CLASS
//
// "rhs_ode" must be defined with these four arguments for it to work in
// in the GSL ODE routines
//
int pkpd_artesunate::rhs_ode(double t, const double y[], double f[], void *pkd_object )
{
    pkpd_artesunate* p = (pkpd_artesunate*) pkd_object;

    // these are the right-hand sides of the derivatives of the six compartments
    
    // this is compartment 1, the gut or fixed dose compartment, i.e. the hypothetical compartment
    // where the drug goes in first
    f[0] =  - p->vprms[i_artesunate_KTR_indiv] * y[0];

    // these are the seven transit compartments
    f[1] = y[0]*p->vprms[i_artesunate_KTR_indiv] - y[1]*p->vprms[i_artesunate_KTR_indiv];
    f[2] = y[1]*p->vprms[i_artesunate_KTR_indiv] - y[2]*p->vprms[i_artesunate_KTR_indiv];
    f[3] = y[2]*p->vprms[i_artesunate_KTR_indiv] - y[3]*p->vprms[i_artesunate_KTR_indiv];
    f[4] = y[3]*p->vprms[i_artesunate_KTR_indiv] - y[4]*p->vprms[i_artesunate_KTR_indiv];
    f[5] = y[4]*p->vprms[i_artesunate_KTR_indiv] - y[5]*p->vprms[i_artesunate_KTR_indiv];
    f[6] = y[5]*p->vprms[i_artesunate_KTR_indiv] - y[6]*p->vprms[i_artesunate_KTR_indiv];
    f[7] = y[6]*p->vprms[i_artesunate_KTR_indiv] - y[7]*p->vprms[i_artesunate_KTR_indiv];
    
    // this is the central compartment (the blood)
    //
    // and the current units here (Aug 7 2024) are simply the total mg of artesunate in the blood
    // NOTE it looks like all of our blood concentrations in these PK classes simply track "total mg of molecule in blood"
    f[8] = y[7]*p->vprms[i_artesunate_KTR_indiv] - y[8]*p->vprms[i_artesunate_k20];
    
    // this is the per/ul parasite population size
    // indiv_central_volume_of_distribution (L) =! patient_blood_volume
    // drug concentration units mg/L, ec50 units ng/microliter, numerically the same
    double a = (-1.0/24.0) * log( 1.0 - (p->pdparam_Pmax * pow((y[8]/p -> vprms[i_artesunate_central_volume_of_distribution_indiv]),p->pdparam_n)) / (pow((y[8]/p -> vprms[i_artesunate_central_volume_of_distribution_indiv]),p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)));

    f[9] = -a * y[9];
    //f[9] = -(pow(a,0.5)) * y[9];

    //f[9] = (-a * p-> immune_killing_rate) * y[9];

    // double current_hour = floor(t);
    // std::filesystem::path folder_kill_art = "parasite_killing_constant_artesunate";
    // std::filesystem::create_directories(folder_kill_art);
    

    // if (current_hour > p->last_logged_hour) {
    //     std::filesystem::path filename_kill_art =  folder_kill_art / ("parasite_killing_constant_" + std::to_string(static_cast<int>(p->patient_weight)) + "kg_" + std::to_string(p->patient_id) + "_artesunate.txt");
    //     //std::string filename_kill_art = "parasite_killing_constant_" + std::to_string(static_cast<int>(p->patient_weight)) + "kg_artesunate.txt";
    //     std::ofstream outputFile_kill_art;
    //     outputFile_kill_art.open(filename_kill_art, std::ios::app);
    //     if (outputFile_kill_art.is_open()) {
    //         // Append data to the file            
    //         outputFile_kill_art << a << "," << t << std::endl;
    //         outputFile_kill_art.close(); 
    //         p->last_logged_hour = current_hour;  
    //     } 
    //     else {
    //     std::cerr << "Error opening" << filename_kill_art <<" for writing." << std::endl;
    //     }
    // } 

    return GSL_SUCCESS;

}

void pkpd_artesunate::give_next_dose_to_patient( double fractional_dose_taken )
{
    if( doses_still_remain_to_be_taken )
    {
        // The KA/KTR is modified by IOV by calling this function
    
        redraw_PK_params_before_newdose(); // these are the dose-specific parameters that you're drawing here
        
        // add the new dose amount to the "dose compartment", i.e. the first compartment
        // There is no IOV in F dha acc. to Tarning 2012, removed IOV between doses on F

        // Implementing F on dose, F has already been adjusted for IIV; refer to Fig 1B in Tarning 2012
        y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_artesunate_bioavailability_F_indiv] * fractional_dose_taken; 
        
        num_doses_given++;

        if( num_doses_given >= v_dosing_amounts.size() ) {
            doses_still_remain_to_be_taken=false;
        }

    }

}


void pkpd_artesunate::predict( double t0, double t1 )
{
    static double last_logged_hour = -1.0;
    gsl_odeiv_system sys = {pkpd_artesunate::rhs_ode, pkpd_artesunate::jac, dim, this};   // the fourth argument is a void pointer that you 
                                                                                          // are supossed to use freely; you normally
                                                                                          // use it to access the paramters
    double t = t0;    
    double h = 1e-6;

    while (t < t1)
    {
        // check if time t is equal to or larger than the next scheduled hour to log
        if( t >= ((double)num_hours_logged)  )
        {

            indiv_central_volume_millilitres = vprms[i_artesunate_central_volume_of_distribution_indiv] * 1000;   // Converting L to ml as Central Volume is in L
            v_concentration_in_blood.push_back( (y0[8] * pow(10, 6)) / indiv_central_volume_millilitres);         // The concentration in the CC is now converted to ng/ml
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



void pkpd_artesunate::initialize_PK_params( void )
{
    
    //WARNING - THE AGE MEMBER VARIABLE MUST BE SET BEFORE YOU CALL THIS FUNCTION

    //NOTE --- in this function alone, THE KTR PARAM HERE HAS INTER-PATIENT VARIABILITY BUT NO INTER-DOSE VARIABILITY
    
    // this is the median weight of the participants whose data were used to estimate the parameters for this study
    // it is used as a relative scaling factor below
    double population_median_weight=48.5;
    
    // initialize these relative dose factors to one (this should be the default behavior if
    // the model is not stochastic or if we decide to remove between-dose and/or between-patient variability

    // There is no IOV in F_thisdose acc to Tarning 2012
    // Leaving this line here in-case we adopt a different model/parameters
    //vprms[i_artesunate_bioavailability_F_thisdose] = 1.0;

    vprms[i_artesunate_bioavailability_F_indiv] = 1.0;
    
    //TVF1 = THETA(4)*(1+THETA(7)*(PARA-3.98)) * (1+THETA(6)*FLAG)
    double THETA7_pe = 0.278;
    double THETA6_pe = -0.375;
    //double THETA4_pe= 1.0;
    
    // To calculate the effect of initial parasitaemia on F_indiv, we need to use initial parasitemia per microliter
    // and not the total parasitemia as the value 3.98 in the equation below is log10(9549.93)
    // which is in parasites per microliter, not total parasitemia, which would be way, way, way higher

    //initial_log10_totalparasitaemia = log10(y0[dim-1]*patient_blood_volume);
    initial_log10_totalparasitaemia = log10(y0[dim-1]);
    double typical_bioavailibility_TVF = 1.0 + THETA7_pe*(initial_log10_totalparasitaemia-3.98);
    if(is_pregnant) typical_bioavailibility_TVF *= (1.0+THETA6_pe);
        
    double indiv_bioavailability_F = typical_bioavailibility_TVF;

    // Applying IIV in relative bioavailability F
    // The variance is from Table 4 of Tarning 2012,  and is calculated as ln((30.3/100)^2+1) = 0.0878 or ln((0.303)^2+1) = 0.0878

    if(pkpd_artesunate::stochastic)
    {
        double ETA4_rv = gsl_ran_gaussian( rng, sqrt(0.08800) );
        indiv_bioavailability_F *= exp(ETA4_rv); 
        // Original function was indiv_bioavailability_F *= ETA4_rv
        // But it has to be indiv_bioavailability_F *= exp(ETA4_rv) acc to Tarning 2012
    }
 
    vprms[i_artesunate_bioavailability_F_indiv] = indiv_bioavailability_F;
    
    // ### ### KTR is the transition rate to, between, and from the seven transit compartments
    double TVMT_pe = 0.982; // this is the point estimate (_pe) for TVMT; there is no need to draw a random variate here
    //double ETA3_rv = gsl_ran_gaussian( rng, sqrt(0.0) );  // _rv means random variate
                                                            // NOTE - the second argument to this function call needs to 
                                                            //        be the STANDARD DEVIATION not the variance
    // ETA3 is fixed at zero in this model
    // therefore, below, we do no add any variation into MT
    //double MT = TVMT_pe; // * exp(ETA3_rv); we think that "MT" here stands for mean time to transition from dose compartment to blood
    vprms[i_artesunate_MT_indiv] = TVMT_pe;


    // NOTE at this point you have an MT value without any effect of dose order (i.e. whether it's dose 1, dose 2, etc.
    // later, you must/may draw another mean-zero normal rv, and multiply by the value above
    // this is done in redraw_params_before_newdose()
    //vprms[i_artesunate_KTR_indiv] = 8.0/MT;
    vprms[i_artesunate_KTR_indiv] = 8.0/vprms[i_artesunate_MT_indiv];


    // ### ### this is the exit rate from the central compartment (the final exit rate in the model)
    double THETA1_pe = 78.0;
    double THETA2_pe = 129.0;
    double typical_clearance_TVCL = THETA1_pe * pow( patient_weight/population_median_weight, 0.75 );  
    
    
    //double ETA1_rv = 0.0; // this is fixed in this model
    //double CL = TVCL * exp(ETA1_rv);
    double indiv_clearance_CL = typical_clearance_TVCL; // just execute this line since ETA1 is fixed at zero above

    double typical_volume_TVV = THETA2_pe * (patient_weight/population_median_weight);  
    double indiv_volume_V = typical_volume_TVV;

    // Applying IIV in Vd
    // The variance is from Table 4 of Tarning 2012,  and is calculated as ln((12.8/100)^2+1) or ln((0.128)^2+1) = 0.01625
    if(pkpd_artesunate::stochastic) 
    
    {
        double ETA2_rv = gsl_ran_gaussian( rng, sqrt(0.0162) );
        indiv_volume_V *= exp(ETA2_rv);
    }

    double indiv_central_volume_of_distribution = indiv_volume_V;

    vprms[i_artesunate_typical_CL] = typical_clearance_TVCL;
    vprms[i_artesunate_CL_indiv] = indiv_clearance_CL;
    vprms[i_artesunate_typical_V] = typical_volume_TVV;
    vprms[i_artesunate_V_indiv] = indiv_volume_V;
    vprms[i_artesunate_central_volume_of_distribution_indiv] = indiv_central_volume_of_distribution;
    vprms[i_artesunate_k20] =  vprms[i_artesunate_CL_indiv]/vprms[i_artesunate_V_indiv];
}

void pkpd_artesunate::initialize_pkpd_artesunate_object() {
    generate_recommended_dosing_schedule();
    initialize_PK_params();
            
}


void pkpd_artesunate::redraw_PK_params_before_newdose()
{
     
    if(pkpd_artesunate::stochastic) 
    {
        // Applying IOV in MTT
        // The variance is from Table 4 of Tarning 2012,  and is calculated as ln((50.9/100)^2+1) or ln((0.509)^2+1) = 0.2303820898
        double ETA_rv = gsl_ran_gaussian( rng, sqrt(0.2303820898) ); 
       
        // Modifying KTR between doses independently and not cumulatively

        // This alters the absorption/transit rate by adding IOV in MTT
        // and not relative bioavailability between doses

        // you need to multiple MT by exp(ETA_rv)
        // BUT:  KTR = 8/MT, so instead, simply multiple KTR by exp(-ETA_rv)
        // Hence the negative sign
        
        vprms[i_artesunate_KTR_thisdose] = 8.0/vprms[i_artesunate_MT_indiv] * exp(-ETA_rv);
        vprms[i_artesunate_KTR_indiv] = vprms[i_artesunate_KTR_thisdose];

    }

}


bool pkpd_artesunate::we_are_past_a_dosing_time( double current_time )
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


void pkpd_artesunate::generate_recommended_dosing_schedule()
{

    /* From WHO: 
    Formulations currently available: 
    A fixed-dose combination in tablets containing 25 + 67.5 mg, 50 + 135 mg or 100 + 270 mg of AS and AQ, respectively
    The target dose (and range) are 4 (2–10) mg/kg bw per day artesunate once a day for 3 days. 
    A total therapeutic dose range of 6–30 mg/kg bw per day artesunate is recommended.
    */
    
    double num_tablets_per_dose;
    
    if( patient_weight >= 5.0 && patient_weight < 9.0 )
    {
        num_tablets_per_dose = 1.0;
    }
    else if( patient_weight >= 9.0 && patient_weight < 18.0  )
    {
        num_tablets_per_dose = 2.0;
    }
    else if( patient_weight >= 18.0 && patient_weight < 36.0  )
    {
        num_tablets_per_dose = 4.0;
    }
    else if (patient_weight >= 36.0)
    {
        num_tablets_per_dose = 8.0;
    } 
    // Adding an error message saying weight not supported
    else {
        std::cerr << "Error: Weight not supported." << std::endl;
    }
   
    // Artesunate given by weight for a total of three days - WHO guidelines, 2024
    // This is the total milligrams of artesunate taken each dose
  
    total_mg_dose_per_occasion = num_tablets_per_dose * 25.0;
    
    v_dosing_times.insert( v_dosing_times.begin(), 3, 0.0 );
    v_dosing_times[0] = 0.0;
    v_dosing_times[1] = 24.0; 
    v_dosing_times[2] = 48.0;

    v_dosing_amounts.insert( v_dosing_amounts.begin(), 3, total_mg_dose_per_occasion );

}