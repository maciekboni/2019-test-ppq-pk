#include "assert.h"
#include "pkpd_artemether.h"
#include <filesystem>
#include <string>
#include <fstream>
#include <iostream>

bool pkpd_artemether::stochastic = true;

// constructor
pkpd_artemether::pkpd_artemether( )
{
    vprms.insert( vprms.begin(), artemether_num_params, 0.0 );
    assert( vprms.size()==artemether_num_params );
    
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

    patient_id = 0;             // Updated in main.cpp
    patient_weight = 54.0;      // default weight of the patient in kg, can be overwritten via command line input
    is_male=false;
    is_pregnant=false;
    patient_age = 25.0;
    patient_blood_volume = 5500000.0; // 5.5L of blood for an adult individual of weight 54kg. 
                                      // Scaled later according to patient_weight in main function
    
    // Individual Patient Dosing Log
    num_doses_given = 0;
    num_hours_logged = 0;  
    last_logged_hour = -1.0; 
    total_mg_dose_per_occasion = -99.0;    // Moved to constructor for uniformity with other classes
    doses_still_remain_to_be_taken = true;
   
    // For testing
    central_volume_exponent = 1;

    // Parasite PD Characteristics
    // the parameters 15, exp( 0.525 * log(2700)), and 0.9 give about a 90% drug efficacy for an initial parasitaemia of 10,000/ul (25yo patient, 54kg)
    
    pdparam_n = 20.0; // default parameter if CLO is not specified
    pdparam_EC50 = 0.1; // default parameter if CLO is not specified, ng/microliter
    pdparam_Pmax = 0.99997; // default parameter if CLO is not specified
    //pdparam_Pmax = 0.983; // here you want to enter the max daily killing rate; it will be converted to hourly later
                            
    rng=NULL;

}

// destructor
pkpd_artemether::~pkpd_artemether()
{
    delete[] y0;
    gsl_odeiv_evolve_free(oe);
    gsl_odeiv_control_free(oc);
    gsl_odeiv_step_free(os);
}

void pkpd_artemether::set_parasitaemia( double parasites_per_ul )
{
    y0[dim-1] = parasites_per_ul; //  the final ODE equation is always the Pf asexual parasitaemia
}

// NOTE MUST REMEMBER THAT THE FUNCTION BELOW IS A STATIC MEMBER FUNCTIONS OF THIS CLASS
//
// "rhs_ode" must be defined with these four arguments for it to work in
// in the GSL ODE routines
//
int pkpd_artemether::rhs_ode(double t, const double y[], double f[], void *pkd_object )
{
    pkpd_artemether* p = (pkpd_artemether*) pkd_object;

    // these are the right-hand sides of the derivatives of the six compartments
    
    // this is compartment 1, the gut or fixed dose compartment, i.e. the hypothetical compartment
    // where the drug goes in first
    f[0] =  - p->vprms[i_artemether_KTR] * y[0];

    // these are the seven transit compartments
    f[1] = y[0]*p->vprms[i_artemether_KTR] - y[1]*p->vprms[i_artemether_KTR];
    f[2] = y[1]*p->vprms[i_artemether_KTR] - y[2]*p->vprms[i_artemether_KTR];
    f[3] = y[2]*p->vprms[i_artemether_KTR] - y[3]*p->vprms[i_artemether_KTR];
    f[4] = y[3]*p->vprms[i_artemether_KTR] - y[4]*p->vprms[i_artemether_KTR];
    f[5] = y[4]*p->vprms[i_artemether_KTR] - y[5]*p->vprms[i_artemether_KTR];
    f[6] = y[5]*p->vprms[i_artemether_KTR] - y[6]*p->vprms[i_artemether_KTR];
    f[7] = y[6]*p->vprms[i_artemether_KTR] - y[7]*p->vprms[i_artemether_KTR];
    
    // this is the central compartment (the blood)
    //
    // and the current units here (Aug 7 2024) are simply the total mg of artemether in the blood
    // NOTE it looks like all of our blood concentrations in these PK classes simply track "total mg of molecule in blood"
    f[8] = y[7]*p->vprms[i_artemether_KTR] - y[8]*p->vprms[i_artemether_k20];
    
    // this is the per/ul parasite population size
    // indiv_central_volume_of_distribution (L) =! patient_blood_volume
    // drug concentration units mg/L, ec50 units ng/microliter, numerically the same
    double a = (-1.0/24.0) * log( 1.0 - (p->pdparam_Pmax * pow((y[8]/p -> vprms[i_artemether_central_volume_of_distribution_indiv]),p->pdparam_n)) / (pow((y[8]/p -> vprms[i_artemether_central_volume_of_distribution_indiv]),p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)));

    // Converting drug amount in blood from mg to ng and dividing it by patient blood volume in microliters
    //double a = (-1.0/24.0) * log( 1.0 - (p->pdparam_Pmax * pow(((y[8] * pow(10, 6))/p -> patient_blood_volume),p->pdparam_n)) / (pow(((y[8] * pow(10, 6))/p -> patient_blood_volume),p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)));

    f[9] = -a * y[9];
    //f[9] = -(pow(a,0.5)) * y[9];

    //f[9] = (-a * p-> immune_killing_rate) * y[9];

    // double current_hour = floor(t);
    // std::filesystem::path folder_kill_art = "parasite_killing_constant_artemether";
    // std::filesystem::create_directories(folder_kill_art);
    

    // if (current_hour > p->last_logged_hour) {
    //     std::filesystem::path filename_kill_art =  folder_kill_art / ("parasite_killing_constant_" + std::to_string(static_cast<int>(p->patient_weight)) + "kg_" + std::to_string(p->patient_id) + "_artemether.txt");
    //     //std::string filename_kill_art = "parasite_killing_constant_" + std::to_string(static_cast<int>(p->patient_weight)) + "kg_artemether.txt";
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

// The function is called using the bioavailability_F_indiv in main.cpp
void pkpd_artemether::give_next_dose_to_patient( double fractional_dose_taken )
{
    if( doses_still_remain_to_be_taken )
    {
        // The KA/KTR is modified by the IOV by calling this function
        // Rewrote function to modify bioavailability_F_thisdose instead

        redraw_params_before_newdose(); // these are the dose-specific parameters that you're drawing here
        
        // add the new dose amount to the "dose compartment", i.e. the first compartment
        
        // After redrawing the dose parameters, the individual bioavailability should be modified by the IOV/bioavailability_F_thisdose
        fractional_dose_taken *= vprms[i_artemether_bioavailability_F_thisdose]; 
        y0[0] +=  v_dosing_amounts[num_doses_given] * fractional_dose_taken; // 
       
        //y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_artemether_bioavailability_F_thisdose];
        // I guess fractional dose taken here is actually the bioavailability of the current dose
        // Need to check how this varies
        // The IOV applied when we redraw parameters is on the rate of absorption/transit between doses
        // And not the modified bioavailability
        // Also, although individual bioavailability is calculated, we don't implement it anywhere
        // Need to see how to apply it, either only to the first dose or if the change in F between doses can be directly
        // applied to F_indiv
        
        num_doses_given++;

        if( num_doses_given >= v_dosing_amounts.size() ) {
            doses_still_remain_to_be_taken=false;
        }

    }

}


void pkpd_artemether::predict( double t0, double t1 )
{
    static double last_logged_hour = -1.0;
    gsl_odeiv_system sys = {pkpd_artemether::rhs_ode, pkpd_artemether::jac, dim, this};   // the fourth argument is a void pointer that you 
                                                                                          // are supossed to use freely; you normally
                                                                                          // use it to access the paramters
    double t = t0;    
    double h = 1e-6;

    while (t < t1)
    {
        // check if time t is equal to or larger than the next scheduled hour to log
        if( t >= ((double)num_hours_logged)  )
        {
            // v_dosing_compartment.push_back( y0[0] );
            // v_transit_compartment1.push_back( y0[1] );
            // v_transit_compartment2.push_back( y0[2] );
            // v_transit_compartment3.push_back( y0[3] );
            // v_transit_compartment4.push_back( y0[4] );
            // v_transit_compartment5.push_back( y0[5] );
            // v_transit_compartment6.push_back( y0[6] );
            // v_transit_compartment7.push_back( y0[7] );

            //v_concentration_in_blood.push_back( y0[8]); // Not the concentration in the blood, but the total mg of artemether in the blood

            indiv_central_volume_millilitres = vprms[i_artemether_central_volume_of_distribution_indiv] * 1000;   // Converting L to ml
            v_concentration_in_blood.push_back( (y0[8] * pow(10, 6)) / indiv_central_volume_millilitres);         // The concentration in the blood, ng/ml
                                                                                                                  // This is only the output! 
                                                                                                                  // The actual hill equation uses drug concentration in mg/L 
                                                                                                                  // mg/L == ng/microliter numerically 
                                                                                                                  // This was done as the unit of ec50 ng/microliter
            
            //v_concentration_in_blood.push_back( (y0[8] * pow(10, 6)) / (patient_blood_volume/1000));            // Reporting drug concentration in the blood as ng/ml


            //v_killing_rate.push_back( y0[9] );

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



void pkpd_artemether::initialize_params( void )
{
    
    //WARNING - THE AGE MEMBER VARIABLE MUST BE SET BEFORE YOU CALL THIS FUNCTION

    //NOTE --- in this function alone, THE KTR PARAM HERE HAS INTER-PATIENT VARIABILITY BUT NO INTER-DOSE VARIABILITY
    
    // this is the median weight of the participants whose data were used to estimate the paramters for this study
    // it is used as a relative scaling factor below
    double median_weight=48.5;
    
    // initialize these relative dose factors to one (this should be the default behavior if
    // the model is not stochastic or if we decide to remove between-dose and/or between-patient variability
    vprms[i_artemether_bioavailability_F_thisdose] = 1.0;
    vprms[i_artemether_bioavailability_F_indiv] = 1.0;
    
    //TVF1 = THETA(4)*(1+THETA(7)*(PARA-3.98)) * (1+THETA(6)*FLAG)
    double THETA7_pe = 0.278;
    double THETA6_pe = -0.375;
    //double THETA4_pe= 1.0;
    
    initial_log10_totalparasitaemia = log10(y0[dim-1]*patient_blood_volume);
    double typical_bioavailibility_TVF = 1.0 + THETA7_pe*(initial_log10_totalparasitaemia-3.98);
    if(is_pregnant) typical_bioavailibility_TVF *= (1.0+THETA6_pe);
        
    double indiv_bioavailability_F = typical_bioavailibility_TVF;
    if(pkpd_artemether::stochastic)
    {
        double ETA4_rv = gsl_ran_gaussian( rng, sqrt(0.08800) );
        indiv_bioavailability_F *= ETA4_rv;
    }
 
    vprms[i_artemether_bioavailability_F_indiv] = indiv_bioavailability_F;
    
    // ### ### KTR is the transition rate to, between, and from the seven transit compartments
    double TVMT_pe = 0.982; // this is the point estimate (_pe) for TVMT; there is no need to draw a random variate here
    //double ETA3_rv = gsl_ran_gaussian( rng, sqrt(0.0) );  // _rv means random variate
                                                            // NOTE - the second argument to this function call needs to 
                                                            //        be the STANDARD DEVIATION not the variance
    // ETA3 is fixed at zero in this model
    // therefore, below, we do no add any variation into MT
    double MT = TVMT_pe; // * exp(ETA3_rv); we think that "MT" here stands for mean time to transition from dose compartment to blood

    // NOTE at this point you have an MT value without any effect of dose order (i.e. whether it's dose 1, dose 2, etc.
    // later, you must/may draw another mean-zero normal rv, and multiply by the value above
    vprms[i_artemether_KTR] = 8.0/MT;


    // ### ### this is the exit rate from the central compartment (the final exit rate in the model)
    double THETA1_pe = 78.0;
    double THETA2_pe = 129.0;
    double typical_clearance_TVCL = THETA1_pe * pow( patient_weight/median_weight, 0.75 );  
    
    
    //double ETA1_rv = 0.0; // this is fixed in this model
    //double CL = TVCL * exp(ETA1_rv);
    double indiv_clearance_CL = typical_clearance_TVCL; // just execute this line since ETA1 is fixed at zero above

    double typical_volume_TVV = THETA2_pe * (patient_weight/median_weight);  
    double indiv_volume_V = typical_volume_TVV;
    double indiv_central_volume_of_distribution = indiv_volume_V;

    if(pkpd_artemether::stochastic) 
    {
        double ETA2_rv = gsl_ran_gaussian( rng, sqrt(0.0162) );
        indiv_volume_V *= exp(ETA2_rv);
    }
    
    vprms[i_artemether_typical_CL] = typical_clearance_TVCL;
    vprms[i_artemether_CL_indiv] = indiv_clearance_CL;
    vprms[i_artemether_typical_V] = typical_volume_TVV;
    vprms[i_artemether_V_indiv] = pow(indiv_volume_V, central_volume_exponent);
    vprms[i_artemether_central_volume_of_distribution_indiv] = pow(indiv_central_volume_of_distribution, central_volume_exponent);
    
    vprms[i_artemether_k20] = indiv_clearance_CL/vprms[i_artemether_V_indiv];
    //vprms[i_artemether_k20] = indiv_clearance_CL/indiv_volume_V;
    
    //vprms[i_artemether_k20] = 0.5973735; // Median value of k20 for 50kg patient using this model
    //vprms[i_artemether_k20] = 0.9042795; // Median value of k20 for 10kg patient using this model
    
    
}

void pkpd_artemether::initialize() {
    generate_recommended_dosing_schedule();
    initialize_params();
            
}

// TODO the function below is a copy-and-paste from the PPQ function; must be modified
// Completed - Venitha, 08/2025
void pkpd_artemether::redraw_params_before_newdose()
{
     
    // ---- first, you don't receive the full dose.  You may receive 80% or 110% of the dose depending
    // on whether you're sitting or standing, whether you've recently had a big meal, whether some gets stuck
    // between your teeth; below, we set the parameter F1 (with some random draws) to adjust this initial dose

    // this is the "dose occasion", i.e. the order of the dose (first, second, etc)
    // This parameter is not used anywhere, perhaps the implementation needs to be re-written?
    double OCC = 1.0 + (double)num_doses_given; // NOTE the RHS here is a class member
    
    if(pkpd_artemether::stochastic) 
    {
        double ETA_rv = gsl_ran_gaussian( rng, sqrt(0.23000) ); // this is ETA6, ETA7, and ETA8
        //vprms[i_artemether_KTR] *= exp(-ETA_rv);                // WARNING this behavior is strange ... check if this is the right way to do it
                                                                // Seems okay - Venitha, April 2025 
                                                                // Actually, this alters the absorption/transit rate and not F_thisdose
                                                                // Modified to alter F_thisdose between occasions
                                                                // This parameter will later modify the individual bioavailability between doses
        vprms[i_artemether_bioavailability_F_thisdose] *= exp(-ETA_rv);
    }
    //double IOV_rv = ETA_rv;
    
    // you need to multiple MT by exp(ETA_rv)
    // BUT:  KTR = 8/MT, so instead, simply multiple KTR by exp(-ETA_rv)
    
    
    //  TVF1 = THETA(4)*(1+THETA(6)*FLAG)*(1+THETA(7)*(PARA-3.98));;
    //  F1   = TVF1*EXP(ETA(4));

    //double THETA4_pe = 1.0; //TODO check if this is really fixed at 1.0 as a point estimate; checked.  It is fixed at 1.0.
    //double ETA4_rv = gsl_ran_gaussian( rng, sqrt(0.08800) );

}


bool pkpd_artemether::we_are_past_a_dosing_time( double current_time )
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


void pkpd_artemether::generate_recommended_dosing_schedule()
{

    // DOSING GUIDELINES SAY    0,  8, 24, 36, 48, 60
    // BUT WE CAN JUST DO       0, 12, 24, 36, 48, 60

    // Dosing updated to match WHO guidelines in the hopes that it will improve efficacy...
    
    double num_tablets_per_dose;
    
    if( patient_weight >= 5.0 && patient_weight < 15.0 )
    {
        num_tablets_per_dose = 1.0;
    }
    else if( patient_weight >= 15.0 && patient_weight < 25.0  )
    {
        num_tablets_per_dose = 2.0;
    }
    else if( patient_weight >= 25.0 && patient_weight < 35.0  )
    {
        num_tablets_per_dose = 3;
    }
    else if (patient_weight >= 35.0)
    {
        num_tablets_per_dose = 4.0;
    } 
    // Adding an error message saying weight not supported
    else {
        std::cerr << "Error: Weight not supported." << std::endl;
    }
   
    // Artemether given by weight, twice daily, for a total of three days - WHO guidelines, 2024
    // This is the total milligrams of artemether taken each dose
    // Fixed-combination come in two doses: 20 mg and 40 mg; switching to 20 mg 'coz why not
    // Also adjusted the doses accordingly 
    total_mg_dose_per_occasion = num_tablets_per_dose * 20.0;
    
    v_dosing_times.insert( v_dosing_times.begin(), 6, 0.0 );
    v_dosing_times[0] = 0.0;
    v_dosing_times[1] = 8.0;
    v_dosing_times[2] = 24.0; 
    v_dosing_times[3] = 36.0;
    v_dosing_times[4] = 48.0;
    v_dosing_times[5] = 60.0;

    v_dosing_amounts.insert( v_dosing_amounts.begin(), 6, total_mg_dose_per_occasion );

}