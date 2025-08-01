#include "assert.h"
#include "pkpd_artemether.h"

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
    y0[dim-1] = 10000.0;
    
    
    const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45;
    os 	= gsl_odeiv_step_alloc(T, dim);
    oc 	= gsl_odeiv_control_y_new (1e-6, 0.0);
    oe 	= gsl_odeiv_evolve_alloc(dim);
    
    patient_weight = 54.0;      // default weight of the patient in kg, can be overwritten via command line input
    median_weight  =  54.0;     // in kilograms
    weight = patient_weight;    // this is the weight that is actually used in the calculations

    num_doses_given = 0;
    num_hours_logged = 0;    
    total_mg_dose_per_occassion = -99.0;    // Moved to constructor for uniformity with other classes
    
    age = 25.0;
    patient_blood_volume = 5500000.0; // 5.5L of blood for an adult individual of weight 54kg. 
                                      // Scaled later according to patient_weight in main function
    is_male=false;
    is_pregnant=false;
    doses_still_remain_to_be_taken = true;

    // the parameters 15, exp( 0.525 * log(2700)), and 0.9 give about a 90% drug efficacy for an initial parasitaemia of 10,000/ul (25yo patient, 54kg)
    pdparam_n = 20.0; // default parameter if CLO is not specified
    pdparam_EC50 = 0.1; // default parameter if CLO is not specified
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
    // double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow(y[8],p->pdparam_n) / (pow(y[8],p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)) );

    // Testing: adjusting the concentration in the central compartment/EC50 by the PATIENT BLOOD VOLUME
    //double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow((y[8]/p -> patient_blood_volume),p->pdparam_n) / (pow((y[8]/p -> patient_blood_volume),p->pdparam_n) + pow((p->pdparam_EC50/p -> patient_blood_volume),p->pdparam_n)));
    double a = (-1.0/24.0) * log( 1.0 - p->pdparam_Pmax * pow((y[8]/p -> patient_blood_volume),p->pdparam_n) / (pow((y[8]/p -> patient_blood_volume),p->pdparam_n) + pow(p->pdparam_EC50,p->pdparam_n)));



    // static double last_logged_hour = -1.0;  //Just a placeholder, will be updated when the first log is written
    // double current_hour = floor(t);

    // if (current_hour > last_logged_hour) {
    //     std::string filename_kill_art = "parasite_killing_constant_" + std::to_string(static_cast<int>(p->patient_weight)) + "kg_artemether.txt";
    //     std::ofstream outputFile_kill_art;
    //     outputFile_kill_art.open(filename_kill_art, std::ios::app);
    //     if (outputFile_kill_art.is_open()) {
    //     // Append data to the file            
    //     outputFile_kill_art << a << "," << t << std::endl;
    //     outputFile_kill_art.close(); 
    //     last_logged_hour = current_hour;  
    //     } 
    //     else {
    //     std::cerr << "Error opening" << filename_kill_art <<" for writing." << std::endl;
    //     }
    // }

    f[9] = -a * y[9];
    
    
    return GSL_SUCCESS;
}

void pkpd_artemether::give_next_dose_to_patient( double fractional_dose_taken )
{
    if( doses_still_remain_to_be_taken )
    {

        redraw_params_before_newdose(); // these are the dose-specific parameters that you're drawing here
        
        y0[0] +=  v_dosing_amounts[num_doses_given] * fractional_dose_taken;
        
        num_doses_given++;

        if( num_doses_given >= v_dosing_amounts.size() ) doses_still_remain_to_be_taken=false;

    }

}


void pkpd_artemether::predict( double t0, double t1 )
{
    gsl_odeiv_system sys = {pkpd_artemether::rhs_ode, pkpd_artemether::jac, dim, this};   // the fourth argument is a void pointer that you 
                                                                            // are supossed to use freely; you normally
                                                                            // use it to access the paramters
    double t = t0;    
    double h = 1e-6;

    while (t < t1)
    {
        // check if there are still doses to give
        if( num_doses_given < v_dosing_times.size() )
        {
            // check if time t is equal to or larger than the next scheduled dose
            if( t >= v_dosing_times[num_doses_given]  )
            {
                redraw_params_before_newdose();
                
                // add the new dose amount to the "dose compartment", i.e. the first compartment
                y0[0] +=  v_dosing_amounts[num_doses_given] * vprms[i_artemether_F1_thisdose];
                
                num_doses_given++;
            }
        }

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
            patient_blood_volume_litres = patient_blood_volume/pow(10,6); // in litres
            v_concentration_in_blood.push_back( y0[8]/ (patient_blood_volume_litres)); // The concentration in the blood, mg/L


            v_killing_rate.push_back( y0[9] );

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
    double mw=48.5;
    
    // initializse these relative dose factors to one (this should be the default behavior if
    // the model is not stochastic or if we decide to remove between-dose and/or between-patient variability
    vprms[i_artemether_F1_thisdose] = 1.0;
    vprms[i_artemether_F1_indiv] = 1.0;
    
    //TVF1 = THETA(4)*(1+THETA(7)*(PARA-3.98)) * (1+THETA(6)*FLAG)
    double THETA7_pe = 0.278;
    double THETA6_pe = -0.375;
    //double THETA4_pe= 1.0;
    
    initial_log10_totalparasitaemia = log10( y0[dim-1]*patient_blood_volume );
    double TVF1 = 1.0 + THETA7_pe*(initial_log10_totalparasitaemia-3.98);
    if(is_pregnant) TVF1 *= (1.0+THETA6_pe);
        
    double F1=TVF1;
    if(pkpd_artemether::stochastic)
    {
        double ETA4_rv = gsl_ran_gaussian( rng, sqrt(0.08800) );
        F1 *= ETA4_rv;
    }
 
    vprms[i_artemether_F1_indiv] = F1;
    
    //vprms[i_artemether_F1_thisdose] = vprms[i_artemether_F1_indiv]; // Just for testing - Venitha, April 2025
 
    
    
    // ### ### KTR is the transition rate to, between, and from the seven transit compartments
    double TVMT_pe = 0.982; // this is the point estimate (_pe) for TVMT; there is no need to draw a random variate here
    //double ETA3_rv = gsl_ran_gaussian( rng, sqrt(0.0) );  // _rv means random variate
                                                            // NOTE - the second argument to this function call needs to 
                                                            //        be the STANDARD DEVIATION not the variance
    // ETA3 is fixed at zero in this model
    // therefore, below, we do no add any variation into MT
    double MT = TVMT_pe; // * exp(ETA3_rv); we think that "MT" here stands for mean time to transition from dose compartment to blood

    // NOTE at this point you have an MT value without any effect of dose order (i.e. whether it's dose 1, dose 2, etc.
    //      later, you must/may draw another mean-zero normal rv, and multiply by the value above
    vprms[i_artemether_KTR] = 8.0/MT;


    // ### ### this is the exit rate from the central compartment (the final exit rate in the model)
    double THETA1_pe = 78.0;
    double THETA2_pe = 129.0;
    double TVCL = THETA1_pe * pow( weight/mw, 0.75 );  
    //double TVCL = THETA1_pe; // For debugging
    //double TVCL = THETA1_pe * (weight/mw); 
    //double TVCL = THETA1_pe * (15.0/mw);  
    
    
    //double ETA1_rv = 0.0; // this is fixed in this model
    //double CL = TVCL * exp(ETA1_rv);
    double CL = TVCL; // just execute this line since ETA1 is fixed at zero above

    double TVV2 = THETA2_pe * (weight/mw);  
    //double TVV2 = THETA2_pe * (15.0/mw); 
    double V2 = TVV2;

    if(pkpd_artemether::stochastic) 
    {
        double ETA2_rv = gsl_ran_gaussian( rng, sqrt(0.0162) );
        V2 *= exp(ETA2_rv);
    }
    
    vprms[i_artemether_k20] = CL/V2;
    
    //vprms[i_artemether_k20] = 0.5973735; // Median value of k20 for 50kg patient using this model
    //vprms[i_artemether_k20] = 0.9042795; // Median value of k20 for 10kg patient using this model

    // Added for debugging by Venitha
    vprms[i_artemether_typical_CL] = TVCL;
    vprms[i_artemether_CL_indiv] = CL;
    vprms[i_artemether_typical_V2] = TVV2;
    vprms[i_artemether_V2_indiv] = V2;
    
    
}

void pkpd_artemether::initialize() {
    generate_recommended_dosing_schedule();
    initialize_params();
            
}

// TODO the function below is a copy-and-paste from the PPQ function; must be modified
void pkpd_artemether::redraw_params_before_newdose()
{
     
    // ---- first, you don't receive the full dose.  You may receive 80% or 110% of the dose depending
    // on whether you're sitting or standing, whether you've recently had a big meal, whether some gets stuck
    // between your teeth; below, we set the parameter F1 (with some random draws) to adjust this initial dose

    // this is the "dose occassion", i.e. the order of the dose (first, second, etc)
    double OCC = 1.0 + (double)num_doses_given; // NOTE the RHS here is a class member
    
    if(pkpd_artemether::stochastic) 
    {
        double ETA_rv = gsl_ran_gaussian( rng, sqrt(0.23000) ); // this is ETA6, ETA7, and ETA8
        vprms[i_artemether_KTR] *= exp(-ETA_rv);       // WARNING this behavior is strange ... check if this is the right way to do it
                                                       // Seems okay - Venitha, April 2025
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

    //printf("The age and weight is %f and %f\n", age, weight);
    
    double num_tablets_per_dose;
    
    if( weight < 15.0 )
    {
        num_tablets_per_dose = 0.50;
    }
    else if( weight < 25.0 )
    {
        num_tablets_per_dose = 1.0;
    }
    else if( weight < 35.0 )
    {
        num_tablets_per_dose = 1.50;
    }
    else
    {
        num_tablets_per_dose = 2.0;
    }
   
    // Artemether given by weight, twice daily, for a total of three days - WHO guidelines, 2024
    total_mg_dose_per_occassion = num_tablets_per_dose * 40.0;
    
    v_dosing_times.insert( v_dosing_times.begin(), 6, 0.0 );
    v_dosing_times[0] = 0.0;
    v_dosing_times[1] = 12.0;
    v_dosing_times[2] = 24.0; 
    v_dosing_times[3] = 36.0;
    v_dosing_times[4] = 48.0;
    v_dosing_times[5] = 60.0;

    v_dosing_amounts.insert( v_dosing_amounts.begin(), 6, total_mg_dose_per_occassion );

}