import numpy as np
from scipy.integrate import solve_ivp
import math
from enum import IntEnum

stochastic_artemether = True

class IndexParameterArtemether(IntEnum):
    index_artemether_KTR = 0
    index_artemether_k20 = 1
    index_artemether_F1_indiv = 2
    index_artemether_F1_thisdose = 3
    index_artemether_typical_CL = 4
    index_artemether_CL_indiv = 5
    index_artemether_typical_V2 = 6
    index_artemether_V2_indiv = 7
    artemether_num_params = 8

class IndexCompartmentArtemether(IntEnum):
     index_compartment_dose = 0
     index_compartment_TR1 = 1
     index_compartment_TR2 = 2
     index_compartment_TR3 = 3
     index_compartment_TR4 = 4
     index_compartment_TR5 = 5
     index_compartment_TR6 = 6
     index_compartment_TR7 = 7
     index_compartment_central_compartment = 8
     index_compartment_parasitaemia = 9

class pkpd_artemether:

    def __init__(self):
          
           # Creating a np.array to hold the parameters using an IntEnum object

           self.vprms = np.zeros(int(IndexParameterArtemether.artemether_num_params), dtype=float)


           # These are the patient instance attributes for each pkpd_artemether object
           
           self.patient_weight = 54.0                   # default weight of the patient in kg, can be overwritten via command line input
           self.median_weight  =  54.0                  # in kilograms
           self.weight = self.patient_weight            # this is the weight that is actually used in the calculations

           self.num_doses_given = 0
           self.num_hours_logged = 0  
           self.total_mg_dose_per_occassion = -99.0      # Moved to constructor for uniformity with other classes
    
           self.age = 25.0
           self.patient_blood_volume = 5500000.0         # 5.5L of blood for an adult individual of weight 54kg. 
                                                         #Scaled later according to patient_weight in main function
           self.is_male = False
           self.is_pregnant = False
           self.doses_still_remain_to_be_taken = True

           # the parameters 15, exp( 0.525 * log(2700)), and 0.9 give about a 90% drug efficacy for an initial parasitaemia of 10,000/ul (25yo patient, 54kg)
           
           self.pdparam_n = 20.0                         # default parameter if CLO is not specified
           self.pdparam_EC50 = 0.1                       # default parameter if CLO is not specified
           self.pdparam_Pmax = 0.99997                   # default parameter if CLO is not specified
           #self.pdparam_Pmax = 0.983;                   # here you want to enter the max daily killing rate; it will be converted to hourly later

           # Now moving onto the state variables

           self.model_dim = 10
           self.model_compartments = np.zeros(self.model_dim, dtype=float)
           self.model_compartments[IndexCompartmentArtemether.index_compartment_parasitaemia] = 10000.0 # initial parasitaemia in parasites per ul
           


    # Now defining the rhs ODE function

    def rhs_ode(self, time, model_compartments):
        
        # Creating a np.array to hold the derivatives of the state variables
        model_derivatives = np.zeros(self.model_dim, dtype=float) 

        # The first compartment is the dosing compartment
        model_derivatives[IndexCompartmentArtemether.index_compartment_dose] = (- self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                                model_compartments[IndexCompartmentArtemether.index_compartment_dose])

        # There are seven transit compartments

        # TR1
        model_derivatives[IndexCompartmentArtemether.index_compartment_TR1] = (self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_dose] 
                                                                               - 
                                                                               self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR1])
        
        # TR2
        model_derivatives[IndexCompartmentArtemether.index_compartment_TR2] = (self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR1] 
                                                                               - 
                                                                               self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR2])
        
        # TR3
        model_derivatives[IndexCompartmentArtemether.index_compartment_TR3] = (self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR2] 
                                                                               - 
                                                                               self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR3])
        
        # TR4
        model_derivatives[IndexCompartmentArtemether.index_compartment_TR4] = (self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR3] 
                                                                               - 
                                                                               self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR4])
        
        # TR5
        model_derivatives[IndexCompartmentArtemether.index_compartment_TR5] = (self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR4] 
                                                                               - 
                                                                               self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR5])
        
        # TR6
        model_derivatives[IndexCompartmentArtemether.index_compartment_TR6] = (self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR5] 
                                                                               - 
                                                                               self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR6])
        
        # TR7
        model_derivatives[IndexCompartmentArtemether.index_compartment_TR7] = (self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR6] 
                                                                               - 
                                                                               self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR7])

        
        # Finally the Central Compartment, that is the blood
        model_derivatives[IndexCompartmentArtemether.index_compartment_central_compartment] = (self.vprms[IndexParameterArtemether.index_artemether_KTR] * 
                                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_TR7] 
                                                                                               -
                                                                                               self.vprms[IndexParameterArtemether.index_artemether_k20] * 
                                                                                               model_compartments[IndexCompartmentArtemether.index_compartment_central_compartment])



        # Now implementing the parasite killing function

        # The killing function is a Hill function, which is a sigmoid function
        # The numerator is pmax * (concentration in blood to the power n)

        hill_numerator  = (self.pdparam_Pmax * (model_compartments[IndexCompartmentArtemether.index_compartment_central_compartment] / self.patient_blood_volume) ** 
                                                 self.pdparam_n)

        # The denominator is (concentration in blood to the power n) + (EC50 to the power n)
        # Splitting the denominator into two parts for clarity

        hill_denominator_part1 = ((model_compartments[IndexCompartmentArtemether.index_compartment_central_compartment] / self.patient_blood_volume) ** 
                                  self.pdparam_n)
        
        hill_denominator_part2 = ((self.pdparam_EC50 / self.patient_blood_volume) ** 
                                  self.pdparam_n)
        
        # Now putting it all together

        a = (-1.0/24.0) * np.log(1.0 - (hill_numerator / (hill_denominator_part1 + hill_denominator_part2)))
        
        # a = ((-1/24) * 
        #      np.log( 1 - (self.pdparam_Pmax *
        #                   ((model_compartments[IndexCompartmentArtemether.index_compartment_central_compartment] /
        #                     self.patient_blood_volume)
        #                     ** self.pdparam_n ) /
        #                     model_compartments[IndexCompartmentArtemether.index_compartment_central_compartment] ** 
        #                     self.pdparam_n +
        #                     (self.pdparam_EC50 / self.patient_blood_volume) ** 
        #                     self.pdparam_n)))
        
        # Now killing the parasites in the blood

        model_derivatives[IndexCompartmentArtemether.index_compartment_parasitaemia] = - a * model_compartments[IndexCompartmentArtemether.index_compartment_parasitaemia]

        return model_derivatives
    

    # Now defining the function to set parameters

    def initialize_parameters_artemether(self):
        
        mw=48.5;
    
        self.vprms[IndexParameterArtemether.index_artemether_F1_thisdose] = 1.0
        self.vprms[IndexParameterArtemether.index_artemether_F1_indiv] = 1.0
    
        THETA7_pe = 0.278
        THETA6_pe = -0.375
        THETA4_pe = 1.0
    
        initial_log10_totalparasitaemia = np.log10(self.model_compartments[IndexCompartmentArtemether.index_compartment_parasitaemia] * self.patient_blood_volume) 
        TVF1 = 1.0 + THETA7_pe*(initial_log10_totalparasitaemia-3.98)
    
        if self.is_pregnant == True:
            TVF1 *= (1.0+THETA6_pe)
        
        F1=TVF1
        
        if self.stochastic == True:
            ETA4_rv = np.random.normal(rng, math.sqrt(0.08800) )
            F1 *= ETA4_rv
    
        vprms[i_artemether_F1_indiv] = F1;
    
        vprms[i_artemether_F1_thisdose] = vprms[i_artemether_F1_indiv]; # Just for testing - Venitha, April 2025
        
        double TVMT_pe = 0.982; # this is the point estimate (_pe) for TVMT; there is no need to draw a random variate here
        //double ETA3_rv = gsl_ran_gaussian( rng, sqrt(0.0) );  # _rv means random variate
                                                                # NOTE - the second argument to this function call needs to be the STANDARD DEVIATION not the variance
         #ETA3 is fixed at zero in this model
         #therefore, below, we do no add any variation into MT
         double MT = TVMT_pe; # * exp(ETA3_rv); we think that "MT" here stands for mean time to transition from dose compartment to blood

        #NOTE at this point you have an MT value without any effect of dose order (i.e. whether it's dose 1, dose 2, etc.
        #      later, you must/may draw another mean-zero normal rv, and multiply by the value above
        # vprms[i_artemether_KTR] = 8.0/MT;


        ### ### this is the exit rate from the central compartment (the final exit rate in the model)
        double THETA1_pe = 78.0;
        double THETA2_pe = 129.0;
        double TVCL = THETA1_pe * pow( weight/mw, 0.75 );  
        
        
        //double ETA1_rv = 0.0; // this is fixed in this model
        //double CL = TVCL * exp(ETA1_rv);
        double CL = TVCL; // just execute this line since ETA1 is fixed at zero above

        double TVV2 = THETA2_pe * (weight/mw);  
        double V2 = TVV2;

        if(pkpd_artemether::stochastic) 
        {
            double ETA2_rv = gsl_ran_gaussian( rng, sqrt(0.0162) );
            V2 *= exp(ETA2_rv);
        }
        
        vprms[i_artemether_k20] = CL/V2;
        
        #vprms[i_artemether_k20] = 0.5973735; # Median value of k20 for 50kg patient using this model
        #vprms[i_artemether_k20] = 0.9042795; # Median value of k20 for 10kg patient using this model

        # Added for debugging by Venitha
        vprms[i_artemether_typical_CL] = TVCL;
        vprms[i_artemether_CL_indiv] = CL;
        vprms[i_artemether_typical_V2] = TVV2;
        vprms[i_artemether_V2_indiv] = V2;
        
        return None

        # Now defining the function to run the ODE solver