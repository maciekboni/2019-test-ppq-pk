#ifndef GENERAL_FUNCTIONS
#define GENERAL_FUNCTIONS

#include "globals.h"
  
class general_functions
{
public:

    // Constructor
    general_functions();
    // Destructor
    ~general_functions();

    // Member function declarations
    // Input/Output functions
    void output_results_monotherapy_dha(int pi, pkpd_dha *dyn);
    void output_results_monotherapy_lum(int pi, pkpd_lum *dyn);

    void output_results_combination_AL(int pi, pkpd_artemether *dyn1, pkpd_lum *dyn2);
    void output_results_combination_DHA_PPQ(int pi, pkpd_dha *dyn1, pkpd_ppq *dyn2);
    void output_results_combination_AS_AQ(int pi, pkpd_artesunate *dyn1, pkpd_adq *dyn2);
    void ParseArgs(int argc, char **argv);

};


#endif // GENERAL_FUNCTIONS




