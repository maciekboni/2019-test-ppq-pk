#include "pkpd_artemether_practice.h"
#include "pkpd_lum_practice.h"


int G_CLO_N = 1;
double G_CLO_AGE = 25.0;
double G_CLO_WEIGHT = 54.0;

main() {
    if( G_CLO_THERAPY == therapy_AL )
    {

        fprintf(stdout, "PID,HOUR,COMP2CONC_ART,COMP2CONC_LUM,PARASITEDENSITY\n" );
        //Its actually not every hour, but the first/second 30min interval doesn't have a major difference, so we just label the 30 min half as an hour

        fprintf(stderr, "\n");
        // pi is patient index
        for(int pi=0; pi < G_CLO_N; pi++)
        {
            auto dyn1 = new pkpd_artemether(G_CLO_AGE, G_CLO_WEIGHT);

            dyn1 -> set_parasitaemia(20000.0); 
           
            //auto dyn2 = new pkpd_lum();
        }
    }
}
