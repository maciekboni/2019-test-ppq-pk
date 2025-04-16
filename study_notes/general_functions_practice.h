#include <stdio.h>
#include <iostream>
#include <string>


void output_results_monotherapy_lum(int pi, pkpd_lum *dyn);

void output_results_monotherapy_art(int pi, pkpd_dha *dyn);

void output_results_combination_AL(int pi, pkpd_artemether *dyn1, pkpd_lum *dyn2);