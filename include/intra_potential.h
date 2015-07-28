//
//  intra_potential.h
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#ifndef ____intra_potential__
#define ____intra_potential__


void intrapot_torsion_well(double (*boltz_ptr), double (*energy_ptr), int (*interaction_hist_ptr), double (*highest_factor_ptr), double (*well_ptr), double (*tor_ptr), double (nT_ptr), const unsigned int ntorangles, int (**poly_ptr), const unsigned int numberofbeads);

void intrapot_torsion_well_scan(double (*boltz_ptr), double (*energy_ptr), int (*interaction_hist_ptr), double (*highest_factor_ptr), double (*well_ptr), double (*tor_ptr), double (nT_ptr), const unsigned int ntorangles, int (**poly_ptr), const unsigned int numberofbeads);


void intrapot_torsion_well_test(double (*boltz_ptr), double (*energy_ptr), int (*interaction_hist_ptr), double(*extra_histogram), double (*highest_factor_ptr), double (*well_ptr), double (*tor_ptr), double (nT_ptr), const unsigned int ntorangles, int (**poly_ptr), const unsigned int numberofbeads);

void intra_binenergy(double energy, double min, double width, int (*histogram));

double intra_loss_of_conf(double (*sum_boltz_ptr), double highest_factor, unsigned long numberofframes);



#endif /* defined(____intra_potential__) */
