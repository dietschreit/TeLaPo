//
//  intra_potential.c
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#include <math.h>

#include "observables.h"
#include "intra_potential.h"


//-----------------------------------------------------------------------------
// nearest neighbor potential + torion potential



void intrapot_torsion_well(double (*boltz_ptr), double (*energy_ptr), int (*interaction_hist_ptr), double (*highest_factor_ptr), double (*well_ptr), double (*tor_ptr), double (nT_ptr), const unsigned int ntorangles, int (**poly_ptr), const unsigned int numberofbeads) {
    
    int interaction_counter = 0;
    int dist2;
    
    for (int dim1=0; dim1<(numberofbeads-8); dim1++) {
        for (int dim2=dim1+8; dim2<(numberofbeads); dim2++) {
            dist2 = int_distance2(poly_ptr[dim1], poly_ptr[dim2]);
            if (dist2<23){
                interaction_counter+=1;
            }
        }
    }
    
    // count interactions in histogram
    interaction_hist_ptr[interaction_counter]++;
    
    double energy;

    // first part transforms nT into nG which will be repulsive
    energy = ( ((double)ntorangles - nT_ptr) * tor_ptr[0]) - (double)interaction_counter*well_ptr[0];
            
    energy_ptr[0] = energy;
    boltz_ptr[0] = exp(-energy);
            
    // The value of the highest factor and the sum of all factors gives estimation for how many conformations are "lost"
    if (boltz_ptr[0]>highest_factor_ptr[0]) {
        highest_factor_ptr[0] = boltz_ptr[0];
    }
}



void intrapot_torsion_well_scan(double (*boltz_ptr), double (*energy_ptr), int (*interaction_hist_ptr), double (*highest_factor_ptr), double (*well_ptr), double (*tor_ptr), double (nT_ptr), const unsigned int ntorangles, int (**poly_ptr), const unsigned int numberofbeads) {
    
    int interaction_counter = 0;
    int dist2;
    int bin;
    
    for (int dim1=0; dim1<(numberofbeads-8); dim1++) {
        for (int dim2=dim1+8; dim2<(numberofbeads); dim2++) {
            dist2 = int_distance2(poly_ptr[dim1], poly_ptr[dim2]);
            if (dist2<23){
                interaction_counter+=1;
            }
        }
    }
    
    // count interactions in histogram
    interaction_hist_ptr[interaction_counter]++;
    
    
    double energy;
    
    for (int dim1=0; dim1<10; dim1++) {// 1-9 nn-well, indices run over the epilons
        for (int dim2=0; dim2<10; dim2++) {//torsion, indices run over the epilons
            // first part transforms nT into nG which will be repulsive
            energy = ( ((double)ntorangles - nT_ptr) * tor_ptr[dim2]) - (double)interaction_counter*well_ptr[dim1];
            
            bin = dim1*10 + dim2;
            
            energy_ptr[bin] = energy;
            boltz_ptr[bin] = exp(-energy);
            
            // The value of the highest factor and the sum of all factors gives estimation for how many conformations are "lost"
            if (boltz_ptr[bin]>highest_factor_ptr[bin]) {
                highest_factor_ptr[bin] = boltz_ptr[bin];
            }
        }
    }
    
    
}


// the lower function was used just to test something maybe it should be removed

void intrapot_torsion_well_test(double (*boltz_ptr), double (*energy_ptr), int (*interaction_hist_ptr), double(*extra_histogram), double (*highest_factor_ptr), double (*well_ptr), double (*tor_ptr), double (nT_ptr), const unsigned int ntorangles, int (**poly_ptr), const unsigned int numberofbeads) {
    
    int interaction_counter = 0;
    int dist2;
    int bin;
    int extra_counter = 0;
    
    
    for (int dim1=0; dim1<(numberofbeads-4); dim1++) {
        for (int dim2=dim1+4; dim2<(numberofbeads); dim2++) {
            dist2 = int_distance2(poly_ptr[dim1], poly_ptr[dim2]);
            if (dist2<23){
                interaction_counter+=1;
            }
        }
    }
    
    // count interactions in histogram
    interaction_hist_ptr[interaction_counter]++;
    
    
    double energy;
    
    for (int dim1=0; dim1<10; dim1++) {// 1-9 nn-well
        for (int dim2=0; dim2<10; dim2++) {//torsion
            // first part transforms nT into nG which will be repulsive
            energy = ( ((double)ntorangles - nT_ptr) * tor_ptr[dim2]) - (double)interaction_counter*well_ptr[dim1];
            
            bin = dim1*10 + dim2;
            
            energy_ptr[bin] = energy;
            boltz_ptr[bin] = exp(-energy);
            
            // The value of the highest factor and the sum of all factors gives estimation for how many conformations are "lost"
            if (boltz_ptr[bin]>highest_factor_ptr[bin]) {
                highest_factor_ptr[bin] = boltz_ptr[bin];
            }
        }
    }
    
    
    // this the part why this ist called test, this is not a regularly implemented method
    /*
    for (int dim2=0; dim2<30; dim2++){
        for (int dim1=(4+dim2); dim1<(numberofbeads-(4+dim2)); dim1++) {
            
            dist2 = int_distance2(poly_ptr[dim1], poly_ptr[dim1-(4+dim2)]);
            if (dist2<23){
                extra_counter+=1;
            }
            dist2 = int_distance2(poly_ptr[dim1], poly_ptr[dim1+(4+dim2)]);
            if (dist2<23){
                extra_counter+=1;
            }
            
            //extra_histogram[extra_counter+(dim2*3)]++;
            extra_histogram[extra_counter+(dim2*3)] += boltz_ptr[46];
            extra_counter = 0;
        }
        
        
    }

    */
    
    
    
}








void intra_binenergy(double energy, double min, double width, int (*histogram)){
    
    int bin;
    
    bin = floor((energy-min)/width);
    
    histogram[bin]++;
    
}


double intra_loss_of_conf(double (*sum_boltz_ptr), double highest_factor, unsigned long numberofframes){
    
    double sum = 0.0;
    
    sum += sum_boltz_ptr[0];
    sum += sum_boltz_ptr[1];
    sum += sum_boltz_ptr[2];
    sum += sum_boltz_ptr[3];
    sum += sum_boltz_ptr[4];
    sum += sum_boltz_ptr[5];
    sum += sum_boltz_ptr[6];
    sum += sum_boltz_ptr[7];
    sum += sum_boltz_ptr[8];
    sum += sum_boltz_ptr[9];
    
    sum /= (double)numberofframes;
    
    sum /= highest_factor;
    
    return sum;
    
}
