//
//  tetra_rw.c
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#include "random_knuth.h"
#include "tetra_rw.h"

void tetra_rw (int (**poly_ptr), const int (move_ptr)[2][4][3], unsigned int numberofbeads) {
    int odd, randn;
    
    for (int bead=1; bead<(numberofbeads); bead++){
        
        odd=bead%2;
        
        
        // random number generation
        do {
            randn = (int) 4.0*rand_knuth(&ARG_randomseed);
        } while (randn > 3);
        
        
        poly_ptr[bead][0] = poly_ptr[bead-1][0] + move_ptr[odd][randn][0];
        poly_ptr[bead][1] = poly_ptr[bead-1][1] + move_ptr[odd][randn][1];
        poly_ptr[bead][2] = poly_ptr[bead-1][2] + move_ptr[odd][randn][2];
        
    }
}

