//
//  tetra_fww.c
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#include "random_knuth.h"
#include "tetra_fww.h"


void tetra_fww (int (**poly_ptr), const int (move_ptr)[2][4][3], const int (saw_ptr)[4][3], unsigned int numberofbeads) {
    
    static int odd, randn, lastmove;
    
    // random number generation
    do {
        randn = (int) 4.0*rand_knuth(&ARG_randomseed);
    } while (randn > 3);
    
    poly_ptr[1][0] = poly_ptr[0][0] + move_ptr[1][randn][0];
    poly_ptr[1][1] = poly_ptr[0][1] + move_ptr[1][randn][1];
    poly_ptr[1][2] = poly_ptr[0][2] + move_ptr[1][randn][2];
    
    lastmove = randn;
    
    for (int bead=2; bead<(numberofbeads); bead++){
        
        odd=bead%2;
        
        // random number generation
        do {
            randn = (int) 3.0*rand_knuth(&ARG_randomseed);
        } while (randn > 2);
        
        poly_ptr[bead][0] = poly_ptr[bead-1][0] + move_ptr[odd][saw_ptr[lastmove][randn]][0];
        poly_ptr[bead][1] = poly_ptr[bead-1][1] + move_ptr[odd][saw_ptr[lastmove][randn]][1];
        poly_ptr[bead][2] = poly_ptr[bead-1][2] + move_ptr[odd][saw_ptr[lastmove][randn]][2];
        
        lastmove = saw_ptr[lastmove][randn];;
        
    }
}
