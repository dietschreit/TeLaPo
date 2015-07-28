//
//  random_knuth.c
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include "random_knuth.h"


/*
 From
 Numerical Recepies, The Art of Scientific Computing
 Press, Teukolsky, Vetterling, Flannery
 Second Edition
 Cambridge University Press, 1992
 
 The primary reference is:
 Knuth, Donald E. 1981, Seminumerical Algorithms, 2nd ed., vol 2 of The Art of Computer Programming (Reading, Mass: Addison-Wesle), §§3.2-3.3
 
 
 According to Knuth, any large MBIG, and any smaller (but still large) MSEED can be substituted
 for the above values.
 
 Returns a uniform random deviate between 0.0 and 1.0. = [0,1], not [0, 1).
 Set idum to any negative value to initialize or reinitialize the sequence.
 */

float rand_knuth(long *idum){
    
    static int inext,inextp;
    static long ma[56];
    /* The value 56 (range ma[1..55]) is special and
     should not be modified; see Knuth. */
    static int iff=0;
    long mj,mk;
    int i,ii,k;
    
    if (*idum < 0 || iff == 0) { // Initialization.
        //printf("Random Number Generator Initiliazaiton. \n\n");
        iff=1;
        mj=labs(MSEED-labs(*idum));  // Initialize ma[55] using the seed idum and the large number MSEED.
        mj %= MBIG;
        ma[55]=mj;
        mk=1;
        for (i=1; i<=54; i++) {
            /* Now initialize the rest of the table,
             in a slightly random order,
             with numbers that are not especially random. */
            
            ii=(21*i) % 55;
            ma[ii]=mk;
            mk=mj-mk;
            
            if (mk < MZ) {
                mk += MBIG;
            }
            
            mj=ma[ii];
        }
        // We randomize them by “warming upthe generator.”
        for (k=1; k<=4;k++){
            for (i=1; i<=55; i++) {
                ma[i] -= ma[1+(i+30) % 55];
                if (ma[i] < MZ) {
                    ma[i] += MBIG;
                }
            }
        }
        
        /*
         Prepare indices for our first generated number. The constant 31 is special; see Knuth.
         */
        
        inext=0;
        inextp=31;
        *idum=1;
    }
    
    // Here is where we start, except on initialization.
    
    if (++inext == 56) {inext=1;} // Increment inext and inextp, wrapping around 56 to 1.
    if (++inextp == 56) {inextp=1;}
    
    mj=ma[inext]-ma[inextp]; //Generate a new random number subtractively.
    if (mj < MZ) {mj += MBIG;} // Be sure that it is in range.
    
    ma[inext]=mj; //Store it, and output the derived uniform deviate.
    
    return mj*FAC;
}

