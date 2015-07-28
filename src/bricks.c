//
//  bricks.c
//  
//
//  Created by Johannes Dietschreit on 29.05.15.
//
//

#include <stdlib.h>
#include <stdbool.h> /* defines things like true, false etc  */
#include <stdio.h>
#include "checks.h"
#include "bricks.h"



int ***make_bricks_saw (int (***bricks_ptr), const int (moves_ptr)[2][4][3], const int (*saw_ptr)[3], int blength_var, int *numbricks_ptr, int strictness) {
    
    int odd_even;
    
    // initialize
    int ***bricks_temp;
    bricks_temp = (int***) calloc(4, sizeof(int**)); //num of bricks, for the beginning
    for (int dim1=0; dim1<4; dim1++){
        // blength
        bricks_temp[dim1] = (int**) calloc((blength_var +1), sizeof(int*));
        for (int dim2=0; dim2<(blength_var + 1); dim2++){
            // x, y, z, lastmove
            bricks_temp[dim1][dim2] = (int*) calloc (4, sizeof(int));
        }
    }
    
    int *good_ones;
    
    int num_bricks;
    int new_bricks;
    
    
    // start of first loop
    for (int dim1=0; dim1<4; dim1++){// num of bricks
        for (int dim2=0; dim2<3; dim2++){ // x, y, z
            bricks_temp[dim1][1][dim2] = moves_ptr[1][dim1][dim2];
        }
        bricks_temp[dim1][1][3] = dim1; // lastmove
    }
    
    // this is only for the following to work
    new_bricks = 4;
    good_ones = (int*) malloc(new_bricks * sizeof(int));
    for (int dim1=0; dim1<new_bricks; dim1++){
        good_ones[dim1] = dim1;
    }
    
    
    // real loop
    for (int atom=2; atom<(blength_var + 1); atom++){
        odd_even = (atom%2);
        // hand over the counter of bricks, which obey the criterion
        num_bricks = new_bricks;
        new_bricks = 0;
        
        bricks_temp = (int ***) realloc(bricks_temp, (3*num_bricks)*sizeof(int**));
        if(NULL == bricks_temp) {
            fprintf(stderr, "No RAM for bricks \n");
            return NULL;
        }
        
        for (int dim1=(num_bricks-1); dim1>=0; dim1--){
            for (int dim11=2; dim11>=0; dim11--) {
                
                if ((3*dim1+dim11)>=num_bricks){
                    bricks_temp[3*dim1+dim11] = (int **) malloc((blength_var +1) * sizeof(int*));
                    if(NULL == bricks_temp[3*dim1+dim11]) {
                        fprintf(stderr, "No RAM for bricks \n");
                        return NULL;
                    }
                    for (int dim2=0; dim2<(blength_var +1); dim2++){
                        // x, y, z, lastmove
                        bricks_temp[3*dim1+dim11][dim2] = (int*) malloc(4 * sizeof(int));
                        if(NULL == bricks_temp[3*dim1+dim11][dim2]) {
                            fprintf(stderr, "No RAM for bricks \n");
                            return NULL;
                        }
                    }
                }
                
                // every entry gets tripled to add all three possible forward moves
                for (int dim2=0; dim2<atom; dim2++) {//copying
                    for (int dim3=0; dim3<4; dim3++){// x, y, z, lastmove
                        bricks_temp[3*dim1+dim11][dim2][dim3] = bricks_temp[good_ones[dim1]][dim2][dim3];
                    }
                }
            }
        }
        
        good_ones = (int*) malloc((3*num_bricks) * sizeof(int));// it has to be this long, since it is unknown how many are lost at each step
        if(NULL == good_ones) {
            fprintf(stderr, "No RAM for bricks \n");
            return NULL;
        }
        
        
        for (int dim1=0; dim1<(num_bricks); dim1++) { // get the new position of each atom
            for (int dim2=0; dim2<3; dim2++){ // saw_moves
                for (int dim3=0; dim3<3; dim3++) {// x, y, z
                    bricks_temp[dim1*3 + dim2][atom][dim3] = bricks_temp[dim1*3 + dim2][atom-1][dim3] + \
                    moves_ptr[odd_even][saw_ptr[bricks_temp[dim1*3 + dim2][atom-1][3]][dim2]][dim3];
                }
                bricks_temp[dim1*3 + dim2][atom][3] = saw_ptr[bricks_temp[dim1*3 + dim2][atom-1][3]][dim2]; //lastmove
                
                // check for self-avoidance
                switch (strictness) {
                    case 1:
                        if(check_saw1(bricks_temp[dim1*3 + dim2], atom)==false){
                            good_ones[new_bricks] = dim1*3 + dim2;
                            new_bricks+=1;
                        }
                        break;
                    case 2:
                        if(check_saw2(bricks_temp[dim1*3 + dim2], atom, bricks_temp[dim1*3 + dim2][atom][3], saw_ptr, moves_ptr, odd_even)==false){
                            good_ones[new_bricks] = dim1*3 + dim2;
                            new_bricks+=1;
                        }
                        break;
                    case 3:
                        if(check_saw3(bricks_temp[dim1*3 + dim2], atom, bricks_temp[dim1*3 + dim2][atom][3], saw_ptr, moves_ptr, odd_even)==false){
                            good_ones[new_bricks] = dim1*3 + dim2;
                            new_bricks+=1;
                        }
                        break;
                        
                    default:
                        break;
                }
                
            }
        }
        
        
    }
    
    
    // handing over the final ensemble to the real one
    bricks_ptr = (int***) calloc(new_bricks, sizeof(int**));
    if(NULL == bricks_ptr) {
        fprintf(stderr, "No RAM for bricks \n");
        return NULL;
    }
    
    for (int dim1=0; dim1<(new_bricks); dim1++){
        
        bricks_ptr[dim1] = (int **) calloc((blength_var), sizeof(int* ));
        if(NULL == bricks_ptr[dim1]) {
            fprintf(stderr, "No RAM for bricks \n");
            return NULL;
        }
        
        for (int dim2=0; dim2<(blength_var); dim2++) {
            
            bricks_ptr[dim1][dim2] = (int*) calloc(4, sizeof(int));
            if(NULL == bricks_ptr[dim1][dim2]) {
                fprintf(stderr, "No RAM for bricks \n");
                return NULL;
            }
            
            for (int dim3=0; dim3<4; dim3++) {
                // the real bricks don't need to start with (0, 0, 0)
                bricks_ptr[dim1][dim2][dim3] = bricks_temp[good_ones[dim1]][dim2+1][dim3];
                //                printf("%i ,", bricks_ptr[dim1][dim2][dim3]);
            }
            //            printf("\n");
        }
        //        printf("\n");
    }
    
    for (int dim1=0; dim1<num_bricks; dim1++){
        for (int dim2=0; dim2<blength_var; dim2++){
            free(bricks_temp[dim1][dim2]);
        }
        free(bricks_temp[dim1]);
    }
    free(bricks_temp);
    free(good_ones);
    
    *numbricks_ptr = new_bricks;
    
    // check that everything went ok
    return bricks_ptr;
    
}


