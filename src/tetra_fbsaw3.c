//
//  tetra_fbsaw3.c
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#include <stdio.h>
#include <stdbool.h>
#include "random_knuth.h"
#include "checks.h"

#include "tetra_fbsaw3.h"



void tetra_fb_saw3 (int (**poly_ptr), const int (move_ptr)[2][4][3], const int (saw_ptr)[4][3], unsigned int numberofbeads, int blength_var, int numbricks_var, int (***bricks_ptr), int flength, unsigned long *tries) {
    
    // call all the variables needed in this function
    int odd, randn, lastmove, lastmove_lastf;
    int insert_start_f, insert_start_b, position;
    
    // these values will be used in every loop and don't change within one run of the programm
    static int nfrags = 0;
    static int nbricks_frag = 0;
    static int nbricks_end = 0;
    
    // determine some important values
    if (nfrags==0 && nbricks_frag==0 && nbricks_end==0){
        
        nfrags = ((numberofbeads-1)/flength);
        nbricks_frag = (flength/blength_var);
        nbricks_end = ((numberofbeads-nfrags*flength-1)/blength_var);
        
        printf("nfrags = %i; nbricks_frag = %i; nbricks_end = %i \n", nfrags, nbricks_frag, nbricks_end);
    }
    
    //--------------------------------------------------------------------------
    
    // set everything for a completely fresh run
try_fb_saw3:
    
    // count attempt for the whole chain
    tries[0] += 1;
    
    lastmove_lastf = -1;
    insert_start_f = 0;
    
    
    for (int fragment=0; fragment<(nfrags); ++fragment) {
        
    try_frag_fb_saw3:
        
        // attempt long fragment
        tries[1] += 1;
        
        lastmove = lastmove_lastf;
        insert_start_b = insert_start_f;
        
        
        for (int brick=0; brick<nbricks_frag; ++brick){
            
            //            insert_start_b = insert_start_f + brick*blength_var;
            
        new_fbrick_fb_saw3:
            
            // random number generation
            do {
                randn = (int) ((double)numbricks_var*rand_knuth(&ARG_randomseed));
            } while (randn > (numbricks_var-1));
            
            
            if (bricks_ptr[randn][0][3]==lastmove){
                goto new_fbrick_fb_saw3;
            }
            
            position = insert_start_b;
            
            for (int atom=0; atom<blength_var; ++atom) {
                
                position += 1;
                odd = position%2;
                
                poly_ptr[position][0] = poly_ptr[insert_start_b][0] + bricks_ptr[randn][atom][0];
                poly_ptr[position][1] = poly_ptr[insert_start_b][1] + bricks_ptr[randn][atom][1];
                poly_ptr[position][2] = poly_ptr[insert_start_b][2] + bricks_ptr[randn][atom][2];
                
                // look at this checking procedure, it probably doesn't suit fb_saw3
                if (check_fb_saw3(poly_ptr, position, insert_start_f, bricks_ptr[randn][atom][3], saw_ptr, move_ptr, odd, atom)==true){
                    goto try_frag_fb_saw3;
                }
                
            }
            
            // this is the direction of the last bond of the latest brick
            lastmove = bricks_ptr[randn][blength_var-1][3];
            insert_start_b += blength_var;
        }
        
        //        insert_start_b += blength_var;
        
        // atoms within the fragment which are not put together by bricks
        for (int bead=(insert_start_b+1); bead<=(insert_start_f+flength); ++bead){
            
            odd=bead%2;
            
            // random number generation
            do {
                randn = (int) 3.0*rand_knuth(&ARG_randomseed);
            } while (randn > 2);
            
            
            poly_ptr[bead][0] = poly_ptr[bead-1][0] + move_ptr[odd][saw_ptr[lastmove][randn]][0];
            poly_ptr[bead][1] = poly_ptr[bead-1][1] + move_ptr[odd][saw_ptr[lastmove][randn]][1];
            poly_ptr[bead][2] = poly_ptr[bead-1][2] + move_ptr[odd][saw_ptr[lastmove][randn]][2];
            
            lastmove = saw_ptr[lastmove][randn];
            
            if (check_fsaw3(poly_ptr, bead, insert_start_f, lastmove, saw_ptr, move_ptr, odd)==true){
                goto try_frag_fb_saw3; // start again with this very fragment
            }
            
            
        }
        
        // fragment completed
        tries[2] += 1;
        
        // now the fragment itself is completed, now is the time to check it against the previous chain
        if (check_whole_fsaw3(poly_ptr, insert_start_f, (insert_start_f+flength), move_ptr)==true){
            goto try_fb_saw3; // start all over again
        }
        
        // the current fragment was added succesfully to the chain
        
        lastmove_lastf = lastmove;
        insert_start_f += flength;
        
        
    }// end of fragments
    
    
    
    // assembling the last fragment/part of the chain
    
try_end_fb_saw3:
    
    // attempt end of chain
    tries[3] += 1;
    
    lastmove = lastmove_lastf;
    insert_start_b = insert_start_f;
    
    
    for (int brick=0; brick<nbricks_end; ++brick){
        
    new_endbrick_fb_saw3:
        
        // random number generation
        do {
            randn = (int) ((double)numbricks_var*rand_knuth(&ARG_randomseed));
        } while (randn > (numbricks_var-1));
        
        
        if (bricks_ptr[randn][0][3]==lastmove){
            goto new_endbrick_fb_saw3;
        }
        
        position = insert_start_b;
        for (int atom=0; atom<blength_var; ++atom) {
            
            position += 1;
            odd = position%2;
            
            poly_ptr[position][0] = poly_ptr[insert_start_b][0] + bricks_ptr[randn][atom][0];
            poly_ptr[position][1] = poly_ptr[insert_start_b][1] + bricks_ptr[randn][atom][1];
            poly_ptr[position][2] = poly_ptr[insert_start_b][2] + bricks_ptr[randn][atom][2];
            
            // look at this checking procedure, it probably doesn't suit fb_saw3
            if (check_fb_saw3(poly_ptr, position, insert_start_f, bricks_ptr[randn][atom][3], saw_ptr, move_ptr, odd, atom)==true){
                goto try_end_fb_saw3;
            }
            
        }
        
        // this is the direction of the last bond of the latest brick
        lastmove = bricks_ptr[randn][blength_var-1][3];
        insert_start_b += blength_var;
    }
    
    
    
    // atoms within the fragment which are not put together by bricks
    for (int bead=(insert_start_b+1); bead<(numberofbeads); ++bead){
        
        odd=bead%2;
        
        // random number generation
        do {
            randn = (int) 3.0*rand_knuth(&ARG_randomseed);
        } while (randn > 2);
        
        
        poly_ptr[bead][0] = poly_ptr[bead-1][0] + move_ptr[odd][saw_ptr[lastmove][randn]][0];
        poly_ptr[bead][1] = poly_ptr[bead-1][1] + move_ptr[odd][saw_ptr[lastmove][randn]][1];
        poly_ptr[bead][2] = poly_ptr[bead-1][2] + move_ptr[odd][saw_ptr[lastmove][randn]][2];
        
        lastmove = saw_ptr[lastmove][randn];
        
        if (check_fsaw3(poly_ptr, bead, insert_start_f, lastmove, saw_ptr, move_ptr, odd)==true){
            goto try_end_fb_saw3; // start again with this very fragment
        }
        
        
    }
    
    // completed the end
    tries[4] += 1;
    
    // now the fragment itself is completed, now is the time to check it against the previous chain
    if (check_whole_fsaw3(poly_ptr, insert_start_f, (numberofbeads-1), move_ptr)==true){
        goto try_fb_saw3; // start all over again
    }
    
    
    //done
    
}





