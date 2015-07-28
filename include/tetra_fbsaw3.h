//
//  tetra_fbsaw3.h
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#ifndef ____tetra_fbsaw3__
#define ____tetra_fbsaw3__

long ARG_randomseed;

void tetra_fb_saw3 (int (**poly_ptr), const int (move_ptr)[2][4][3], const int (saw_ptr)[4][3], unsigned int numberofbeads, int blength_var, int numbricks_var, int (***bricks_ptr), int flength, unsigned long *tries);

#endif /* defined(____tetra_fbsaw3__) */
