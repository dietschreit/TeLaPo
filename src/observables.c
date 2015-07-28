//
//  observables.c
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#include <stdlib.h>
#include <math.h> /* for math stuff like cos, log, etc. */
#include "main.h"
#include "observables.h"

//=========================================================================================================
// recasts the diamond lattice polymer with a lattice constant of sqrt(3) to any given bond length
void recast (double (**double_poly_ptr), int (**int_poly_ptr), unsigned int numberofbeads , double *recast_ptr) {
    
    for (int dim1=0; dim1<numberofbeads; dim1++){
        
        double_poly_ptr[dim1][0] = *recast_ptr * (double)int_poly_ptr[dim1][0];
        double_poly_ptr[dim1][1] = *recast_ptr * (double)int_poly_ptr[dim1][1];
        double_poly_ptr[dim1][2] = *recast_ptr * (double)int_poly_ptr[dim1][2];
    }
    
}



//=========================================================================================================
// returns a distance squared
double double_distance2 (double (*atom1_ptr), double (*atom2_ptr)) {
    
    return ((atom1_ptr[0]-atom2_ptr[0])*(atom1_ptr[0]-atom2_ptr[0]) + \
            (atom1_ptr[1]-atom2_ptr[1])*(atom1_ptr[1]-atom2_ptr[1]) + \
            (atom1_ptr[2]-atom2_ptr[2])*(atom1_ptr[2]-atom2_ptr[2]));
    
}

// returns a distance squared
int int_distance2 (int (*atom1_ptr), int (*atom2_ptr)) {
    
    return ((atom1_ptr[0]-atom2_ptr[0])*(atom1_ptr[0]-atom2_ptr[0]) + \
            (atom1_ptr[1]-atom2_ptr[1])*(atom1_ptr[1]-atom2_ptr[1]) + \
            (atom1_ptr[2]-atom2_ptr[2])*(atom1_ptr[2]-atom2_ptr[2]));
    
}


//=========================================================================================================
// returns the square of the radius of gyration of the conformation
double radius_of_gyration (double (**poly_ptr), unsigned int numberofbeads) {
    
    double *centroid;
    centroid = calloc(3, sizeof(double));
    
    double rgyr2 = 0.0;
    
    for (int dim1=0; dim1<(numberofbeads); dim1++){
        centroid[0] += poly_ptr[dim1][0];
        centroid[1] += poly_ptr[dim1][1];
        centroid[2] += poly_ptr[dim1][2];
    }
    
    // calculated the position of the centroid
    centroid[0] /= (double)(numberofbeads);
    centroid[1] /= (double)(numberofbeads);
    centroid[2] /= (double)(numberofbeads);
    
    for (int dim1=0; dim1<(numberofbeads); dim1++){
        rgyr2 += (poly_ptr[dim1][0]-centroid[0])*(poly_ptr[dim1][0]-centroid[0]) + \
        (poly_ptr[dim1][1]-centroid[1])*(poly_ptr[dim1][1]-centroid[1]) + \
        (poly_ptr[dim1][2]-centroid[2])*(poly_ptr[dim1][2]-centroid[2]);
        
    }
    
    // calculate mean of square
    rgyr2 /= (double)(numberofbeads);
    
    free(centroid);
    
    return rgyr2;
    
}



//=========================================================================================================
// returns torsion angle between the four consecutive atoms
// in the old programm this is much easier since there is this vector < Vec3 > "name" stuffe, which has .dot(), .corss(), .normalize() (struct Vec3), I will try to write something like this here, avoiding vector, which occupies more discspace than a normal array, it's probably also slower
double torsion_angle(double (**atoms_ptr), int start_position){
    
    // initialize bond vectors and fill them
    double **bond_vectors;
    bond_vectors = malloc(3 * sizeof(double *));
    
    // make the following explicit, because short for-loop are hard to optimize
    for (int dim1=0; dim1<3; dim1++){
        bond_vectors[dim1] = malloc(3 * sizeof(double *));
        
        for (int dim2=0; dim2<3; dim2++){// difference in x, y, z
            
            bond_vectors[dim1][dim2] = atoms_ptr[start_position + dim1 +1][dim2] - atoms_ptr[start_position + dim1][dim2];
            
        }
    }
    
    // bond normals
    double **bond_normals;
    bond_normals = malloc(3 * sizeof(double *));
    for (int dim1=0; dim1<3; dim1++){
        bond_normals[dim1] = malloc(3 * sizeof(double *));
    }
    
    // cross products
    bond_normals[0][0] = bond_vectors[0][1]*bond_vectors[1][2] - bond_vectors[0][2]*bond_vectors[1][1];
    bond_normals[0][1] = bond_vectors[0][2]*bond_vectors[1][0] - bond_vectors[0][0]*bond_vectors[1][2];
    bond_normals[0][2] = bond_vectors[0][0]*bond_vectors[1][1] - bond_vectors[0][1]*bond_vectors[1][0];
    
    bond_normals[1][0] = bond_vectors[1][1]*bond_vectors[2][2] - bond_vectors[1][2]*bond_vectors[2][1];
    bond_normals[1][1] = bond_vectors[1][2]*bond_vectors[2][0] - bond_vectors[1][0]*bond_vectors[2][2];
    bond_normals[1][2] = bond_vectors[1][0]*bond_vectors[2][1] - bond_vectors[1][1]*bond_vectors[2][0];
    
    
    //normalizing the two normal vectors and the second bond vector
    double vec_length;
    
    vec_length= sqrt( bond_normals[0][0]*bond_normals[0][0] + bond_normals[0][1]*bond_normals[0][1] + bond_normals[0][2]*bond_normals[0][2] );
    bond_normals[0][0]/=vec_length;
    bond_normals[0][1]/=vec_length;
    bond_normals[0][2]/=vec_length;
    
    vec_length= sqrt( bond_normals[1][0]*bond_normals[1][0] + bond_normals[1][1]*bond_normals[1][1] + bond_normals[1][2]*bond_normals[1][2] );
    bond_normals[1][0]/=vec_length;
    bond_normals[1][1]/=vec_length;
    bond_normals[1][2]/=vec_length;
    
    vec_length= sqrt( bond_vectors[1][0]*bond_vectors[1][0] + bond_vectors[1][1]*bond_vectors[1][1] + bond_vectors[1][2]*bond_vectors[1][2] );
    bond_vectors[1][0]/=vec_length;
    bond_vectors[1][1]/=vec_length;
    bond_vectors[1][2]/=vec_length;
    
    //this istn't a normal to two bonds but forms an orthogonal frame with 2nd bond and 1st bond normal
    //it gives the orientation around the bond to distinguish -60 from +60 degrees
    bond_normals[2][0] = bond_normals[0][1]*bond_vectors[1][2] - bond_normals[0][2]*bond_vectors[1][1];
    bond_normals[2][1] = bond_normals[0][2]*bond_vectors[1][0] - bond_normals[0][0]*bond_vectors[1][2];
    bond_normals[2][2] = bond_normals[0][0]*bond_vectors[1][1] - bond_normals[0][1]*bond_vectors[1][0];
    
    
    // these x and y components are needed for atan2
    double x_torsion = bond_normals[2][0]*bond_normals[1][0] + bond_normals[2][1]*bond_normals[1][1] + bond_normals[2][2]*bond_normals[1][2];
    double y_torsion = bond_normals[0][0]*bond_normals[1][0] + bond_normals[0][1]*bond_normals[1][1] + bond_normals[0][2]*bond_normals[1][2];
    
    double angle;
    angle = 180.0*atan2(x_torsion, y_torsion)/pi;
    
    // free the allocated memory
    free(bond_vectors[0]);
    free(bond_vectors[1]);
    free(bond_vectors[2]);
    free(bond_vectors);
    free(bond_normals[0]);
    free(bond_normals[1]);
    free(bond_normals[2]);
    free(bond_normals);
    
    return angle;
    
}


double get_nT(double (**atoms_ptr), unsigned int numtorsions){
    
    int count_T = 0;
    
    for (int dim1=0; dim1<(numtorsions); dim1++){
        
        if (fabs(torsion_angle(atoms_ptr, dim1))>120.0){
            count_T += 1;
        }
        
    }
    
    return ((double) count_T);
}


//=========================================================================================================
// returns the total number of attempts to estimate S
// this total numer includes all attempts to grow a chain, also those which are "hidden" in fragments or bricks

double recalc_attempts(unsigned long *attempts_successes, unsigned int *n_frag_brick, int n_bricks, int typeofrun){
    
    double attempts=1.0;
    
    if (typeofrun==10 || typeofrun==20 || typeofrun==30) {
        attempts = log((double)attempts_successes[0]);
    }
    else if (typeofrun==11 || typeofrun==21 || typeofrun==31) {
        attempts = log((double)attempts_successes[0]) + \
        (double)n_frag_brick[0]*log((double)attempts_successes[1]/(double)attempts_successes[2]) + \
        log((double)attempts_successes[3]/(double)attempts_successes[4]);
    }
    else if (typeofrun==12 || typeofrun==22 || typeofrun==32) {
        attempts = 0.0;
    }
    else if (typeofrun==13 || typeofrun==23 || typeofrun==33) {
        attempts = log((double)attempts_successes[0]) + \
        (double)n_frag_brick[3]*(log(4.0) + (double)(ARG_blength-1)*log(3.0) - log((double)n_bricks)) + \
        log((double)attempts_successes[1]/(double)attempts_successes[2]);
    }
    else if (typeofrun==14 || typeofrun==24 || typeofrun==34) {
        attempts = log((double)attempts_successes[0]) + \
        (double)(n_frag_brick[3]*n_frag_brick[0]+n_frag_brick[4])*(log(4.0) + (double)(ARG_blength-1)*log(3.0) - log((double)n_bricks)) + \
        (double)n_frag_brick[0]*log((double)attempts_successes[1]/(double)attempts_successes[2]) + \
        log((double)attempts_successes[3]/(double)attempts_successes[4]);
    }
    
    
    
    
    return attempts;
    
}



