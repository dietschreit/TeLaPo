//
//  observables.h
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#ifndef ____observables__
#define ____observables__

// global variables
int ARG_blength;
double pi;

typedef struct {
    double ee2;
    double maxee;
    double rgyr;
    double pt;
    double dsasa;
    double *pair_corr;
    double S;
    
    
    double *err_ee2;
    double *err_rgyr;
    double *err_pt;
    double *err_dsasa;
    double err_S;
    
    
}OBSERVABLE;


// predeclaration of functions
void recast (double (**double_poly_ptr), int (**int_poly_ptr), unsigned int numberofbeads , double *recast_ptr);

double double_distance2 (double (*atom1_ptr), double (*atom2_ptr));

int int_distance2 (int (*atom1_ptr), int (*atom2_ptr));

double radius_of_gyration (double (**poly_ptr), unsigned int numberofbeads);

double torsion_angle(double (**atoms_ptr), int start_position);

double get_nT(double (**atoms_ptr), unsigned int numtorsions);

double recalc_attempts(unsigned long *attempts_successes, unsigned int *n_frag_brick, int n_bricks, int typeofrun);






#endif /* defined(____observables__) */
