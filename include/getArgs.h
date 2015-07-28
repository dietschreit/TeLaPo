//
//  getArgs.h
//  
//
//  Created by Johannes Dietschreit on 29.05.15.
//
//

#ifndef ____getArgs__
#define ____getArgs__

// pre-declarations
void usage_error(void);
void getArgs(FILE *output_ptr, int argc, char **argv);
void print_set_options(FILE *output_ptr, int typeofrun, int flength, int fflength, int blength, int numberofbeads, unsigned long long numberofframes, long seed, double bond, int torsion, int intramolecular, double intra_eps1[2], double intra_eps2[2]);


// initializing the global variables, which will be set in this routine
int ARG_typeofrun;      // sets method
int ARG_strictness;     // degree of self-avoidance
int ARG_blength;        // sets bricklength
int ARG_flength;        // sets fragment length

long ARG_randomseed;

unsigned int ARG_numberofbeads;     // number of atoms in chain
unsigned long ARG_numberofframes;    // number  of chains that will be generated successfully

//optional variables
double ARG_bondlength;  // bond length
int ARG_torsion;        // switch for torsional analysis
int ARG_pair_correlation;
int ARG_intra_potential; // switch for intra molecular potential
double ARG_intra_parameter1[2];
double ARG_intra_parameter2[2];


char* dcdFileName; // name of the dcd file which is to be loaded instead of a generation run


#endif /* defined(____getArgs__) */
