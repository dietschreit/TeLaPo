//
//  main.c
//  
//
//  Created by Johannes Dietschreit on 29.05.15.
//
//

// general libraries
#include <stdio.h>  /* standard input/output like printf  */
#include <stdlib.h> /* standard libraries  */
#include <stdbool.h> /* defines things like true, false etc  */
#include <string.h> /* standard library for string functions, like strcmp */
#include <math.h> /* standard math stuff like cos, acos, etc  */

//-----------------------------------------------
// self-written files

// general functions
#include "output.h"
#include "getArgs.h"
#include "read_dcd.h"
#include "time.h"
#include "observables.h"
#include "statistics.h"
#include "entropy.h"
#include "pair_correlation_function.h"

// functions that are particularly involved in chain creation
#include "bricks.h"
#include "tetra_rw.h"
#include "tetra_fww.h"
#include "tetra_saw1.h"
#include "tetra_saw2.h"
#include "tetra_saw3.h"
#include "tetra_fsaw3.h"
#include "tetra_bsaw3.h"
#include "tetra_fbsaw3.h"

#include "intra_potential.h"
#include "SASA.h"

#include "main.h"




/* define important variables */
#define RANDMAX_4 ((RAND_MAX / 4 ) * 4)
#define RANDMAX_3 ((RAND_MAX / 3 ) * 3)

//---------------------------------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------------------------------

/*
 global variables, mostly defined by command line
 the settings are mostly default values
 */

// obligatory (always and addtional to some methods)
int ARG_typeofrun = -1;             // which method
int ARG_strictness = 0;             // used for various runs with building bricks to distinguish between saw1, 2, 3 ... and other checks
int ARG_blength = 1;                // number of bonds in of bricks
int ARG_flength = 1;                // number of bonds in fragments
int ARG_fflength = 1;               // number of bonds in small fragments (not yet implemented in this version)

unsigned int ARG_numberofbeads = 1; // how many atoms
unsigned long ARG_numberofframes = 10000;         // how many chains, default 10^4

long ARG_randomseed = 0;            // seed for the random number generator, if none is given, it will be taken from the system clock

char* dcdFileName;                  // name of the dcd file which is to be loaded instead of a generation run


// optional
double ARG_bondlength = 1;          // bond length
int ARG_torsion = 0;                // torsional analysis, default off
int ARG_pair_correlation = 0;       // switches on the calculatio of pair correlation
int ARG_intra_potential = 0;        // intramolecular potentials, default off
double ARG_intra_parameter1[2] = {0.0, 0.0};     // paramters
double ARG_intra_parameter2[2] = {0.0, 0.0};
int ARG_sasa = 0;                    // switches on the computation of SASA computation


//--------------------------------------------------------------
// physical constants
double pi;







int main(int argc, char **argv) {
    
    /* the times recorded are: 
     start of programm, start of chain generation,
     end of chain generation, end of programm */
    double time[4];
    time[0] = time_of_day();
    
    
    //-----------------------------------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------
    // open all neccessary files
    
    // general output file, contains basically everything which was written to the terminal
    FILE * output_file;
    output_file = fopen("output_file.dat", "w+");
    
    FILE * conv_obs_file;
    conv_obs_file = fopen("convergence_observables.dat", "w+");
    
    // pair correlation function can be weighted with potentials
    FILE * pair_correl_file;
    pair_correl_file = fopen("pair_correlation_function.dat", "w+");
    
    //------------------------------------------------------------------------------
    // all files connected to intramolecular interactions

    FILE * intra_pot_file;
    intra_pot_file = fopen("intramolecular_potential_file.dat", "w+");
    
    FILE * intra_boltzman_factors_hist;
    intra_boltzman_factors_hist = fopen("intramolecular_factors_hist.dat", "w+");
    
    FILE * intra_interactions_file;
    intra_interactions_file = fopen("intramolecular_interactions_hist.dat", "w+");
    
 //   FILE * intra_interactions_testfile;
 //   intra_interactions_testfile = fopen("intramolecular_interactions_testhist.dat", "w+");
    
    FILE * convergence_intraweights;
    convergence_intraweights = fopen("intramolecular_convergence_Z.dat", "w+");
    fprintf(convergence_intraweights, "### convergence_point -- sum-of-boltzman-factors \n");
    
    FILE * conv_intraobs_file;
    conv_intraobs_file = fopen("intramolecular_convergence_obs.dat", "w+");

    
    //-----------------------------------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------
    // input and first output
    
    
    // hello message of the programm, always displayed!
    
    log_out(output_file, "\n-----------------------------------------------------------------------------------\n");
    log_out(output_file, "TeLaPo , version 1.1\n");
    log_out(output_file, "\tQuestions and bug reports to: dietschreit@zedat.fu-berlin.de\n");
    log_out(output_file, "\tPlease include in published work based on TeLaPo:\n");
    log_out(output_file, "\t\tDietschreit, J. C. B.; Diestler, D. J.; Knapp, E. W.,\n");
    log_out(output_file, "\t\tModels for Self-Avoiding Polymer Chains on the Tetrahedral Lattice.\n");
    log_out(output_file, "\t\tMacromol. Theory Simul. 2014, 23, 452-463\n");
    log_out(output_file, "\tFor more information use the flag -info.\n");
    log_out(output_file, "-----------------------------------------------------------------------------------\n");
    
    
    
    /* get the arguments from the comand line */
    getArgs(output_file, argc, argv);
    
    
    //------------------------------------------------------------------------------------
    // for reading DCD-files
    // this has to be early in the code, because it sets the variables ARG_numberofbeads and ARG_numberofframes!
    
    // variables concerning reading dcd
    molfile_timestep_t timestep;
    void *v;
    dcdhandle *dcd;
    int natoms;
    float sizeMB =0.0, totalMB = 0.0;
    // reading the dcd-file and setting global variables accordingly
    if (ARG_typeofrun==0) {
        natoms = 0;
        v = open_dcd_read(dcdFileName, "dcd", &natoms);
        if (!v) {
            fprintf(stderr, "ERROR: open_dcd_read failed for file %s\n", dcdFileName);
            return EXIT_FAILURE;
        }
        
        dcd = (dcdhandle *)v;
        sizeMB = ((natoms * 3.0) * dcd->nsets * 4.0) / (1024.0 * 1024.0);
        totalMB += sizeMB;
        
        log_out(output_file, "Read DCD: %d atoms, %d frames, size: %6.1fMB\n", natoms, dcd->nsets, sizeMB);
        
        timestep.coords = (float *)malloc(3*sizeof(float)*natoms);
        
        ARG_numberofbeads=dcd->natoms;
        ARG_numberofframes=dcd->nsets;
    }
    //------------------------------------------------------------------------------------
    

    // print all the options to the screen so one can check whether the right thing gets computed
    print_set_options(output_file, ARG_typeofrun, ARG_flength, ARG_fflength, ARG_blength, ARG_numberofbeads, ARG_numberofframes, ARG_randomseed, ARG_bondlength, ARG_torsion, ARG_intra_potential, ARG_intra_parameter1, ARG_intra_parameter2);
    
    
    
    //-----------------------------------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------
    // physical constants
    
    pi = acos(-1.0);
    
    //-----------------------------------------------------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------------------------------------------------
    
    
    /* Initialize the most important variables */
    
    // stuff with bond lengths
    const double inv_sqrt3 = 1.0/sqrt(3.0);
    const double bondlength = ARG_bondlength;
    double recast_factor;
    if (ARG_typeofrun<40){// this recasts walk on diamond lattice
        recast_factor = inv_sqrt3*bondlength;
    }
    else {// this is for runs on simple cubic lattice
        recast_factor = bondlength;
    }
    
    // variables with atom numbers etc
    const int last_atom = ARG_numberofbeads -1;
    const unsigned int number_of_torsions = ARG_numberofbeads -3;
    // number of frag, fragfags, endfrags, fragbricks, endbricks
    unsigned int numof_frags_bricks[5] = {0};
    
    
    // fractions of the number of frames
    const unsigned long permill_frames = ARG_numberofframes / 1000; // used for convergence
    const unsigned long percent_frames = ARG_numberofframes / 100; // used for ramining time
    const unsigned long tenth_frames = ARG_numberofframes / 10; // used for error estimation
    int cent;
    int tenth;
    
    int counter; // counter which can be used at any parts of the main programm, should only be used locally in a loop
    
    
    // basic moves on the tetrahedral lattice, back and forth
    const int move[2][4][3] = {
        {
            {-1, -1, -1},
            {1, 1, -1},
            {1, -1, 1},
            {-1, 1, 1}
        },
        {
            {1, 1, 1},
            {-1, -1, 1},
            {-1, 1, -1},
            {1, -1, -1}
        }
    };
    
    // moves possible in SAW (no walking back)
    const int sawmoves[4][3] = {
        {1, 2, 3},
        {0, 2, 3},
        {0, 1, 3},
        {0, 1, 2}
    };
    
    
    //-----------------------------------------------------------------------
    // building bricks are used to put parts together which are pre-checked
    int ***building_bricks=NULL;
    
    int numberofbricks;
    
    if (ARG_typeofrun==13 || ARG_typeofrun==14 || ARG_typeofrun==23 || ARG_typeofrun==24 || ARG_typeofrun==33 || ARG_typeofrun==34){
        
        // thise generates the bricks
        building_bricks = make_bricks_saw(building_bricks, move, sawmoves, ARG_blength, &numberofbricks, ARG_strictness);
        if (NULL==building_bricks[0][0]){
            return EXIT_FAILURE;
        }
        
        log_out(output_file, "%d bricks were generated, with a length of %d \n", numberofbricks, ARG_blength);
        
    }
    
    
    
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    /* Initialize beads vector and set all values to zero */
    int **int_polymer;
    int_polymer = calloc(ARG_numberofbeads, sizeof(int *));
    
    double **double_polymer;
    double_polymer = calloc(ARG_numberofbeads, sizeof(double *));
    
    for (int dim1=0; dim1<ARG_numberofbeads;dim1++){
        
        int_polymer[dim1] = calloc(3, sizeof(int *));
        double_polymer[dim1] = calloc(3, sizeof(double *));
        
    }
    
    
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    // observables
    
    OBSERVABLE normal_obs;
    
    normal_obs.maxee = 0.0; // maximal stretch of polymer
    
    OBSERVABLE intra_obs[100];

    //-------------------------------------------------------
    // initialization of all the OBSERVABLE variables
    //-------------------------------------------------------
    
    normal_obs.err_ee2 = (double *) calloc(11, sizeof(double));
    if(NULL == normal_obs.err_ee2) {
        fprintf(stderr, "Allocation of ee2 variable failed! \n");
        return EXIT_FAILURE;
    }
    normal_obs.err_rgyr = (double *) calloc(11, sizeof(double));
    if(NULL == normal_obs.err_rgyr) {
        fprintf(stderr, "Allocation of rgyr variable failed! \n");
        return EXIT_FAILURE;
    }
    
    if (ARG_intra_potential>0) {
        for (int dim1=0; dim1<100; dim1++) {
            intra_obs[dim1].err_ee2 =(double *) calloc(11, sizeof(double));
            intra_obs[dim1].err_rgyr =(double *) calloc(11, sizeof(double));
            if(NULL == intra_obs[dim1].err_ee2 || NULL == intra_obs[dim1].err_rgyr) {
                fprintf(stderr, "Allocation of ee2 or rgyr variable (intramolecular) failed! \n");
                return EXIT_FAILURE;
            }
        }
    }
    
    
    // initializes the observables for torsional analysis (optional)
    if (ARG_torsion==1) {
        normal_obs.err_pt = (double *) calloc(11, sizeof(double));
        if(NULL == normal_obs.err_pt) {
            fprintf(stderr, "Allocation of torsion variable failed! \n");
            return EXIT_FAILURE;
        }
        
        if (ARG_intra_potential>0) {
            for (int dim1=0; dim1<100; dim1++) {
                intra_obs[dim1].pt = 0.0;
                intra_obs[dim1].err_pt = (double *) calloc(11, sizeof(double));
                if(NULL == intra_obs[dim1].err_pt) {
                    fprintf(stderr, "Allocation of torsion variable (intramolecular) failed! \n");
                    return EXIT_FAILURE;
                }
            }
        }
    }
    
    // initializes the observables for loss of solven accessible surface area analysis (optional)
    if (ARG_sasa==1){
        normal_obs.err_dsasa = (double *) calloc(11, sizeof(double));
        if(NULL == normal_obs.err_dsasa) {
            fprintf(stderr, "Allocation of D-SASA variable failed! \n");
            return EXIT_FAILURE;
        }
        
        if (ARG_intra_potential>0) {
            for (int dim1=0; dim1<100; dim1++) {
                intra_obs[dim1].dsasa = 0.0;
                intra_obs[dim1].err_dsasa = (double *) calloc(11, sizeof(double));
                if(NULL == intra_obs[dim1].err_dsasa) {
                    fprintf(stderr, "Allocation of D-SASA variable (intramolecular) failed! \n");
                    return EXIT_FAILURE;
                }
            }
        }
    }
    
    // this is needed for the pair correlation function
    double *pair_correlation_obs;
    if (ARG_pair_correlation==1) {
        pair_correlation_obs = (double*) calloc(2*ARG_numberofbeads, sizeof(double));
        if(NULL == pair_correlation_obs) {
            fprintf(stderr, "Allocation of pair_correlation_obs failed! \n");
            return EXIT_FAILURE;
        }
        
        normal_obs.pair_corr = (double *) calloc(2*ARG_numberofbeads, sizeof(double));
        if(NULL == normal_obs.pair_corr) {
            fprintf(stderr, "Allocation of pair_corr failed! \n");
            return EXIT_FAILURE;
        }
        
        if (ARG_intra_potential>0) {
            for (int dim1=0; dim1<100; dim1++) {
                intra_obs[dim1].pair_corr = (double *) calloc(2*ARG_numberofbeads, sizeof(double));
                if(NULL == intra_obs[dim1].pair_corr) {
                    fprintf(stderr, "Allocation of pair_corr (intramolecular) failed! \n");
                    return EXIT_FAILURE;
                }
            }
        }
    }
    //-------------------------
    
    
    // this will provide a measure for entropy loss calculation
    unsigned long *attempts_successes;
    double *log_attempts;
    log_attempts = (double*) calloc(11, sizeof(double));
    
    // set the number of entries in this list, it depends on the typeofrun, but not the saw-type
    // last entry is the recast number of attempts which provides a measure for the entropy loss
    
    //-------------------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------------------
    // calculate variables which depend on the typeofrun
    
    
    switch (ARG_typeofrun) {
            // normal SAWX
        case 10:
        case 20:
        case 30:
            attempts_successes = calloc(2, sizeof(unsigned long));
            break;
            // fSAWX
        case 11:
        case 21:
        case 31:
            attempts_successes = calloc(5, sizeof(unsigned long));
            numof_frags_bricks[0] = (ARG_numberofbeads-1)/ARG_flength;
            break;
            // bSAWX
        case 13:
        case 23:
        case 33:
            attempts_successes = calloc(3, sizeof(unsigned long));
            numof_frags_bricks[3] = (ARG_numberofbeads-1)/ARG_blength;
            break;
            // fb_SAWX
        case 14:
        case 24:
        case 34:
            attempts_successes = calloc(5, sizeof(unsigned long));
            numof_frags_bricks[0] = (ARG_numberofbeads-1)/ARG_flength;
            numof_frags_bricks[3] = ARG_flength/ARG_blength;
            numof_frags_bricks[4] = (ARG_numberofbeads-1-ARG_flength*numof_frags_bricks[3])/ARG_blength;
            break;
            
        default:
            break;
    }
    
    
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    // intra molecular interactions
    
    
    // boltzman-factors, the intramolecular energy, the enrgy time the boltzmanfactor (entropy)
    double *intra_boltzman_factors;
    double *intra_highest_boltzmanfactor;
    double *intra_energy;
    double **intra_sum_of_boltzfactors;
    double **intra_sum_of_enrgyboltz;
    int *intra_interactions_hist;
 //   double *intra_interactions_testhist;
    
    double intra_max_factor = 0.0;
    double intra_min_factor = 0.0;
    
    // binning the energies of the intra-factors
    double intra_binmin;
    double intra_binmax;
    double intra_binwidth;
    int **intra_energybin;
    
    
    
    // these will be the parameter of the intra molecular force
    double *intra_parameter1;
    double *intra_parameter2;
    
    
    
    if (ARG_intra_potential > 0){
        switch (ARG_intra_potential) {
                // 1-10 potential well + torsional potential, given energy values
            case 1:
                intra_boltzman_factors = (double*) calloc(1, sizeof(double));
                intra_highest_boltzmanfactor = (double*) calloc(1, sizeof(double));
                intra_energy = (double*) calloc(1, sizeof(double));
                intra_sum_of_boltzfactors = (double**) calloc(1, sizeof(double));
                intra_sum_of_enrgyboltz = (double**) calloc(1, sizeof(double*));
                intra_sum_of_boltzfactors[0] = (double*) calloc(10, sizeof(double));
                intra_sum_of_enrgyboltz[0]  = (double*) calloc(10, sizeof(double));
                
                intra_parameter1 = (double*) malloc(1 * sizeof(double));
                intra_parameter2 = (double*) malloc(1 * sizeof(double));
                intra_interactions_hist = calloc(ARG_numberofbeads, sizeof(int));
 //               intra_interactions_testhist = calloc(100, sizeof(double));
                intra_parameter1[0] = ARG_intra_parameter1[0]; // nearest neighbor potential
                intra_parameter2[0] = ARG_intra_parameter2[0]; // torsion potential
                intra_binmin = -100.0;
                intra_binmax = 100.0;
                intra_binwidth = 1.0;
                intra_energybin = (int**) calloc(1, sizeof(int *));
                intra_energybin[0] = (int*) calloc(((intra_binmax-intra_binmin)/intra_binwidth), sizeof(int));
                break;
                
                // 1-10 potential well + torsional potential, given energy value range! 10x10
            case 2:
                intra_boltzman_factors = (double*) calloc(100, sizeof(double));
                intra_highest_boltzmanfactor = (double*) calloc(100, sizeof(double));
                intra_energy = (double*) calloc(100, sizeof(double));
                intra_sum_of_boltzfactors = (double**) calloc(100, sizeof(double));
                intra_sum_of_enrgyboltz = (double**) calloc(100, sizeof(double*));
                
                for (int dim1=0; dim1<100; ++dim1) {
                    intra_sum_of_boltzfactors[dim1] = (double*) calloc(10, sizeof(double));
                    intra_sum_of_enrgyboltz[dim1]  = (double*) calloc(10, sizeof(double));
                }
                
                intra_parameter1 = (double*) malloc(10 * sizeof(double));
                intra_parameter2 = (double*) malloc(10 * sizeof(double));
                intra_interactions_hist = calloc(ARG_numberofbeads, sizeof(int));
 //               intra_interactions_testhist = calloc(100, sizeof(double));
                for (int dim1=0; dim1<10; dim1++) {
                    intra_parameter1[dim1] = ARG_intra_parameter1[0]+ ARG_intra_parameter1[1]*(double)dim1; // nearest neighbor potential
                    intra_parameter2[dim1] = ARG_intra_parameter2[0]+ ARG_intra_parameter2[1]*(double)dim1; // torsion potential
                }
                intra_binmin = -100.0;
                intra_binmax = 100.0;
                intra_binwidth = 1.0;
                intra_energybin = (int**) calloc(100, sizeof(int *));
                for (int dim1=0; dim1<100; ++dim1){
                    intra_energybin[dim1] = (int*) calloc(((intra_binmax-intra_binmin)/intra_binwidth), sizeof(int));
                }
                break;
                
            default:
                fprintf(stderr, "This intramolecular potential doesn't exist! \n");
                break;
        }

    }

    
    
//----------------------------------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------------------------------

    
    
    
    
    /* start of chain generation*/
    time[1] = time_of_day();
    
    for (unsigned long frame=0; frame<ARG_numberofframes; frame++){
        
        tenth = frame/tenth_frames;
        
        switch(ARG_typeofrun){
                
            case 0:
                dcd_to_polymer(double_polymer, v, natoms, &timestep, dcd, frame, ARG_numberofbeads);
                break;
                
            case 1:
                tetra_rw(int_polymer, move, ARG_numberofbeads);
                break;
                
            case 2:
                tetra_fww(int_polymer, move, sawmoves, ARG_numberofbeads);
                break;
                
            case 10:
                tetra_saw1(int_polymer, move, sawmoves, ARG_numberofbeads, attempts_successes);
                break;
                
            case 20:
                tetra_saw2(int_polymer, move, sawmoves, ARG_numberofbeads, attempts_successes);
                break;
                
            case 30:
                tetra_saw3(int_polymer, move, sawmoves, ARG_numberofbeads, attempts_successes);
                break;

            case 31:
                tetra_fsaw3(int_polymer, move, sawmoves, ARG_numberofbeads, ARG_flength, attempts_successes);
                break;
                
            case 33:
                tetra_bsaw3(int_polymer, move, sawmoves, ARG_numberofbeads, ARG_blength, numberofbricks, building_bricks, attempts_successes);
                break;
                
            case 34:
                tetra_fb_saw3(int_polymer, move, sawmoves, ARG_numberofbeads, ARG_blength, numberofbricks, building_bricks, ARG_flength, attempts_successes);
                break;

                
            default:
                log_out(output_file, "ERROR: ARG_typeofrun = %d isn't recognized by the main part of the programm!\n", ARG_typeofrun);
                usage_error();
                
        } // end of chain generaiton
        
        
        // copies the lattice polymer into an array with doubles so that the chosen bond length can be used.
        if (ARG_typeofrun>0) {
            recast(double_polymer, int_polymer, ARG_numberofbeads, &recast_factor);
        }
        
        
        //-----------------------------------------------------------------------------
        //-----------------------------------------------------------------------------
        // get most important observables
        
        // end-to-end distance^2
        normal_obs.ee2 = double_distance2(double_polymer[last_atom], double_polymer[0]);
        normal_obs.err_ee2[tenth] += normal_obs.ee2;
        // radius of gyration
        normal_obs.rgyr = radius_of_gyration(double_polymer, ARG_numberofbeads);
        normal_obs.err_rgyr[tenth] += normal_obs.rgyr;
        // pT
        if (ARG_torsion==1){
            normal_obs.pt = get_nT(double_polymer, number_of_torsions);
            normal_obs.err_pt[tenth] += normal_obs.pt;
        }
        if (ARG_sasa==1){
            //observables[4] = get_asa(double_polymer, ARG_numberofbeads, (sqrt(16.0/3.0)*bondlength/2.0));
            normal_obs.dsasa = delta_asa(double_polymer, ARG_numberofbeads, bondlength);
            normal_obs.err_dsasa[tenth] += normal_obs.dsasa;
        }
        if (ARG_pair_correlation==1) {
            if (false==pair_correlation_fct(pair_correlation_obs, double_polymer, ARG_numberofbeads)){
                return EXIT_FAILURE;
            }
            for (int pairs=0; pairs<(2*ARG_numberofbeads); pairs++) {
                normal_obs.pair_corr[pairs] += pair_correlation_obs[pairs];
            }
        }
        
        
        

        
        
        
        
        
        //-----------------------------------------------------------------------------
        //-----------------------------------------------------------------------------
        // calculate boltzman factors if intramolecular forces are switched on
        
        switch (ARG_intra_potential) {
                
            // torsion + square well potential
            case 1:
                intrapot_torsion_well(intra_boltzman_factors, intra_energy, intra_interactions_hist, intra_highest_boltzmanfactor, intra_parameter1, intra_parameter2, normal_obs.pt, number_of_torsions, int_polymer, ARG_numberofbeads);
                
                intra_sum_of_enrgyboltz[0][tenth] += (intra_boltzman_factors[0]*intra_energy[0]);
                intra_sum_of_boltzfactors[0][tenth] += intra_boltzman_factors[0];
                
                intra_obs[0].err_ee2[tenth] += normal_obs.ee2 * intra_boltzman_factors[0];
                intra_obs[0].err_rgyr[tenth] += normal_obs.rgyr * intra_boltzman_factors[0];
                intra_obs[0].err_pt[tenth] += normal_obs.pt * intra_boltzman_factors[0];
                
                intra_binenergy(intra_energy[0], intra_binmin, intra_binwidth, intra_energybin[0]);
                
                // pair correlation function
                if (ARG_pair_correlation==1){
                    for (int pairs=0; pairs<(2*ARG_numberofbeads); pairs++) {
                        intra_obs[0].pair_corr[pairs] += (pair_correlation_obs[pairs]*intra_boltzman_factors[0]);
                    }
                }
                break;
                
            // torsion + square well potential, range of energy values
            case 2:
                intrapot_torsion_well_scan(intra_boltzman_factors, intra_energy, intra_interactions_hist, intra_highest_boltzmanfactor, intra_parameter1, intra_parameter2, normal_obs.pt, number_of_torsions, int_polymer, ARG_numberofbeads);
                //intrapot_torsion_well_test(intra_boltzman_factors, intra_energy, intra_interactions_hist, intra_interactions_testhist,intra_highest_boltzmanfactor, intra_parameter1, intra_parameter2, normal_obs.pt, number_of_torsions, int_polymer, ARG_numberofbeads);
                
                for (int dim1=0; dim1<100; ++dim1) {
                    intra_sum_of_enrgyboltz[dim1][tenth] += (intra_boltzman_factors[dim1]*intra_energy[dim1]);
                    intra_sum_of_boltzfactors[dim1][tenth] += intra_boltzman_factors[dim1];
                    intra_obs[dim1].err_ee2[tenth] += normal_obs.ee2 * intra_boltzman_factors[dim1];
                    intra_obs[dim1].err_rgyr[tenth] += normal_obs.rgyr * intra_boltzman_factors[dim1];
                    intra_obs[dim1].err_pt[tenth] += normal_obs.pt * intra_boltzman_factors[dim1];
                    intra_binenergy(intra_energy[dim1], intra_binmin, intra_binwidth, intra_energybin[dim1]);
                    // pair correlation function
                    if (ARG_pair_correlation==1){
                        for (int pairs=0; pairs<(2*ARG_numberofbeads); pairs++) {
                            intra_obs[dim1].pair_corr[pairs] += (pair_correlation_obs[pairs]*intra_boltzman_factors[dim1]);
                        }
                    }
                }
                
                break;
                
            default:
                break;
                
        }// boltzman factors and energies have been determined
            // end of anything related to intramolecular potentials
        
        
        
        //--------------------------------------------------------------------------------------------
        //--------------------------------------------------------------------------------------------
        // everything has been calculated, now is the opportunity to look at convergence
        
        if ((frame+1)%permill_frames==0) {
            // convergence of variables with equal weights
            convergence(&normal_obs, (double)number_of_torsions, (frame+1), conv_obs_file);
            
            
            switch (ARG_intra_potential) {
                    
                case 1:
                    weighted_convergence(intra_obs, 1,intra_sum_of_boltzfactors, (double)number_of_torsions, (frame+1), conv_intraobs_file);
                    weights_growth(intra_sum_of_boltzfactors, 1, (frame+1), convergence_intraweights);
                    break;
                    
                    // convergence of weighted ensemble
                case 2:
                    weighted_convergence(intra_obs, 100, intra_sum_of_boltzfactors, (double)number_of_torsions, (frame+1), conv_intraobs_file);
                    weights_growth(intra_sum_of_boltzfactors, 100, (frame+1), convergence_intraweights);
                    break;
                    
                default:
                    break;
            }
            
            
            if ((frame+1)%percent_frames==0) {
                
                cent = (frame+1)/percent_frames;
                
                log_out(output_file, "Finished %i%%\t...remaining time: %f seconds \n", (cent), ((time_of_day()-time[1])*(100-cent)/(cent)));
                
                if ((frame+1)%tenth_frames==0) {
                    // every bin after the first will also include the attempts in the previous bin!
                    log_attempts[tenth] = recalc_attempts(attempts_successes, numof_frags_bricks, numberofbricks, ARG_typeofrun);
                }
                
                
            }
        }
        
    } // end of loop over number of frames
    
    
    

    
    
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    /* Post-Process */
    time[2] = time_of_day();
    
    // Maybe there should be a function, which returns means and errors; it calls these subroutines ....
    
    log_attempts[10] = recalc_attempts(attempts_successes, numof_frags_bricks, numberofbricks, ARG_typeofrun);
    for (int dim1=9; dim1>0; --dim1) {
        log_attempts[dim1] = log( exp(log_attempts[dim1]) - exp(log_attempts[dim1-1]) );
    }
    
    // get means and errors
    normal_obs.ee2 = sqrt( average(normal_obs.err_ee2, ARG_numberofframes) );
    normal_obs.err_ee2[10] = error_sq_ten(normal_obs.ee2, normal_obs.err_ee2, ARG_numberofframes);
    normal_obs.rgyr = sqrt( average(normal_obs.err_rgyr, ARG_numberofframes) );
    normal_obs.err_rgyr[10] = error_sq_ten(normal_obs.rgyr, normal_obs.err_rgyr, ARG_numberofframes);
    normal_obs.S = (log((double)ARG_numberofframes) - log_attempts[10]);
    
    if (ARG_torsion==1) {
        normal_obs.pt = average(normal_obs.err_pt, ARG_numberofframes);
        normal_obs.err_pt[10] = error_ten(normal_obs.pt, normal_obs.err_pt, ARG_numberofframes);
    }
    if (ARG_sasa==1) {
        normal_obs.dsasa = average(normal_obs.err_dsasa, ARG_numberofframes);
        normal_obs.err_dsasa[10] = error_ten(normal_obs.dsasa, normal_obs.err_dsasa, ARG_numberofframes);
    }
    
    
    
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    /* Output */
    time[3] = time_of_day();
    
    // print the output to screen and to the output_file
    
    log_out(output_file, "\n\n");
    log_out(output_file, "Final Output\n");
    
    log_out(output_file, "Attempts to compute ensemble: %e  \n", exp(log_attempts[10]));
    log_out(output_file, "\tDelta S / kB (FWW -> SAWn): %f \n",  normal_obs.S);
    
    log_out(output_file, "Flory Radius: %f +- %f \n", normal_obs.ee2, normal_obs.err_ee2[10]);
    log_out(output_file, "Radius of Gyration: %f +- %f \n", normal_obs.rgyr, normal_obs.err_rgyr[10]);
    if (ARG_torsion==1) {
        log_out(output_file, "Probability of trans = %f +- %f\n", (normal_obs.pt/(double)number_of_torsions), (normal_obs.err_pt[10]/(double)number_of_torsions));

    }
    if (ARG_sasa==1) {
        log_out(output_file, "Delta SASA = %f +- %f\n", normal_obs.dsasa, normal_obs.err_dsasa[10]);
    }

    
    
    
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    
    /* Output to file */
    
    // pair correlation function
    if (ARG_pair_correlation==1){
        
        if (ARG_intra_potential==0) {
            fprintf(pair_correl_file,"# r_0-r_1 pair_correlation \n");
            for (int dim1=0; dim1<(2*ARG_numberofbeads); dim1++){
                fprintf(pair_correl_file, "%i %e \n", dim1, (normal_obs.pair_corr[dim1]/(double)ARG_numberofframes));
            }
        }
        else if (ARG_intra_potential==2){
            fprintf(pair_correl_file,"# r_0-r_1 pair_correlation(unweigthed) pair_correlation(weigthed)\n");
            for (int dim1=0; dim1<(2*ARG_numberofbeads); dim1++){
                fprintf(pair_correl_file, "%i %e %e \n", dim1, (normal_obs.pair_corr[dim1]/(double)ARG_numberofframes), (intra_obs[46].pair_corr[dim1]/average(intra_sum_of_boltzfactors[46], 1)));
            }
        }
        
        
    }
    
    
    
    switch (ARG_intra_potential) {
            
        case 1:
            // prints the weighted observables to file with error estimate
            fprintf(intra_pot_file, "### 1-9 well, torsion eps, Rf, Rf_err, Rg, Rg_err, pT, pT_err, DS, DS_err, <exp(-E/kT)>/exp(-Emax/kT) \n" );

            intra_obs[0].ee2 = sqrt(weighted_average(intra_obs[0].err_ee2, intra_sum_of_boltzfactors[0]));
            intra_obs[0].err_ee2[10] = weighted_error_sq_ten(intra_obs[0].ee2, intra_obs[0].err_ee2, intra_sum_of_boltzfactors[0]);
            
            intra_obs[0].rgyr = sqrt(weighted_average(intra_obs[0].err_rgyr, intra_sum_of_boltzfactors[0]));
            intra_obs[0].err_rgyr[10] = weighted_error_sq_ten(intra_obs[0].rgyr, intra_obs[0].err_rgyr, intra_sum_of_boltzfactors[0]);
            
            intra_obs[0].pt = weighted_average(intra_obs[0].err_pt, intra_sum_of_boltzfactors[0]);
            intra_obs[0].err_pt[10] = weighted_error_ten(intra_obs[0].pt, intra_obs[0].err_pt, intra_sum_of_boltzfactors[0]);
            
            intra_obs[0].S = intra_entropy(intra_sum_of_boltzfactors[0], intra_sum_of_enrgyboltz[0], log_attempts[10]);
            intra_obs[0].err_S = intra_entropy_error_ten(intra_obs[0].S, intra_sum_of_boltzfactors[0], intra_sum_of_enrgyboltz[0], log_attempts);
            
            fprintf(intra_pot_file, "%f %f %e %e %e %e %e %e %e %e %e \n", intra_parameter1[0], intra_parameter2[0], intra_obs[0].ee2, intra_obs[0].err_ee2[10], intra_obs[0].rgyr, intra_obs[0].err_rgyr[10], (intra_obs[0].pt/(double)number_of_torsions), (intra_obs[0].err_pt[10]/(double)number_of_torsions), intra_obs[0].S, intra_obs[0].err_S, intra_loss_of_conf(intra_sum_of_boltzfactors[0], intra_highest_boltzmanfactor[0], ARG_numberofframes));

            
            // histogram over 1-9 interacitons
            fprintf(intra_interactions_file, "# number-of-nn-interactions, occurence \n");
            for (int dim1=0; dim1<ARG_numberofbeads; ++dim1) {
                fprintf(intra_interactions_file, "%i %i \n", dim1, intra_interactions_hist[dim1]);
            }
            
            // test histogram
            //      fprintf(intra_interactions_testfile, "# sperating-bonds 0_contacts 1_contacts 2_contacts\n");
            //      for (int dim1=0; dim1<98; dim1+=3) {
            //          fprintf(intra_interactions_testfile, "%i %e %e %e \n", (dim1/3+4), intra_interactions_testhist[dim1]/average(intra_sum_of_boltzfactors[46], 1), intra_interactions_testhist[dim1+1]/average(intra_sum_of_boltzfactors[46], 1), intra_interactions_testhist[dim1+2]/average(intra_sum_of_boltzfactors[46], 1));
            //     }
            
            // histogram of intra energies
            fprintf(intra_boltzman_factors_hist, "# energy, bincount-for-these-parameters");
            for (int dim1=0; dim1<((intra_binmax-intra_binmin)/intra_binwidth); ++dim1) {
                fprintf(intra_boltzman_factors_hist, "%f %i \n", (intra_binmin+(double)dim1*intra_binwidth), intra_energybin[0][dim1]);
            }
            break;
            
        case 2:
            // prints the weighted observables to file with error estimate
            fprintf(intra_pot_file, "### 1-9 well, torsion eps, Rf, Rf_err, Rg, Rg_err, pT, pT_err, DS, DS_err, <exp(-E/kT)>/exp(-Emax/kT) \n" );
            for (int dim1=0; dim1<10; dim1++) {
                for (int dim2=0; dim2<10; dim2++) {
                    
                    counter = dim1*10 + dim2;
                    
                    intra_obs[counter].ee2 = sqrt(weighted_average(intra_obs[counter].err_ee2, intra_sum_of_boltzfactors[counter]));
                    intra_obs[counter].err_ee2[10] = weighted_error_sq_ten(intra_obs[counter].ee2, intra_obs[counter].err_ee2, intra_sum_of_boltzfactors[counter]);
                    
                    intra_obs[counter].rgyr = sqrt(weighted_average(intra_obs[counter].err_rgyr, intra_sum_of_boltzfactors[counter]));
                    intra_obs[counter].err_rgyr[10] = weighted_error_sq_ten(intra_obs[counter].rgyr, intra_obs[counter].err_rgyr, intra_sum_of_boltzfactors[counter]);
                    
                    intra_obs[counter].pt = weighted_average(intra_obs[counter].err_pt, intra_sum_of_boltzfactors[counter]);
                    intra_obs[counter].err_pt[10] = weighted_error_ten(intra_obs[counter].pt, intra_obs[counter].err_pt, intra_sum_of_boltzfactors[counter]);
                    
                    intra_obs[counter].S = intra_entropy(intra_sum_of_boltzfactors[counter], intra_sum_of_enrgyboltz[counter], log_attempts[10]);
                    intra_obs[counter].err_S = intra_entropy_error_ten(intra_obs[counter].S, intra_sum_of_boltzfactors[counter], intra_sum_of_enrgyboltz[counter], log_attempts);
                    
                    fprintf(intra_pot_file, "%f %f %e %e %e %e %e %e %e %e %e \n", intra_parameter1[dim1], intra_parameter2[dim2], intra_obs[counter].ee2, intra_obs[counter].err_ee2[10], intra_obs[counter].rgyr, intra_obs[counter].err_rgyr[10], (intra_obs[counter].pt/(double)number_of_torsions), (intra_obs[counter].err_pt[10]/(double)number_of_torsions), intra_obs[counter].S, intra_obs[counter].err_S, intra_loss_of_conf(intra_sum_of_boltzfactors[counter], intra_highest_boltzmanfactor[counter], ARG_numberofframes));
                }
            }
            
            // histogram over 1-9 interacitons
            fprintf(intra_interactions_file, "# number-of-nn-interactions, occurence \n");
            for (int dim1=0; dim1<ARG_numberofbeads; ++dim1) {
                fprintf(intra_interactions_file, "%i %i \n", dim1, intra_interactions_hist[dim1]);
            }
            
            // test histogram
            //      fprintf(intra_interactions_testfile, "# sperating-bonds 0_contacts 1_contacts 2_contacts\n");
            //      for (int dim1=0; dim1<98; dim1+=3) {
            //          fprintf(intra_interactions_testfile, "%i %e %e %e \n", (dim1/3+4), intra_interactions_testhist[dim1]/average(intra_sum_of_boltzfactors[46], 1), intra_interactions_testhist[dim1+1]/average(intra_sum_of_boltzfactors[46], 1), intra_interactions_testhist[dim1+2]/average(intra_sum_of_boltzfactors[46], 1));
            //     }
            
            // histogram of intra energies
            fprintf(intra_boltzman_factors_hist, "# energy, bincount-for-these-parameters");
            for (int dim1=0; dim1<((intra_binmax-intra_binmin)/intra_binwidth); ++dim1) {
                fprintf(intra_boltzman_factors_hist, "%f ", (intra_binmin+(double)dim1*intra_binwidth));
                for (int dim2=0; dim2<100; ++dim2) {
                    fprintf(intra_boltzman_factors_hist, "%i ", intra_energybin[dim2][dim1]);
                }
                fprintf(intra_boltzman_factors_hist, "\n");
            }

            break;
            
        default:
            break;
    }


    
    
    
    
    log_out(output_file, "\nComputation of Chains: %f seconds \nTotal Runtime of program: %fseconds \n", (time[2]-time[1]), (time[3]-time[0]));
    log_out(output_file, "End of Programm!\n------------------------------------------------\n");
    
    // make sure everything gets printed
    fflush(stdout);
    
    
    
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------------------------------------------------------------------
    // close every remaining file and free remaining pointers!
    // cleaning up ;-)
    
    //close all files
    if (ARG_typeofrun==0){// dcd-file
        close_file_read(v);
    }
    
    fclose(conv_obs_file);
    fclose(pair_correl_file);
    fclose(intra_pot_file);
    fclose(intra_boltzman_factors_hist);
    fclose(convergence_intraweights);
    fclose(intra_interactions_file);
    fclose(conv_intraobs_file);
    
    // very last file to close
    fclose(output_file);
    
    
    
    //-----------------------------
    // free pointers
    
    
    // free general variables
    free(normal_obs.err_ee2);
    free(normal_obs.err_rgyr);
    
    // torsion angles
    if (ARG_torsion==1) {
        free(normal_obs.err_pt);
        if (ARG_intra_potential>0) {
            for (int dim1=0; dim1<100; dim1++) {
                free(intra_obs[dim1].err_pt);
            }
        }
    }
    
    // D-SASA
    if (ARG_sasa==1) {
        free(normal_obs.err_dsasa);
        if (ARG_intra_potential>0) {
            for (int dim1=0; dim1<100; dim1++) {
                free(intra_obs[dim1].err_dsasa);
            }
        }
    }
    
    // pair correlation function
    if (ARG_pair_correlation==1) {
        free(pair_correlation_obs);
        free(normal_obs.pair_corr);
        if (ARG_intra_potential>0) {
            for (int dim1=0; dim1<100; dim1++) {
                free(intra_obs[dim1].pair_corr);
            }
        }
    }
    
    free(attempts_successes);
    
    // free arrays which held polymer coordinates
    for (int dim1=0; dim1<ARG_numberofbeads;dim1++){
        free(int_polymer[dim1]);
        free(double_polymer[dim1]);
    }
    free(int_polymer);
    free(double_polymer);
    
    
    
    // free building blocks
    if (ARG_typeofrun==13 || ARG_typeofrun==14 || ARG_typeofrun==23 || ARG_typeofrun==24 || ARG_typeofrun==33 || ARG_typeofrun==34){
        for (int dim1=0; dim1<numberofbricks; ++dim1) {
            for (int dim2=0; dim2<ARG_blength; ++dim2) {
                free(building_bricks[dim1][dim2]);
            }
            free(building_bricks[dim1]);
        }
        free(building_bricks);
    }
    
    
    
    // free potential variables
    if (ARG_intra_potential > 0){
        switch (ARG_intra_potential) {
            case 1:
                free(intra_sum_of_boltzfactors[0]);
                break;
                
            case 2:
                for (int dim1=0; dim1<100; dim1++){
                    free(intra_sum_of_boltzfactors[dim1]);
                }
                break;

            default:
                break;
        }

        free(intra_energy);
        free(intra_boltzman_factors);
        free(intra_sum_of_boltzfactors);
        free(intra_sum_of_enrgyboltz);
        free(intra_interactions_hist);
  //      free(intra_interactions_testhist);
        free(intra_parameter1);
        free(intra_parameter2);
    }
    
    
    



    
    // end of programm
    return EXIT_SUCCESS;
}

    
    
    
    
    
    
    
    
    
    
    

