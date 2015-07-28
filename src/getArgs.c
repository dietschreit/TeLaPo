//
//  getArgs.c
//  
//
//  Created by Johannes Dietschreit on 29.05.15.
//
//

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "output.h"
#include "time.h"
#include "getArgs.h"



void usage_error(void){
    // this needs to be filled with all the options one can choose from
    
    printf("--------------------------------------------------------------------------\n");
    printf("Usage: \n\n");
    printf("Runs \n");
    printf("\t-tetra_rw      [NumberOfAtoms] : simple random walk on tetrahedral lattice \n");
    printf("\t-tetra_fww     [NumberOfAtoms] : forward random walk on tetrahedral lattice \n");
    printf("\t-tetra_saw0    [NumberOfAtoms] : forward random walk on tetrahedral lattice \n");
    printf("\t-tetra_saw1    [NumberOfAtoms] : self-avoiding walk on tetrahedral lattice \n");
    printf("\t-tetra_saw2    [NumberOfAtoms] : self-avoiding walk 2 on tetrahedral lattice \n");
    printf("\t-tetra_saw3    [NumberOfAtoms] : self-avoiding walk 3 on tetrahedral lattice \n");
    printf("\t-tetra_fsaw3   [NumberOfAtoms] : tetra_saw3 built from fragments (needs -flength or -autocalc) \n");
    printf("\t-tetra_bsaw3   [NumberOfAtoms] : tetra_saw3 built from bricks (needs -blength or -autocalc) \n");
    printf("\t-tetra_f+bsaw3 [NumberOfAtoms] : tetra_saw3 built from fragments and bricks (needs -flength [] and -blength []) \n");
    printf("\n");
    printf("mandatory parameters \n");
    printf("\t-num [NumberOfConformers] : number of successfully generated chains \n");
    printf("\t-flength [FragmentLength] : number of bonds in a fragment \n");
    printf("\t-blength    [BrickLength] : number of bonds in a brick (check RAM!) \n");
    printf("\t-randseed    [Randomseed] : seed for pseudo random number generator (positive integer) \n");
    printf("\t-bond        [BondLength] : bond length in Angstrom (grid constant) \n");
    printf("\n");
    printf("further options \n");
    printf("\t-torsion        : switches torsional analysis on \n");
    printf("\t-pair_correl    : switches pair correlation function around central atom on \n");
    printf("\t-intra_tor_well : intramolecular potential with torsional and nearest neighbor part (needs -well_eps [] and -tor_eps []) \n");
    printf("\t-intra_tor_well_scan : same potential, scans over energy range in one run (needs -well_eps [][] and -tor_eps [][]) \n");
    printf("\t-well_eps [Energy] [StepWidth]: Energy parameter in kT, the [StepWidth] is not always necessary \n");
    printf("\t-tor_eps  [Energy] [StepWidth]: Energy parameter in kT, the [StepWidth] is not always necessary \n");
    printf("--------------------------------------------------------------------------\n");




 

    exit(1);
    
}


//----------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------------------------


void getArgs(FILE *output_ptr, int argc, char **argv)
{
    int i=0;
    
    if (argc < 3){
        log_out(output_ptr, "ERROR: The number of arguments is insufficient for any kind of run! \n\n");
        usage_error();
    }
    
    //-----------------------------------------------------------------
    
    // Get DCD trajectory file name
    for( i = 1; i < argc; i++) {
        if (strcmp("-dcd",argv[i])==0) {
            ARG_typeofrun=0;
            if ( (i+1)<=(argc-1)) {
                dcdFileName=argv[i+1];
            }
            else {
                log_out(output_ptr, "ERROR: Cannot read the specified DCD-file! \n\n");
                usage_error(); }
            break;

        }
    }
    
    
    //------------------------------------------------------------------
    // define the walk and number of atoms
    
    
    // RW on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_rw",argv[i])==0) {
            ARG_typeofrun=1;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // FWW on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_fww",argv[i])==0 || strcmp("-tetra_saw0",argv[i])==0) {
            ARG_typeofrun=2;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    //-----------------------------------------------------------------
    // all SAW1s
    
    // SAW1 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_saw1",argv[i])==0) {
            ARG_typeofrun=10;
            ARG_strictness=1;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW1 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_fsaw1",argv[i])==0) {
            ARG_typeofrun=11;
            ARG_strictness=1;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW1 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_ffsaw1",argv[i])==0) {
            ARG_typeofrun=12;
            ARG_strictness=1;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW1 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_bsaw1",argv[i])==0) {
            ARG_typeofrun=13;
            ARG_strictness=1;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW1 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_f+bsaw1",argv[i])==0) {
            ARG_typeofrun=14;
            ARG_strictness=1;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    //-----------------------------------------------------------------
    // all SAW2s
    
    // SAW2 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_saw2",argv[i])==0) {
            ARG_typeofrun=20;
            ARG_strictness=2;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW2 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_fsaw2",argv[i])==0) {
            ARG_typeofrun=21;
            ARG_strictness=2;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW2 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_ffsaw2",argv[i])==0) {
            ARG_typeofrun=22;
            ARG_strictness=2;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW2 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_bsaw2",argv[i])==0) {
            ARG_typeofrun=23;
            ARG_strictness=2;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW2 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_f+bsaw2",argv[i])==0) {
            ARG_typeofrun=24;
            ARG_strictness=2;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    //-----------------------------------------------------------------
    // all SAW3s
    
    // SAW3 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_saw3",argv[i])==0) {
            ARG_typeofrun=30;
            ARG_strictness=3;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW3 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_fsaw3",argv[i])==0) {
            ARG_typeofrun=31;
            ARG_strictness=3;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW3 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_ffsaw3",argv[i])==0) {
            ARG_typeofrun=32;
            ARG_strictness=3;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW3 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_bsaw3",argv[i])==0) {
            ARG_typeofrun=33;
            ARG_strictness=3;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // SAW3 on tetrahedral lattice
    for( i = 1; i < argc; i++) {
        if (strcmp("-tetra_f+bsaw3",argv[i])==0) {
            ARG_typeofrun=34;
            ARG_strictness=3;
            if ( (i+1)<=(argc-1)) {
                ARG_numberofbeads=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your chosen generation method or number of atoms! \n\n");
                usage_error(); }
            break;
        }
    }
    
    //-----------------------------------------------------------------
    // sets number of frames
    for( i = 1; i < argc; i++) {
        if (strcmp("-num",argv[i])==0) {
            if ( (i+1)<=(argc-1)) {
                ARG_numberofframes=atoll(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with number of frames! \n\n");
                usage_error(); }
            break;
        }
    }
    
    for( i = 1; i < argc; i++) {
        if (strcmp("-blength",argv[i])==0) {
            if ( (i+1)<=(argc-1)) {
                ARG_blength=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with blength! \n\n");
                usage_error(); }
            break;
        }
    }
    
    for( i = 1; i < argc; i++) {
        if (strcmp("-flength",argv[i])==0) {
            if ( (i+1)<=(argc-1)) {
                ARG_flength=atoi(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with flength! \n\n");
                usage_error(); }
            break;
        }
    }
    
    
    // sets number of frames
    for( i = 1; i < argc; i++) {
        if (strcmp("-randseed",argv[i])==0) {
            if ( (i+1)<=(argc-1)) {
                ARG_randomseed=atol(argv[i+1]);
                ARG_randomseed*=(-1);
            }
            else {
                log_out(output_ptr, "ERROR: Problem with your randomseed! \n\n");
                usage_error(); }
            break;
        }
    }
    if (ARG_randomseed==0){
        ARG_randomseed = - (long) time_of_day();
    }
    
    
    
    //-----------------------------------------------------------------
    // additional options, never mandatory
    
    // bond length
    for( i = 1; i < argc; i++) {
        if (strcmp("-bond",argv[i])==0) {
            if ( (i+1)<=(argc-1)) {
                ARG_bondlength = atof(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: There is a mistake in your choice of bond length! \n\n");
                usage_error(); }
            break;
        }
    }
    
    // torsional analysis: on
    for( i = 1; i < argc; i++) {
        if (strcmp("-torsion",argv[i])==0) {
            ARG_torsion = 1;
            break;
        }
    }
    
    // torsional analysis: on
    for( i = 1; i < argc; i++) {
        if (strcmp("-pair_correl",argv[i])==0) {
            ARG_pair_correlation = 1;
            break;
        }
    }
    
    
    //-----------------------------------------------------------------
    //intramolecular forces
    
    //--------------------------------------------
    // torsion potential + 1-9 well
    for( i = 1; i < argc; i++) {
        if (strcmp("-intra_tor_well",argv[i])==0) {
            ARG_intra_potential = 1;
            ARG_torsion = 1;
            break;
        }
    }
    
    for( i = 1; i < argc; i++) {
        if (strcmp("-intra_tor_well_scan",argv[i])==0) {
            ARG_intra_potential = 2;
            ARG_torsion = 1;
            break;
        }
    }
    
    for( i = 1; i < argc; i++) {
        if (strcmp("-well_eps",argv[i])==0) {
            if ( (i+1)<=(argc-1)) {
                ARG_intra_parameter1[0] = atof(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: There is a mistake in your choice of -well_eps! \n\n");
                usage_error(); }
            if (ARG_intra_potential==2){
                if ( (i+2)<=(argc-1)) {
                    ARG_intra_parameter1[1] = atof(argv[i+2]);
                }
                else {
                    log_out(output_ptr, "ERROR: There is a mistake in your choice of -well_eps! \n\n");
                    usage_error(); }
            }
            break;
        }
    }
    for( i = 1; i < argc; i++) {
        if (strcmp("-tor_eps",argv[i])==0) {
            if ( (i+1)<=(argc-1)) {
                ARG_intra_parameter2[0] = atof(argv[i+1]);
            }
            else {
                log_out(output_ptr, "ERROR: There is a mistake in your choice of -tor_eps! \n\n");
                usage_error(); }
            if (ARG_intra_potential==2){
                if ( (i+2)<=(argc-1)) {
                    ARG_intra_parameter2[1] = atof(argv[i+2]);
                }
                else {
                    log_out(output_ptr, "ERROR: There is a mistake in your choice of -tor_eps! \n\n");
                    usage_error(); }
            }
            break;
        }
    }
    //--------------------------------------------
    
    
    
    //------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------
    //intramolecular forces
    
    // torsion potential + 1-9 well
    for( i = 1; i < argc; i++) {
        if (strcmp("-info",argv[i])==0) {
            log_out(output_ptr, "\nFull information on program history:\n\n");
            log_out(output_ptr, "\tVersion 1.0 (chain creation and evaluation) published in:\n");
            log_out(output_ptr, "\t\tDietschreit, J. C. B.; Diestler, D. J.; Knapp, E. W.,\n");
            log_out(output_ptr, "\t\tModels for Self-Avoiding Polymer Chains on the Tetrahedral Lattice.\n");
            log_out(output_ptr, "\t\tMacromol. Theory Simul., 2014, 23, 452-463\n");
            log_out(output_ptr, "\tVersion 1.1 (intramolecular potential) published in:\n");
            log_out(output_ptr, "\t\tDietschreit, J. C. B.; Diestler, D. J.; Knapp, E. W.,\n");
            log_out(output_ptr, "\t\tTitle.\n");
            log_out(output_ptr, "\t\tJournal. 2015, ???\n");
            log_out(output_ptr, "\tVersion 1.2 (interactions with a surface) uses algorithms from:\n");
            log_out(output_ptr, "\t\tNumata, J.; Juneja, A.; Diestler, D. J.; Knapp, E. W.,\n");
            log_out(output_ptr, "\t\tInfluence of Spacer-Receptor Interactions on the Stability of Bivalent Ligand-Receptor Complexes.\n");
            log_out(output_ptr, "\t\tJ. Phys. Chem. B, 2012, 116 (8), 2595-2604\n");
            exit(1);
        }
    }

    
    
    
}


//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------------------------------
//

void print_set_options(FILE *output_ptr, int typeofrun, int flength, int fflength, int blength, int numberofbeads, unsigned long long numberofframes, long seed, double bond, int torsion, int intramolecular, double intra_eps1[2], double intra_eps2[2]){
    
    
    
    
    //-----------------------------------------------------------------
    // sanity checks
    //-----------------------------------------------------------------
    
    if (ARG_typeofrun<0){
        log_out(output_ptr, "ERROR: The chosen generation method does not exist! \n\n");
        usage_error();
    }
    if (ARG_numberofbeads<3){
        log_out(output_ptr, "ERROR: The given number of beads makes no sense! \n\n");
        usage_error();
    }
    if (ARG_numberofframes<100){
        log_out(output_ptr, "ERROR: Please set the number of frames at least to 100! The ensemble is otherwise statistically not relevant! \n\n");
        usage_error();
    }
    if ((ARG_numberofframes%100)!=0){
        log_out(output_ptr, "ERROR: Please set the number of frames to a multiple of 100, this makes error estimation easier! \n\n");
        usage_error();
    }
    if ((ARG_typeofrun==11 || ARG_typeofrun==12 || ARG_typeofrun==14 || ARG_typeofrun==21 || ARG_typeofrun==22 || ARG_typeofrun==24 || ARG_typeofrun==31 || ARG_typeofrun==32 || ARG_typeofrun==34) && (((ARG_flength%2)!=0) || ARG_flength>(ARG_numberofbeads-1))){
        log_out(output_ptr, "ERROR: If you want to use fragments, set it to an even number! It should be smaller than the number of atoms! \n\n");
        usage_error();
    }
    if ((ARG_typeofrun==13 || ARG_typeofrun==14 || ARG_typeofrun==23 || ARG_typeofrun==24 || ARG_typeofrun==33 || ARG_typeofrun==34) && ((ARG_blength%2)!=0 || ARG_blength>16)){
        log_out(output_ptr, "ERROR: If you want to use bricks, set it to an even number! \n\n");
        usage_error();
    }
    if (ARG_randomseed>0){
        log_out(output_ptr, "ERROR: The chosen random seed is faulty! \n");
        usage_error();
    }
    

    
    
    
    
    
    /*
     this prints all the options that were set;
     this will give a good overview at the beginning whether the programm does what the user wants
     */

    log_out(output_ptr, "The chosen options for this computation are: \n\n");
    log_out(output_ptr, "mandatory\n");
    
    switch (typeofrun) {
        case 1:
            log_out(output_ptr, "\tRun: Random Walk on Tetrahedral Lattice (-tetra_rw) \n");
            break;
        case 2:
            log_out(output_ptr, "\tRun: Foreward Random Walk on Tetrahedral Lattice (-tetra_fww) \n");
            break;
        case 10:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 1 on Tetrahedral Lattice (-tetra_saw1) \n");
            break;
        case 11:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 1 on Tetrahedral Lattice, fragmented (-tetra_fsaw1) \n");
            log_out(output_ptr, "\t\tFragment Length: %i \n", flength);
            break;
        case 12:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 1 on Tetrahedral Lattice, doubly fragmented (-tetra_ffsaw1) \n");
            log_out(output_ptr, "\t\tFragment Length (long):  %i \n", flength);
            log_out(output_ptr, "\t\tFragment Length (short): %i \n", fflength);
            break;
        case 13:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 1 on Tetrahedral Lattice, bricks (-tetra_bsaw1) \n");
            log_out(output_ptr, "\t\tBrick Length: %i \n", blength);
            break;
        case 14:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 1 on Tetrahedral Lattice, fragments and bricks (-tetra_f+bsaw1) \n");
            log_out(output_ptr, "\t\tFragment Length: %i \n", flength);
            log_out(output_ptr, "\t\tBrick Length: %i \n", blength);
            break;
        case 20:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 2 on Tetrahedral Lattice (-tetra_saw2) \n");
            break;
        case 21:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 2 on Tetrahedral Lattice, fragmented (-tetra_fsaw2) \n");
            log_out(output_ptr, "\t\tFragment Length: %i \n", flength);
            break;
        case 22:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 2 on Tetrahedral Lattice, doubly fragmented (-tetra_ffsaw2) \n");
            log_out(output_ptr, "\t\tFragment Length (long):  %i \n", flength);
            log_out(output_ptr, "\t\tFragment Length (short): %i \n", fflength);
            break;
        case 23:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 2 on Tetrahedral Lattice, bricks (-tetra_bsaw2) \n");
            log_out(output_ptr, "\t\tBrick Length: %i \n", blength);
            break;
        case 24:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 2 on Tetrahedral Lattice, fragments and bricks (-tetra_f+bsaw2) \n");
            log_out(output_ptr, "\t\tFragment Length: %i \n", flength);
            log_out(output_ptr, "\t\tBrick Length: %i \n", blength);
            break;
        case 30:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 3 on Tetrahedral Lattice (-tetra_saw3) \n");
            break;
        case 31:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 3 on Tetrahedral Lattice, fragmented (-tetra_fsaw3) \n");
            log_out(output_ptr, "\t\tFragment Length: %i \n", flength);
            break;
        case 32:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 3 on Tetrahedral Lattice, doubly fragmented (-tetra_ffsaw3) \n");
            log_out(output_ptr, "\t\tFragment Length (long):  %i \n", flength);
            log_out(output_ptr, "\t\tFragment Length (short): %i \n", fflength);
            break;
        case 33:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 3 on Tetrahedral Lattice, bricks (-tetra_bsaw3) \n");
            log_out(output_ptr, "\t\tBrick Length: %i \n", blength);
            break;
        case 34:
            log_out(output_ptr, "\tRun: Self-avoiding Walk 3 on Tetrahedral Lattice, fragments and bricks (-tetra_f+bsaw3) \n");
            log_out(output_ptr, "\t\tFragment Length: %i \n", flength);
            log_out(output_ptr, "\t\tBrick Length: %i \n", blength);
            break;
            
        default:
            break;
    }
    
    log_out(output_ptr, "\tChain Length: %i \n", numberofbeads);
    log_out(output_ptr, "\tEnsemble Size: %llu \n", numberofframes);
    log_out(output_ptr, "\tRandom Seed: %li \n", -seed);
    
    //---------------------------------------------------------------------
    log_out(output_ptr, "\noptional\n");
    
    log_out(output_ptr, "\tBond length: %f \n", bond);
    
    if (torsion==1) {
        log_out(output_ptr, "\tAnalysis of Dihedrals: on \n");
    }
    else {
        log_out(output_ptr, "\tAnalysis of Dihedrals: off \n");
    }
    
    switch (intramolecular) {
        case 1:
            log_out(output_ptr, "\tIntramolecular Potential: nearest neighbor + torsion \n");
            log_out(output_ptr, "\t\tattractive well  : %f kT\n", intra_eps1[0]);
            log_out(output_ptr, "\t\ttorsion parameter: %f kT\n", intra_eps2[0]);
            break;
            
        case 2:
            log_out(output_ptr, "\tIntramolecular Potential (scan): nearest neighbor + torsion \n");
            log_out(output_ptr, "\t\tattractive well  : %f (stepwidth = %f) kT\n", intra_eps1[0], intra_eps1[1]);
            log_out(output_ptr, "\t\ttorsion parameter: %f (stepwidth = %f) kT\n", intra_eps2[0], intra_eps2[1]);
            break;
            
        default:
            log_out(output_ptr, "\tIntramolecular Potential: off \n");
            break;
    }
    
    log_out(output_ptr, "\n-----------------------------");
    log_out(output_ptr, "\nStart of the computation!\n\n");
    
    
    
    
}
















