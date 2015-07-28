//
//  read_dcd.h
//  
//
//  Created by Johannes Dietschreit on 16.06.15.
//
//

#ifndef ____read_dcd__
#define ____read_dcd__

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661922
#endif

#define RECSCALE32BIT 1
#define RECSCALE64BIT 2
#define RECSCALEMAX   2
#define     LINE_LENGTH         256
#define     STRING_LENGTH       10


#include <string.h>
#include "fastio.h"
#include "molfile_plugin.h"

typedef struct {
    fio_fd fd;
    int natoms;
    int nsets;
    int setsread;
    int istart;
    int nsavc;
    double delta;
    int nfixed;
    float *x, *y, *z;
    int *freeind;
    float *fixedcoords;
    int reverse;
    int charmm;
    int first;
    int with_unitcell;
} dcdhandle;

typedef char LINE[LINE_LENGTH];
typedef char STRING[STRING_LENGTH];

typedef struct{
    int     NumberOfAtom;	// Total number of atoms
    STRING 	segid;		// Segment ID
    STRING	atomtype;	// Charmm atom type
    double 	mass;		// Mass in a.m.u.
    double  sqrtmass;       // Square root of the mass in a.m.u.^(1/2) for covariance matrix
}  CharmmPSF;


// pre-declarations

void subString(const char *source, int begin, int end, char *destination);

static void print_dcderror(const char *func, int errcode);

static int read_dcdheader(fio_fd fd, int *N, int *NSET, int *ISTART,
                          int *NSAVC, double *DELTA, int *NAMNF,
                          int **FREEINDEXES, float **fixedcoords, int *reverseEndian,
                          int *charmm);

static int read_charmm_extrablock(fio_fd fd, int charmm, int reverseEndian,
                                  float *unitcell);

static int read_fixed_atoms(fio_fd fd, int N, int num_free, const int *indexes,
                            int reverseEndian, const float *fixedcoords,
                            float *freeatoms, float *pos, int charmm);

static int read_charmm_4dim(fio_fd fd, int charmm, int reverseEndian);

static int read_dcdstep(fio_fd fd, int N, float *X, float *Y, float *Z,
                        float *unitcell, int num_fixed,
                        int first, int *indexes, float *fixedcoords,
                        int reverseEndian, int charmm);

static int skip_dcdstep(fio_fd fd, int natoms, int nfixed, int charmm);

// static int write_dcdstep(fio_fd fd, int curframe, int curstep, int N, const float *X, const float *Y, const float *Z, const double *unitcell, int charmm);

// static int write_dcdheader(fio_fd fd, const char *remarks, int N, int ISTART, int NSAVC, double DELTA, int with_unitcell, int charmm);

static void close_dcd_read(int *indexes, float *fixedcoords);

void *open_dcd_read(const char *path, const char *filetype,
                           int *natoms);

static int read_next_timestep(void *v, int natoms, molfile_timestep_t *ts);

void close_file_read(void *v);

// static void *open_dcd_write(const char *path, const char *filetype, int natoms);

// static int write_timestep(void *v, const molfile_timestep_t *ts);

// static void close_file_write(void *v);

void dcd_to_polymer(double (**polymer), void *v, int natoms, molfile_timestep_t *ts, dcdhandle *dcd_hadle, unsigned long i,int numberofbeads);


#endif /* defined(____read_dcd__) */
