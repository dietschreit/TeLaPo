//
//  statistics.h
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#ifndef ____statistics__
#define ____statistics__

double average(double (*obs_ptr), unsigned long numberofframes);

double weighted_average(double (*obs_ptr), double (*weights_ptr));

double error_ten(double mean, double (*obs_ptr), unsigned long numberofframes);

double error_sq_ten(double mean, double (*obs_ptr), unsigned long numberofframes);

double weighted_error_ten(double mean, double (*means_ptr), double (*weights_ptr));

double weighted_error_sq_ten(double mean, double (*means_ptr), double (*weights_ptr));

int histogram(int numbins, double (*obs_ptr), unsigned long numberofframes, FILE *file_ptr);

int histogram_boltz(int numbins, double (**obs_ptr), unsigned long numberofframes, int nfactors, FILE *file_ptr);

void convergence(OBSERVABLE *obs_ptr, double ntor, unsigned long convpoint, FILE *file_ptr);

void weighted_convergence(OBSERVABLE *obs_ptr, int num_eps, double (**weights_ptr), double ntor, unsigned long convpoint, FILE *file_ptr);

void weights_growth(double (**obs_ptr), int num_eps, unsigned long convpoint, FILE *file_ptr);

#endif /* defined(____statistics__) */
