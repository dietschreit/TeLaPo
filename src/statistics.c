//
//  statistics.c
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "observables.h"
#include "statistics.h"


// averages

double average(double (*obs_ptr), unsigned long numberofframes){
    
    double mean = 0.0;
    
    for (int dim1=0; dim1<10; dim1++){
        
        mean += obs_ptr[dim1];
        
    }
    
    mean /= ((double)numberofframes);
    
    return mean;
    
}


double weighted_average(double (*obs_ptr), double (*weights_ptr)){
    
    double mean = 0.0;
    double weight_sum = 0.0;
    
    for (int dim1=0; dim1<10; dim1++){
        
        mean += obs_ptr[dim1];
        weight_sum += weights_ptr[dim1];
        
    }
    
    mean /= weight_sum;
    
    return mean;
    
}



//=========================================================================================================
// error estimates

double error_ten(double mean, double (*obs_ptr), unsigned long numberofframes){
    
    double *means;
    means = calloc(10, sizeof(double));
    
    double error = 0.0;
    
    unsigned long decade = numberofframes/10;
    
    for (int dim1=0; dim1<10; ++dim1){
        
        // get mean of that part
        means[dim1]= obs_ptr[dim1]/(double)decade;
        
        // start calculating the variance
        error += ( (means[dim1]-mean)*(means[dim1]-mean) );
        
    }
    
    // this is now the variance
    error /= 10.0;
    // standard deviation
    error = sqrt(error);
    
    // free the pointers
    free(means);
    
    // return standard deviation
    return error;
    
}

// this slight different version is needed, since Rf and Rgyr are properties derived as root mean square values
double error_sq_ten(double mean, double (*obs_ptr), unsigned long numberofframes){
    
    double *means;
    means = calloc(10, sizeof(double));
    
    double error = 0.0;
    
    unsigned long decade = numberofframes/10;
    
    for (int dim1=0; dim1<10; ++dim1){
        
        // get mean of that part
        means[dim1]= obs_ptr[dim1] / (double)decade;
        
        // start calculating the variance
        error += ( (sqrt(means[dim1])-mean)*(sqrt(means[dim1])-mean) );
        
    }
    
    // this is now the variance
    error /= 10.0;
    // standard deviation
    error = sqrt(error);
    
    // free the pointers
    free(means);
    
    // return standard deviation
    return error;
    
}

//---------------------------------------------------------------------------------------------------------
// error for variables with boltzman factors

double weighted_error_ten(double mean, double (*means_ptr), double (*weights_ptr)){
    
    double *means;
    means = calloc(10, sizeof(double));
    
    double error = 0.0;
    
    for (int dim1=0; dim1<10; dim1++){
        
        // get mean of that part
        means[dim1]= (means_ptr[dim1] / weights_ptr[dim1]);
        
        // start calculating the variance
        error += ( (means[dim1]-mean)*(means[dim1]-mean) );
        
    }
    
    // this is now the variance
    error /= 10.0;
    // standard deviation
    error = sqrt(error);
    
    // free pointers
    free(means);
    
    // return standard deviation
    return error;
    
}


// this slight different version is needed, since Rf and Rgyr are properties derived as root mean square values
double weighted_error_sq_ten(double mean, double (*means_ptr), double (*weights_ptr)){
    
    double *means;
    means = calloc(10 ,sizeof(double));
    
    double error = 0.0;
    
    for (int dim1=0; dim1<10; dim1++){
        
        
        // get mean of that part
        means[dim1]= (means_ptr[dim1] / weights_ptr[dim1]);
        
        // start calculating the variance
        error += (sqrt(means[dim1])-mean)*(sqrt(means[dim1])-mean);
        
    }
    
    // this is now the variance
    error /= 10.0;
    // standard deviation
    error = sqrt(error);
    
    // free pointers
    free(means);
    
    // return standard deviation
    return error;
    
}


//=========================================================================================================
// binning

int histogram(int numbins, double (*obs_ptr), unsigned long numberofframes, FILE *file_ptr){
    
    // to give them any value
    double min = obs_ptr[0];
    double max = obs_ptr[0];
    
    // determining min and max
    for (int dim1=1; dim1<numberofframes; dim1++){
        if (min>obs_ptr[dim1]){
            min = obs_ptr[dim1];
        }
        if (max<obs_ptr[dim1]){
            max = obs_ptr[dim1];
        }
    }
    
    int *hist;
    hist = calloc(numbins, sizeof(int));
    
    int bin;
    double width = (max-min)/(double)numbins;
    
    //    printf("min:  %f; max: %f; width: %f; numbins: %i\n", min, max, width, numbins);
    
    // doing the actual binning
    for (int dim1=0; dim1<numberofframes; dim1++){
        bin = floor((obs_ptr[dim1]-min)/width);
        if (bin == numbins){
            bin-=1;
        }
        else if ((bin > numbins) || (bin < 0)){
            fprintf(stderr, "The bin is out of range, bin=%i \n", bin);
            return false;
        }
        
        // raise the counter in this bin
        hist[bin]++;
        
    }
    
    
    // printing the histogram to file
    fprintf(file_ptr, "### bin -- counter -- normalized \n");
    for (int dim1=0; dim1<numbins; dim1++){
        fprintf(file_ptr, "%f %i %f \n", (min + (0.5+(double)dim1)*width), hist[dim1], ((double)hist[dim1]/((double)numbins)*width) );
    }
    
    free(hist);
    
    
    return true;
    
}


int histogram_boltz(int numbins, double (**obs_ptr), unsigned long numberofframes, int nfactors, FILE *file_ptr){
    
    double *min = malloc(nfactors * sizeof(double));
    double *max = malloc(nfactors * sizeof(double));
    
    // to give them any value
    for (int dim1=0; dim1<nfactors; ++dim1) {
        min[dim1] = obs_ptr[dim1][0];
        max[dim1] = obs_ptr[dim1][0];
        
        // determining min and max
        for (int dim2=1; dim2<numberofframes; ++dim2){
            if (min[dim1]>obs_ptr[dim1][dim2]){
                min[dim1] = obs_ptr[dim1][dim2];
            }
            if (max[dim1]<obs_ptr[dim1][dim2]){
                max[dim1] = obs_ptr[dim1][dim2];
            }
        }
        
    }
    
    // the histogram
    int **hist;
    hist = malloc(nfactors * sizeof(int *));
    
    int bin;
    double *width;
    width = malloc(nfactors * sizeof(double));
    
    // initializing the histograms and differnt bin widths
    for (int dim1=0; dim1<nfactors; ++dim1) {
        hist[dim1]= calloc(numbins, sizeof(int));
        width[dim1] = (max[dim1]-min[dim1])/(double)numbins;
        
        
        // doing the actual binning
        for (int dim2=0; dim2<numberofframes; ++dim2){
            bin = floor((obs_ptr[dim1][dim2]-min[dim1])/width[dim1]);
            if (bin == numbins){
                bin-=1;
            }
            else if ((bin > numbins) || (bin < 0)){
                fprintf(stderr, "The bin is out of range, bin=%i \n", bin);
                return false;
            }
            
            // raise the counter in this bin
            hist[dim1][bin]++;
        }
    }
    
    
    // printing the histogram to file
    fprintf(file_ptr, "### bin -- counter -- normalized \n");
    
    
    for (int dim1=0; dim1<numbins; dim1++){
        for (int dim2=0; dim2<nfactors; ++dim2) {
            
            fprintf(file_ptr, "%f %i %f   ", (min[dim2] + (0.5+(double)dim1)*width[dim2]), hist[dim2][dim1], ((double)hist[dim2][dim1]/((double)numbins)*width[dim2]) );
        }
        fprintf(file_ptr, "\n");
        
    }
    
    free(min);
    free(max);
    free(width);
    for (int dim1=0; dim1<nfactors; ++dim1) {
        free(hist[dim1]);
    }
    free(hist);
    
    
    return true;
    
}




//=========================================================================================================
// convergence,
//      this is programmed on the basis that numframes is always a multiple of at least 100!



// this is for normal variables like pT or distances
void convergence(OBSERVABLE *obs_ptr, double ntor, unsigned long convpoint, FILE *file_ptr){
    
    double sum[3]={0.0};
    // 0: Rf
    // 1: rgyr
    // 2: pT
    
    fprintf(file_ptr, "%lu ", convpoint);
    
    for (int dim1=0; dim1<10; ++dim1){
        
        sum[0] += obs_ptr->err_ee2[dim1];
        sum[1] += obs_ptr->err_rgyr[dim1];
        sum[2] += obs_ptr->err_pt[dim1];
        
    }
    
    fprintf(file_ptr, "%f %f %f \n", sqrt(sum[0]/(double)convpoint), sqrt(sum[1]/(double)convpoint), sum[2]/(ntor*(double)convpoint));
    fflush(stdout);
}




void weighted_convergence(OBSERVABLE *obs_ptr, int num_eps, double (**weights_ptr), double ntor, unsigned long convpoint, FILE *file_ptr){
    
    double sum[3];
    // 0: Rf
    // 1: rgyr
    // 2: pT
    double weight;
    
    fprintf(file_ptr, "%lu ", convpoint);
    
    for (int dim1=0; dim1<num_eps; ++dim1){
        
        sum[0] = 0.0;
        sum[1] = 0.0;
        sum[2] = 0.0;
        weight = 0.0;
        
        for (int dim2=0; dim2<10; ++dim2) {
            
            sum[0] += obs_ptr[dim1].err_ee2[dim2];
            sum[1] += obs_ptr[dim1].err_rgyr[dim2];
            sum[2] += obs_ptr[dim1].err_pt[dim2];
            
            weight += weights_ptr[dim1][dim2];
        }
        
        fprintf(file_ptr, "%f %f %f ", sqrt(sum[0]/weight), sqrt(sum[1]/weight), sum[2]/(ntor*weight));
    }
    
    fprintf(file_ptr, "\n");
    fflush(stdout);
}



// running sum over boltzman factors, a jump indicates optically an outlier! Best is a smooth line.

void weights_growth(double (**obs_ptr), int num_eps, unsigned long convpoint, FILE *file_ptr){
    
    double sum;
    
    fprintf(file_ptr, "%lu ", convpoint);
    
    for (int dim1=0; dim1<num_eps; ++dim1) {
        
        sum = 0.0;
        
        for (int dim2=0; dim2<10; ++dim2){
            
            sum += obs_ptr[dim1][dim2];
        }
        
        fprintf(file_ptr, "%e ", sum);
    }
    
    fprintf(file_ptr, "\n");
    
    fflush(stdout);
    
}


