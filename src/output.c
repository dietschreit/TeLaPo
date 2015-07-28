//
//  output.c
//  
//
//  Created by Johannes Dietschreit on 29.05.15.
//
//

#include <stdio.h>
#include <stdarg.h>
#include "output.h"


void log_out(FILE *logfile, const char * fmt, ... ){
    
    va_list args1;                   // get all arguments from the variable list
    va_list args2;
    
    va_start(args1, fmt);            //initialiazes the variable argument list
    va_start(args2, fmt);
    
    vprintf(fmt, args1);             //prints the va-list to screen
    vfprintf(logfile, fmt, args2);   //prints the va-list to file
    
    
    va_end(args1);
    va_end(args2);                   //ends va-list
    
    fflush(logfile);
    fflush(stdout);
    
}

