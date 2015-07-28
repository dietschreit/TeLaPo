//
//  SASA.h
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#ifndef ____SASA__
#define ____SASA__

double pi;

double get_asa(double (**chain), unsigned int numberofbeads, double r);

double delta_asa(double (**chain), unsigned int numberofbeads, double b);

#endif /* defined(____SASA__) */
