//
//  entropy.h
//  
//
//  Created by Johannes Dietschreit on 11.06.15.
//
//

#ifndef ____entropy__
#define ____entropy__

double intra_entropy(double (*boltz_ptr), double (*boltzener_ptr), double logatt);

double intra_entropy_error_ten(double entropy, double (*boltz_ptr), double (*boltzener_ptr), double (*logatt));

#endif /* defined(____entropy__) */
