/** Useful small functions used multiple places in program. */

#ifndef SRC_UTILITIES_MISC_H_
#define SRC_UTILITIES_MISC_H_ 

#include <random> 
#include <string>

#include "../molecules/bead.h"

/** Bead-bead dist func that takes bead ptrs. from MC12 & earlier. */
double getDist(Bead&, Bead&, double[], int);
/** Bead-bead dist func that takes refs to beads. this is the MC13 way. */
double getDist(Bead*, Bead*, double[], int);
/** The distance vector from bead1 to bead2. The vector components will be
    stored in double (&)[3]. */
void GetDistVector(Bead&, Bead&, double[], int, double (&)[3]);
/** Same function using the current position. */
void GetDistVectorC(Bead&, Bead&, double[], int, double (&)[3]);
/** Move all atoms to the positive central cell and do NOT apply minimum
    image convention. */
void GetDistVectorConsistent(Bead&, Bead&, double[], int, double (&)[3]);
/** Fills in a vector with random point on unit sphere. */
void randSphere(double[], std::mt19937&);
/** Returns rand deviate from gaussian distribution with mean 0 and stdev 1.
    multiplied /shifted to give the requested distribution. */
double gasdev(double, double, std::mt19937&);
/** Print bool as yes or no. */
string YesOrNo(bool);
int factorial(int);
/** This is a simple least square fit procedure. */
double Interpolate(double *, double *, int, double);
double Interpolate2(double *, double *, int, double);

#endif


