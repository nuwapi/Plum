#ifndef SRC_UTILITY_CONSTANTS_H_
#define SRC_UTILITY_CONSTANTS_H_

/** A very small number. */
const double kVerySmallNumber = 1E-12;
/** A small number. */
const double kSmallNumber = 1E-6;
/** A medium small number. */
const double kMedSmallNumber = 1E-4;
/** For Ewald sum. */
const double kEwaldCutoff = 1E-7;
/** In unit length. */
const double kVerySmallDistance = 1E-8;
/** In kBT. */
const double kVeryLargeEnergy = 1E+8;
/** For detecting adsorption, in unit length. */
const double kAdsorbCutoff = 1.13;  // 1.5 is slightly larger than r_m.

/** Constant pi. */
const double kPi = 3.14159265359;
/** kB in J/K. */
const double kKB = 1.38064853E-23;
/** kBT at T = 273.15K, in J. 
    1.38064853E-23 * 273.15 J/K * K
*/
const double kKBT0 = 3.7712415E-21;
/** Coulomb's constant, 1/(4*pi*eps0), in J*A/e2
    Derivation:
      8.9875517873681764E+9  (N*m2/C2 = kg*m3/s2C2 = J*m/C2)
    =(8.9875517873681764E+9)*(1E+10)/(6.24150934E+18)^2 (J*A/e2)
 */
const double kKe = 2.3070774E-18;
/** Coulomb's constant in eV*A/e^2.
    Derivation:
      8.9875517873681764E+9  (N*m2/C2 = kg*m3/s2C2 = J*m/C2)
    =(8.9875517873681764E+9)*(1E+10)/(6.24150934E+18)^2 (J*A/e2)
    =(8.9875517873681764E+9)*(1E+10)*(6.2415093E+18)/(6.24150934E+18)^2(eV*A/e2)
    =14.3996447657 (eV*A/e2)
 */
const double kKe_2 = 14.3996447657;
/** Vacuum permitivity in e2/(eV*A).
    Derivation:
     (1/ke in e2/(eV*A)) * (Eps0 in F/m)/(1/ke in F/m)
    =(1/14.3996447657) * (8.854187817E-12)/(1/(8.9875517873681764E+9))
 */
const double kEps0 = 0.00552634963;

// 1 / (2 * 4 * pi * 0.00552634963 * (1/42))
// = 302.392540167
// this is "water" dielectric constant

/** The number of different Monte Carlo moves implemented in the code. */
const int kNoMoveType = 5;

const double kDz = 0.00001;
/** Scale the box this number of times when using dipole correction. */
const int kDiCorrection = 3;

#endif


