#ifndef SRC_FORCE_FIELD_POTENTIAL_HARD_SPHERE_H_
#define SRC_FORCE_FIELD_POTENTIAL_HARD_SPHERE_H_

#include <iostream>
#include <map>
#include <string>

#include "../molecules/bead.h"
#include "potential_pair.h"

using namespace std; 

class PotentialHardSphere : public PotentialPair {
 private:
  /** The radius for the beads in each chain type. */
  map<string, double> radii; 

 public: 
  // Initialization functions.
  PotentialHardSphere(string);
  /** Read parameters from file, can read in the sigmas and epsilons from
      multiple chain types (symbols). */
  void ReadParameters();

  // Energy functions.
  /** Energy between two beads. */
  double PairEnergy(Bead&, Bead&, double[], int);   

}; 

#endif


