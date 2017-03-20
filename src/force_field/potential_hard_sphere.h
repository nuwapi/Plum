#ifndef SRC_FORCE_FIELD_POTENTIAL_HARD_SPHERE_H_
#define SRC_FORCE_FIELD_POTENTIAL_HARD_SPHERE_H_

#include <iostream>
#include <map>
#include <string>

#include "../molecules/bead.h"
#include "potential_pair.h"

using namespace std; 

/** Hard sphere potential is self-explanatory. */
class PotentialHardSphere : public PotentialPair {
 private:
  /** The radius for different types of beads. */
  map<string, double> radii; 

 public: 
  /////////////////////
  // Initialization. //
  /////////////////////
  PotentialHardSphere(string);
  void ReadParameters();

  ///////////////////////
  // Energy and force. //
  ///////////////////////
  double PairEnergy(Bead&, Bead&, double[], int);   
  double PairForce(Bead&, Bead&, double[], int);

}; 

#endif


