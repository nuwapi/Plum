#ifndef SRC_FORCE_FIELD_POTENTIAL_TRUNCATED_LJ_H_
#define SRC_FORCE_FIELD_POTENTIAL_TRUNCATED_LJ_H_

#include <iostream>
#include <map>
#include <string>

#include "../molecules/bead.h"
#include "potential_pair.h"

using namespace std; 

/** The Lennard-Jones interaction pair potential truncated at a given cutoff.
    It uses arithmetic mixing for sigmas and the geometric mixing for epsilons.
    */
class PotentialTruncatedLJ : public PotentialPair {
 private:
  /** The cutoff for Lennard-Jones potential. */
  double lj_cutoff;
  /** Sigma values for different types of beads. */
  map<string, double> sigmas; 
  /** Epsilon values for different types of beads. */
  map<string, double> epsilons; 

 public:
  /////////////////////
  // Initialization. //
  /////////////////////
  PotentialTruncatedLJ(string);
  void ReadParameters();

  ///////////////////////
  // Energy and force. //
  ///////////////////////
  double PairEnergy(Bead&, Bead&, double[], int);   
  double PairForce(Bead&, Bead&, double[], int);

}; 

#endif


