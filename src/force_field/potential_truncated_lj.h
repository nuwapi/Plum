/** The basic Lennard-Jones interaction pair potential truncated at a given
    cutoff, it also uses arithmetic mixing for sigmas and the  geometric mixing
    for epsilons.
  */

#ifndef SRC_FORCE_FIELD_POTENTIAL_TRUNCATED_LJ_H_
#define SRC_FORCE_FIELD_POTENTIAL_TRUNCATED_LJ_H_

#include <iostream>
#include <map>
#include <string>

#include "../molecules/bead.h"
#include "potential_pair.h"

using namespace std; 

class PotentialTruncatedLJ : public PotentialPair {
 private:
  /** The cutoff for Lennard-Jones potential. */
  double lj_cutoff;
  /** The sigma values for the beads in each chain type. */
  map<string, double> sigmas; 
  /** The epsilon values for the beads in each chain type. */
  map<string, double> epsilons; 

 public: 
  // Initialization functions.
  PotentialTruncatedLJ(string);
  /** Read parameters from file, can read in the sigmas and epsilons from
      multiple chain types (symbols). */
  void ReadParameters();

  // Energy functions.
  /** Energy between two beads. */
  double PairEnergy(Bead&, Bead&, double[], int);   

}; 

#endif


