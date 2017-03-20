#ifndef SRC_FORCE_FIELD_POTENTIAL_WELL_WALL_H_
#define SRC_FORCE_FIELD_POTENTIAL_WELL_WALL_H_

#include <iostream>
#include <iomanip>
#include <numeric>
#include <string>
#include <map>
#include <vector>

#include "../molecules/bead.h"
#include "../molecules/molecule.h"
#include "potential_external.h"

using namespace std;

/** Square well + hard wall. */
class PotentialWellWall : public PotentialExternal {
 private:
  /** The diameter of the monomer. */
  double well_width;
  /** Has to be negative for an attractive well. */
  double well_depth;
  /** The radius for the beads in each chain type. */
  map<string, double> radii;

 public: 
  /////////////////////
  // Initialization. //
  /////////////////////
  PotentialWellWall(string); 
  void ReadParameters(); 

  ///////////////////////
  // Energy and force. //
  ///////////////////////
  double BeadEnergy(Bead&, double[]);
  double BeadForceOnWall(Bead&, double[]);

};

#endif


