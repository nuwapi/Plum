#ifndef SRC_FORCE_FIELD_POTENTIAL_HARD_WALL_H_
#define SRC_FORCE_FIELD_POTENTIAL_HARD_WALL_H_

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

/** This external potential places two walls perpendicular to the z-axis one at
    z=0, the other at z=box_l[2]. Bead energy is 0 if it is not overlapping with
    the walls and infinite otherwise. */
class PotentialHardWall : public PotentialExternal {
 private:
  /** The radius for different bead types. */
  map<string, double> radii;

 public:
  /////////////////////
  // Initialization. //
  /////////////////////
  PotentialHardWall(string); 
  void ReadParameters(); 

  ///////////////////////
  // Energy and force. //
  ///////////////////////
  double BeadEnergy(Bead&, double[]);
  double BeadForceOnWall(Bead&, double[]); 

};

#endif


