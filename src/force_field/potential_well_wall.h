/** This external potential places two walls perpendicular to the z-axis one at
    z=0, one at z=box_l[2] (z-dimension). Each wall has the same
    Lennard-Jones parameters and cutoff. The potentials are truncated at the
    cutoff distance from the corresponding wall. The LJ wall potential uses the
    mixing rules the same as the LJ pair potential. 
*/

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
  // Initialization functions.
  PotentialWellWall(string); 
  void ReadParameters(); 

  // Energy and force.
  /** Calculate the energy between a bead and the LJ walls. */
  double BeadEnergy(Bead&, double[]);
  double CalculateForce(vector < Molecule >&, double[]); 

};

#endif


