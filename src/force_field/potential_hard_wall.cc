#include "potential_hard_wall.h" 

#include "../utilities/constants.h"

using namespace std; 

PotentialHardWall::PotentialHardWall(string potential_name)
                        : PotentialExternal(potential_name) {
  ReadParameters(); 

}

void PotentialHardWall::ReadParameters() {
  cout << setw(35) << "[eP] External potential type: " << "Hard wall"
       << endl;

  string flag;
  string symbol;
  double radius;

  while (true) {
    cin >> flag >> symbol;
    if (symbol == "end") {
      break;
    }
    else {
      cin >> flag >> radius;
      radii[symbol] = radius;
      cout << setw(35) << "[eP] Radius for bead type   : " << symbol << " - "
           << radius << endl;
    }
  }

}

double PotentialHardWall::BeadEnergy(Bead& bead, double box_l[]) {
  double energy = 0;
  // Get the z coordinate from the trial coordinate array.
  double z = bead.GetCrd(1, 2);
  // If the bead falls outside of the box in the z-direction.
  if (z <= radii[bead.Symbol()] || z+radii[bead.Symbol()] >= box_l[2]) {
    energy = kVeryLargeEnergy;
  }

  return energy;

}

double PotentialHardWall::BeadForceOnWall(Bead& bead, double box_l[]) {
  // 0 force since this is a hard wall.
  return 0; 

}


