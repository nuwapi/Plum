#include "potential_well_wall.h" 

#include "../utilities/constants.h"

using namespace std; 

PotentialWellWall::PotentialWellWall(string potential_name)
                        : PotentialExternal(potential_name) {
  ReadParameters(); 

}

void PotentialWellWall::ReadParameters() {
  cout << setw(35) << "External potential type     : "
       << "Square well wall" << endl;

  string flag;
  string symbol;
  double radius;

  cin >> flag >> well_width;
  cin >> flag >> well_depth;
  cout << setw(31) << "Well witdh in unit length   = " << well_width << endl;
  cout << setw(31) << "Well depth in kT            = " << well_depth << endl;

  while (true) {
    cin >> flag >> symbol;
    if (symbol == "end") {
      break;
    }
    else {
      cin >> flag >> radius;
      radii[symbol] = radius;
      cout << setw(31) << "Radius for bead type    " << setw(3) << symbol
           << " = " << radius << endl;
    }
  }

}

double PotentialWellWall::BeadEnergy(Bead& bead, double box_l[]) {
  double energy = 0;
  // Get the z coordinate from the trial coordinate array.
  double z = bead.GetCrd(1, 2);
  // If the bead falls outside of the box in the z-direction.
  if (z <= radii[bead.Symbol()] || z >= box_l[2] - radii[bead.Symbol()]) {
    energy = kVeryLargeEnergy;
  }
  else if (z < well_width || z >  box_l[2] - well_width)  {
    energy = well_depth;
  }

  return energy;

}

// Calculate force PER UNIT AREA due to molecule-wall interaction.
double PotentialWellWall::CalculateForce(vector < Molecule >& mols, double box_l[]){
  // 0 force since this is a hard wall.
  return 0;

}


