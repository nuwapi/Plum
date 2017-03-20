/** Hard sphere potential is 0 if there is no bead-bead overlap, infinite
    otherwise. */

#include "potential_hard_sphere.h"

#include <iomanip>
#include <vector>

#include "../molecules/molecule.h"
#include "../utilities/constants.h"

PotentialHardSphere::PotentialHardSphere(string potential_name) :
                     PotentialPair(potential_name) {
  ReadParameters(); 

}

void PotentialHardSphere::ReadParameters() {
  cout << setw(35) << "[PP] Pair potential type    : " << "Hard sphere" << endl;

  string flag;
  string symbol;
  double radius;

  while (true) {
    cin >> flag >> symbol;
    if (symbol == "end")  break;
    else {
      cin >> flag >> radius;
      radii[symbol] = radius;
      cout << setw(35) << "[PP] Bead type and r (ul)   : "  << symbol << " - "
           << radius << endl;
    }
  }

}

double PotentialHardSphere::PairEnergy(Bead& bead1, Bead& bead2,
                                       double box_l[], int npbc) {
  double r = bead1.BBDist(bead2, box_l, npbc);
  double allowed_r = radii[bead1.Symbol()] + radii[bead2.Symbol()];
  if (r <= allowed_r) {
    return kVeryLargeEnergy;
  }
  else {
    return 0;
  }

}

double PotentialHardSphere::PairForce(Bead& bead1, Bead& bead2,
                                      double box_l[], int npbc) {
  // Force for hard sphere potential is not well-defined. This function should
  // not be use to generate meaningful data.
  return 0;

}


