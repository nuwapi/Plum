/** The repulsive Lennard-Jones and full Lennard-Jones are implemented together
    through the use of cutoff values: 1) If cutoff < 0, use repulsive Lennard-
    Jones. 2) If cutoff is non-negative, use it as the cutoff for full Lennard-
    Jones. */

#include "potential_truncated_lj.h"

#include <iomanip>
#include <vector>

#include "../molecules/molecule.h"
#include "../utilities/constants.h"

PotentialTruncatedLJ::PotentialTruncatedLJ(string potential_name) :
                      PotentialPair(potential_name) {
  ReadParameters();

}

void PotentialTruncatedLJ::ReadParameters() {
  cout << setw(35) << "[PP] Pair potential type    : " << "Truncated LJ" 
       << endl;

  string flag; 
  string symbol;
  double sigma;
  double epsilon;

  cin >> flag >> lj_cutoff;
  while (true) {
    cin >> flag >> symbol;
    if (symbol == "end") {
      break;
    }
    else {
      cin >> flag >> sigma >> flag >> epsilon;
      sigmas[symbol] = sigma; 
      epsilons[symbol] = epsilon; 
      cout << setw(35) << "[PP] Sigma for bead type    : " << symbol << " - "
           << sigma << endl;
      cout << setw(35) << "[PP] Epsilon for bead type  : " << symbol << " - "
           << epsilon << endl;
    }
  }

}

double PotentialTruncatedLJ::PairEnergy(Bead& bead1, Bead& bead2,
                                        double box_l[], int npbc) {
  double energy = 0;
  double r = 0;
  r = bead1.BBDist(bead2, box_l, npbc);

  double sigma = (sigmas[bead1.Symbol()] + sigmas[bead2.Symbol()]) / 2;
  double epsilon = sqrt((epsilons[bead1.Symbol()] * epsilons[bead2.Symbol()]));
  double r6;

  if (r > 0)
    r6 = pow((sigma/r), 6);
  else
    r6 = -1;

  if (r6 == -1) {
    energy = kVeryLargeEnergy;
  }
  // Repulsive Lennard-Jones.
  else if (lj_cutoff < 0) {
    if (r < sigma) {
      energy = 4 * epsilon * (r6*r6 - r6);
    }
  }
  // Full Lennard-Jones.
  else if (r < lj_cutoff) {
    energy = 4 * epsilon * (r6*r6 - r6);
  }

  return energy;

}

double PotentialTruncatedLJ::PairForce(Bead& bead1, Bead& bead2,
                                       double box_l[], int npbc) {
  double force = 0;
  double r = 0;
  r = bead1.BBDist(bead2, box_l, npbc);

  double sigma = (sigmas[bead1.Symbol()] + sigmas[bead2.Symbol()]) / 2;
  double epsilon = sqrt((epsilons[bead1.Symbol()] * epsilons[bead2.Symbol()]));
  double r6;

  if (r > 0)
    r6 = pow((sigma/r), 6);
  else
    r6 = -1;

  if (r6 == -1) {
    force = kVeryLargeEnergy;
  }
  // Repulsive Lennard-Jones.
  else if (lj_cutoff < 0) {
    if (r < sigma) {
      force = 4 * epsilon * (r6*r6*12/r - r6*6/r);
      // With this repulsive LJ formula, the force at the cutoff point is not
      // well-defined. This should be adjusted in the future.
    }
  }
  // Full Lennard-Jones.
  else if (r < lj_cutoff) {
    force = 4 * epsilon * (r6*r6*12/r - r6*6/r);
  }

  return force;

}


