#include "potential_truncated_lj.h"

#include <iomanip>
#include <vector>

#include "../molecules/molecule.h"
#include "../utilities/constants.h"

PotentialTruncatedLJ::PotentialTruncatedLJ(string potential_name)
                                  : PotentialPair(potential_name) {
  ReadParameters();

}

void PotentialTruncatedLJ::ReadParameters() {
  cout << setw(35) << "Pair potential type         : " << "Truncated LJ" 
       << endl;

  string flag; 
  string symbol;
  double sigma;
  double epsilon;

  // Read simple pot parameters.
  cin >> flag >> lj_cutoff;

  // Read sigma + epsilon associated with each symbol.
  while (true) {
    cin >> flag >> symbol;
    if (symbol == "end") {
      break;
    }
    else {
      cin >> flag >> sigma
          >> flag >> epsilon; 
      sigmas[symbol] = sigma; 
      epsilons[symbol] = epsilon; 
      cout << setw(31) << "Sigma for bead type     " << setw(3) << symbol 
           << " = " << sigma << endl;
      cout << setw(31) << "Epsilon for bead type   " << setw(3) << symbol
           << " = " << epsilon << endl;
    }
  }

}

double PotentialTruncatedLJ::PairEnergy(Bead& bead1, Bead& bead2,
                                        double box_l[], int npbc) {
  double energy = 0;
  double r = 0;
  r = bead1.BBDist(bead2, box_l, npbc);

  // Arithmetic mixing.
  double sigma = (sigmas[bead1.Symbol()] + sigmas[bead2.Symbol()]) / 2;
  // Geometric mixing.
  double epsilon = sqrt((epsilons[bead1.Symbol()] * epsilons[bead2.Symbol()]));
  double r6;

  if (r > 0)
    r6 = pow((sigma/r), 6);
  else
    r6 = -1;

  // Use purely repulsive potential.
  if (r6 == -1) {
    energy = kVeryLargeEnergy;
  }
  else if (lj_cutoff < 0) {
    if (r < sigma) {
      energy = 4 * epsilon * (r6*r6 - r6);
    }
    // Else, energy = 0.
  }
  // Use normal LJ.
  else if (r < lj_cutoff) {
    energy = 4 * epsilon * (r6*r6 - r6);
  }

  return energy;

}


