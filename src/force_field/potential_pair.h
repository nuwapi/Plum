//Base object for pair potential. This must be defined for every simulation. 
//To avoid recalculating after each move (very expensive)
//exisiting energies are kept in a map with keys consisting of a pair containing 
//the two bead IDs. 

#ifndef PAIRPOT_H
#define PAIRPOT_H

#include <string> 
#include <vector> 
#include <map> 
#include <iostream> 
#include <numeric> 
#include <iomanip> 

#include "../molecules/bead.h"
#include "../molecules/molecule.h"
#include "../utilities/misc.h"

using namespace std; 

class PotentialPair {
 private:
  string name;
  map<pair<int,int>, double> current_energy_map;
  map<pair<int,int>, double> trial_energy_map;
  double E_tot;
  double dE;

 public:
  // Initialization functions.
  PotentialPair(string);
  virtual void ReadParameters() = 0;

  // Different manipulations of energies.
  /** Get pair energy. */
  virtual double PairEnergy(Bead&, Bead&, double[], int) = 0;
  /** Set energy between a specific pair in the designated energy map. */
  void SetE(int, int, int, double);
  /** Set energy between a specific pair in both energy maps. */
  void SetEBothMaps(int, int, double);
  /** Get energy between a specific pair in the designated energy map. */
  double GetE(int, int, int);
  /** Return total energy. */
  double GetTotalEnergy();
  /** Fill in all energy maps with values for initial configuration.\n
      Assume that the unit cell is big enough such that the beads don't interact
      with its periodic image. */
  void EnergyInitialization(vector<Molecule>&, double[], int); 
  /** Initialize a new molecule (requires that the beads have proper IDs).
      Assumes mol is inserted at end of mol array! */
  void EnergyInitForLastMol(vector<Molecule>&, int, double, double[], int);
  /** Calculate the dE due to a MC move. */
  double EnergyDifference(vector<Molecule>&, int, double[], int);
  /** Update all energy maps after the decision of a MC move is made. */
  void FinalizeEnergyBothMaps(vector<Molecule>&, int, bool);
  /** Erase the energy of the deleted molecule. */
  void AdjustEnergyUponMolDeletion(vector<Molecule>&, int); 
  string PotentialName();

}; 

#endif


