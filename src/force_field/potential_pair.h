#ifndef SRC_FORCE_FIELD_POTENTIAL_PAIR_H_
#define SRC_FORCE_FIELD_POTENTIAL_PAIR_H_

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

/** "Pair potential" are those short-ranged potentials, such as Lennard-Jones
    potential, whose interaction range is not longer than the dimension of the
    simulation box. The dimensions of the simulation box and the strength of
    pair potentials need to be compatible. Typical simulation box sizes will
    usually fulfill this requirement.  */
class PotentialPair {
 private:
  string name;
  map<pair<int,int>, double> current_energy_map;
  map<pair<int,int>, double> trial_energy_map;
  /** Total pair energy of the system. */
  double E_tot;
  /** Pair energy difference upon MC move */
  double dE;

 public:
  /////////////////////
  // Initialization. //
  /////////////////////
  PotentialPair(string);
  virtual void ReadParameters() = 0;
  /** Fill in trial and current energy maps for the initial configuration. */
  void EnergyInitialization(vector<Molecule>&, double[], int);

  ////////////////////////////////////
  // Energy and force calculations. //
  ////////////////////////////////////
  /** Calculate pair energy. */
  virtual double PairEnergy(Bead&, Bead&, double[], int) = 0;
  /** Calculate pair force (scalar as a function of r). */
  virtual double PairForce(Bead&, Bead&, double[], int) = 0;

  ///////////////////////////////////
  // Reading and storing energies. //
  ///////////////////////////////////
  /** Set the energy between a pair in the designated energy map. */
  void SetE(int, int, int, double);
  /** Set the energy between a pair in both trial and current energy maps. */
  void SetEBothMaps(int, int, double);
  /** Get the energy between a pair from the designated energy map. */
  double GetE(int, int, int);
  /** Return total pair energy of the system. */
  double GetTotalEnergy();

  ////////////////////
  // MC utitlities. //
  ////////////////////
  /** Initialize the pair energy for a new molecule for the GCMC routines
      (requires that the beads have proper IDs). This assumes that the molecule
      is inserted at end of molecule array, which is how the GCMC routine does
      it. */
  void EnergyInitForLastMol(vector<Molecule>&, int, double, double[], int);
  /** Erase the energy of the deleted molecule for the GCMC routine. */
  void AdjustEnergyUponMolDeletion(vector<Molecule>&, int);
  /** Calculate the pair dE of the system due to a MC move. */
  double EnergyDifference(vector<Molecule>&, int, double[], int);
  /** Update all energy maps after the decision of a MC move is made. */
  void FinalizeEnergyBothMaps(vector<Molecule>&, int, bool);

  ////////////
  // Other. //
  ////////////
  string PotentialName();

}; 

#endif


