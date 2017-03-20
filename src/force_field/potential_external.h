#ifndef SRC_FORCE_FIELD_POTENTIAL_EXTERNAL_H_
#define SRC_FORCE_FIELD_POTENTIAL_EXTERNAL_H_

#include <vector>
#include <map> 
#include <iostream>
#include <iomanip> 
#include <numeric> 

#include "../molecules/bead.h"
#include "../molecules/molecule.h" 

using namespace std;

/** "External potential" is only used for z-direction confined systems and it
    includes uniform wall potentials, such as uniform Lennard-Jones wall,
    charged wall or hard wall. */
class PotentialExternal{
 private:
  string name;
  map<int, double> current_energy_map; 
  map<int, double> trial_energy_map; 
  double E_tot;
  double dE;

 public:
  /////////////////////
  // Initialization. //
  /////////////////////
  PotentialExternal(string);
  virtual void ReadParameters() = 0;
  void EnergyInitialization(vector < Molecule >& mols, double[], int);

  ///////////////////////////////
  // Compute energy and force. //
  ///////////////////////////////
  virtual double BeadEnergy(Bead&, double[]) = 0;
  virtual double BeadForceOnWall(Bead&, double[]) = 0;

  /////////////////////////////////
  // Reading and storing energy. //
  /////////////////////////////////
  void SetE(int, int, double);
  double GetE(int, int);
  void SetEBothMaps(int, double);
  double GetTotalEnergy();
  double CalcTrialTotalEnergy(vector < Molecule >& mols, double[], int);

  ///////////////////
  // MC utilities. //
  ///////////////////
  void EnergyInitForLastMol(vector < Molecule >& mols, int, double, double[], int);
  void AdjustEnergyUponMolDeletion(vector < Molecule >& mols, int);
  double EnergyDifference(vector < Molecule >& mols, int, double[], int);
  void FinalizeEnergyBothMaps(vector < Molecule >& mols, int, bool);

  ////////////
  // Other. //
  ////////////
  string PotentialName();

};

#endif


