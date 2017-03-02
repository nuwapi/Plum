//Base class for external potentials. 
//"external potential" is any pot where the energy for a bead
//depends only on the bead's own coordinates. 
//useful for hard walls, electric fields, etc. 
//the extpot is optional; only one may be defined per simulation
#ifndef EXTPOT_H
#define EXTPOT_H

#include <vector>
#include <map> 
#include <iostream>
#include <iomanip> 
#include <numeric> 
#include "../molecules/bead.h"
#include "../molecules/molecule.h" 
using namespace std;

class PotentialExternal{
 private:
  string name;
  map<int, double> current_energy_map; 
  map<int, double> trial_energy_map; 
  double E_tot;
  double dE;

 public: 
  PotentialExternal(string);
  void SetE(int, int, double);
  double GetE(int, int);
  void SetEBothMaps(int, double);
  void EnergyInitialization(vector < Molecule >& mols, double[], int);
  void EnergyInitForLastMol(vector < Molecule >& mols, int, double, double[], int);
  void AdjustEnergyUponMolDeletion(vector < Molecule >& mols, int);
  double EnergyDifference(vector < Molecule >& mols, int, double[], int);
  void FinalizeEnergyBothMaps(vector < Molecule >& mols, int, bool);
  double GetTotalEnergy();
  double CalcTrialTotalEnergy(vector < Molecule >& mols, double[], int);
  string PotentialName();

  virtual void ReadParameters() = 0;
  virtual double BeadEnergy(Bead&, double[]) = 0;
  virtual double CalculateForce(vector < Molecule >&, double[]) = 0;
  
};

#endif


