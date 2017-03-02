//A base class for bonded potentials
//Derived classes may contain any combination of bond, angle, and dihedral energy funcs
//May be any potential that depends only upon the positions of beads within a single molecule
//Exactly one PotentialBond must be defined in a given simulation, even if no bonds are listed in the 
//topology file 

#ifndef BONDPOT_H 
#define BONDPOT_H

#include <iostream>
#include <iomanip> 
#include <numeric> 
#include <random> 
#include <vector>

#include "../molecules/molecule.h" 

using namespace std; 

class PotentialBond{
 private:
  string name;
  /** The bond energy array is save for one molecule at a time. 
      int is molecule ID instead of bond ID. */
  vector<double> current_energy_array;
  /** The bond energy array is save for one molecule at a time. */
  vector<double> trial_energy_array;

 protected: 
  double E_tot;
  double dE;

 public:
  PotentialBond(int, string);
  void EnergyInitialization(vector <Molecule>&, double[], int);  // args are mol array, length, npbc
  void EnergyInitForLastMol(vector <Molecule>&, double[], int);  // init one new mol being added 
  void AdjustEnergyUponMolDeletion(int);  // delete stored energies for a molecule being deleted 
  void SetE(int, int, double);  // same array, index, value args 
  double GetE(int, int);
  void FinalizeEnergy(int, bool); 
  double GetTotalEnergy(); 
  string PotentialName();

  virtual double MoleculeEnergy(Molecule&, double[], int) = 0; 
  virtual double EnergyDifference(vector<Molecule>&, double[], int, int) = 0; 
  virtual void ReadParameters() = 0;
  virtual double RandomBondLen(double, mt19937&) = 0; //for growing mols for gc insertion
  /** Returning the equilibrium bond length. */
  virtual double EqBondLen() = 0;
  double CalcTrialTotalEnergy(vector<Molecule>&, double[], int); //for use in pressure calculation 

}; 

#endif


