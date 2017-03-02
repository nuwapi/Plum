#ifndef SRC_FORCE_FIELD_POTENTIAL_EWALD_H_
#define SRC_FORCE_FIELD_POTENTIAL_EWALD_H_

#include <iostream>
#include <map>
#include <vector> 

#include "../molecules/bead.h"
#include "../molecules/molecule.h"

using namespace std; 

/** Note: 1. repl = reciprocal*/
class PotentialEwald {
 private:
  string name;
  // Energy maps.
  /** Real part of the Ewald summation of the current configuration. */
  map<pair<int,int>, double> current_real_energy_map;
  //// The reciprocal part could be stored in a 1D array. To be optimized in the
  //// future.
  /** Reciprocal part of the Ewald summation of the current configuration. */
  map<pair<int,int>, double> current_repl_energy_map;
  /** Self part of the Ewald summation of the current configuration. */
  map<int,           double> current_self_energy_map;
  // Surface correction of the current configuration.
  // double current_J;
  /** Real part of the Ewald summation of the MC trial configuration. */
  map<pair<int,int>, double> trial_real_energy_map;
  /** Reciprocal part of the Ewald summation of the MC trial configuration. */
  map<pair<int,int>, double> trial_repl_energy_map;
  /** Self part of the Ewald summation of the MC trial configuration. */
  map<int,           double> trial_self_energy_map;

  map<pair<int,int>, double> current_real_pphi_map;
  map<pair<int,int>, double> current_repl_pphi_map;
  map<pair<int,int>, double> trial_real_pphi_map;
  map<pair<int,int>, double> trial_repl_pphi_map;

  /** Temporarily store the energies for the trial chain. */
  double * trial_chain_e;

  // Energy sums.
  double E_tot;
  double current_real_E;
  double current_repl_E;
  double current_self_E;
  double current_dipl_E;
  double trial_dipl_E;
  double dE;

  bool calc_pphi;

 protected:
  /** The volume of the simulation unit cell, in unit length. */
  double box_vol;
  double box_l[3];
  /** When there is a confining wall, we also use dipole correction. */
  bool dipole_correction;

 public:
  PotentialEwald(string, double[3]);
  ~PotentialEwald();
  virtual void ReadParameters() = 0;

  // Different manipulations of energies.
  // Part 1.
  /** Get real space pair energy. */
  virtual double PairEnergyReal(Bead&, Bead&, int) = 0;
  /** Get reciprocal space pair energy. */
  virtual double PairEnergyRepl(Bead&, Bead&, int) = 0;
  /** Get the self energy for each bead. */
  virtual double SelfEnergy(Bead&) = 0;
  // For pressure calculations.
  virtual double PairDForceReal(Bead&, Bead&, Bead&, Bead&, int) = 0;
  virtual double PairDForceRepl(Bead&, Bead&, Bead&, Bead&, int) = 0;
  virtual double DipoleE(vector<Molecule>&) = 0;
  virtual double DipoleE(vector<Molecule>&, vector<Bead>&) = 0;
  virtual double DipoleE(vector<Molecule>&, int, int) = 0;
  // Partial U partial V.
  double PUPV(vector<Molecule>&, double, int);
  double RDotF(vector<Molecule>&, double, int);

  /** Set energy between a specific pair in the designated energy maps.\n
      Input arguments:\n
        int - 0 sets current real and repl maps, 1 sets the trial ones.\n
        int - Index for bead 1.\n
        int - Index for bead 2.\n
        double - The real energy value.\n
        double - The repl energy value. */
  void SetERealRepl(int, int, int, double, double);
  void SetPPhiRealRepl(int, int, int, double, double);
  /** Similar to above, set self-energy in the designated map. */
  void SetESelf(int, int, double);
  /** Set energy between a specific pair in both trial and current, real and
      repl energy maps. */
  void SetEBoth2DMaps(int, int, double, double);
  void SetPPhiBoth2DMaps(int, int, double, double);
  /** Set energy between a specific pair in both trial and current self-energy
      maps. */
  void SetEBothSelfMaps(int, double);
  /** Get real energy + repl energy between a specific pair in the designated
      energy map. */
  double GetERealRepl(int, int, int);
  /** Get the self-energy for a bead. */
  double GetESelf(int, int);
  /** Return total energy. */
  double GetTotalEnergy();
  // Part 2.
  /** Fill in all energy maps with values for initial configuration. */
  void EnergyInitialization(vector<Molecule>&, int); 
  /** Calculate the dE due to a MC move. */
  double EnergyDifference(vector<Molecule>&, int, int);
  /** Update all energy maps after the decision of a MC move is made. */
  void FinalizeEnergyBothMaps(vector<Molecule>&, int, bool);

  /** Calculate the energy between the CBMC trial chain and the rest of the
      system. */
  double TrialChainEnergy(vector<Molecule>&, vector<Bead>&, int, int, int);
  /** Initialize a new molecule (requires that the beads have proper IDs).
      Assumes mol is inserted at end of mol array! */
  void EnergyInitForLastMol(vector<Molecule>&, int, double, int);
  /** Erase the energy of the deleted molecule. */
  void AdjustEnergyUponMolDeletion(vector<Molecule>&, int); 
  string PotentialName();

  /** This is a bad inefficient function, only for debug purposes. */
  void UpdateEnergyComponents(vector<Molecule>&, int);
  /** Only for debug purposes. */
  double GetRealEnergy();
  /** Only for debug purposes. */
  double GetReplEnergy();
  /** Only for debug purposes. */
  double GetSelfEnergy();
  virtual double GetlB() = 0;

}; 

#endif


