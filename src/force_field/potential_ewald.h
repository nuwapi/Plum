#ifndef SRC_FORCE_FIELD_POTENTIAL_EWALD_H_
#define SRC_FORCE_FIELD_POTENTIAL_EWALD_H_

#include <iostream>
#include <map>
#include <vector> 

#include "../molecules/bead.h"
#include "../molecules/molecule.h"

/** The "ewald potential" here refer to long-ranged potentials that interact
    across multiple periodic boxes. Electrostatic potentials are commonly
    calculated by Ewald summation but techinically, any r^n potential can be
    computed by the Ewald method. */
using namespace std; 

class PotentialEwald {
 private:
  /** Name of the potential. */
  string name;

  ////////////////////////////
  // Energy and force maps. //
  ////////////////////////////
  /** Real energy of the current configuration. */
  map<pair<int,int>, double> current_real_energy_map;
  /** Reciprocal energy of the current configuration. */
  map<pair<int,int>, double> current_repl_energy_map;
  /** Self energy of the current configuration. */
  map<int,           double> current_self_energy_map;
  /** Real energy of the MC trial configuration. */
  map<pair<int,int>, double> trial_real_energy_map;
  /** Reciprocal energy of the MC trial configuration. */
  map<pair<int,int>, double> trial_repl_energy_map;
  /** Self energy of the MC trial configuration. (unchanged) */
  map<int,           double> trial_self_energy_map;
  // The corresponding force maps for the current and the trial configurations.
  map<pair<int,int>, double> current_real_pphi_map;
  map<pair<int,int>, double> current_repl_pphi_map;
  map<pair<int,int>, double> trial_real_pphi_map;
  map<pair<int,int>, double> trial_repl_pphi_map;

  //////////////////
  // Energy sums. //
  //////////////////
  /** Total current energy of the system. */
  double E_tot;
  /** Total current real energy of the system. */
  double current_real_E;
  /** Total current reciprocal energy of the system. */
  double current_repl_E;
  /** Total current self energy of the system. */
  double current_self_E;
  /** Total current dipole correction energy of the system. */
  double current_dipl_E;
  /** Total trial dipole correction energy of the system. */
  double trial_dipl_E;
  /** Total energy change upon MC move. */
  double dE;

  ////////////
  // Other. //
  ////////////
  /** Temporarily store the energies for the trial chain for the GCMC routine to
      avoid recalculation (Actually, this has not been fully utilized, the
      current simulation code is still doing repetative calcultions). */
  double * trial_chain_e;
  /** Decide whether forces are calculated when energies are calculated. */
  bool calc_pphi;

 protected:
  /** The dimesion of the simulation unit cell. They will be the padded
      dimensions if dipole correction is used. This is also why the Ewald
      object stores its own box dimensions. */
  double box_l[3];
  /** The volume of the simulation unit cell. This volume will be the padded
      volume if dipole correction is used. */
  double box_vol;
  /** The scaled volume of the unit cell, used in the volume perturbation method
      to calculate pressure. */
  double box_vol_forP;
  /** Decide whether to use dipole correction. Dipole correction should be used
      for systems confined along the z-direction. */
  bool dipole_correction;

 public:
  /////////////////////
  // Initialization. //
  /////////////////////
  PotentialEwald(string, double[3]);
  ~PotentialEwald();
  virtual void ReadParameters() = 0;

  ///////////////////////
  // Energy and force. //
  ///////////////////////
  /** Compute real pair energy. */
  virtual double PairEnergyReal(Bead&, Bead&, int) = 0;
  /** Compute reciprocal pair energy. */
  virtual double PairEnergyRepl(Bead&, Bead&, int) = 0;
  /** Compute the self energy for each bead. */
  virtual double SelfEnergy(Bead&) = 0;
  /** Compute real pair energy for scaled volume in pressure calculations. */
  virtual double PairEnergyRealForP(Bead&, Bead&, int) = 0;
  /** Compute reciprocalpair energy for scaled volume in pressure calculations.
      */
  virtual double PairEnergyReplForP(Bead&, Bead&, int) = 0;
  /** Compute real pair force along Z direction, for pressure calculations. */
  virtual double PairForceZReal(Bead&, Bead&, int) = 0;
  /** Compute reciprocal pair force along Z direction, for pressure
      calculations. */
  virtual double PairForceZRepl(Bead&, Bead&, int) = 0;
  /***/
  virtual double ForceZDipole(Bead&, double) = 0;
  virtual double PairDForceReal(Bead&, Bead&, Bead&, Bead&, int) = 0;
  virtual double PairDForceRepl(Bead&, Bead&, Bead&, Bead&, int) = 0;
  virtual double DipoleE(vector<Molecule>&) = 0;
  virtual double DipoleE(vector<Molecule>&, vector<Bead>&) = 0;
  virtual double DipoleE(vector<Molecule>&, int, int) = 0;
  virtual double DipoleEDiff(vector<Molecule>&, vector<Bead>&, Bead&, Bead&,
                             int, int, double, int) = 0;
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
  double GetEReal(int, int, int);
  double GetERepl(int, int, int);
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
  bool UseDipoleCorrection();
}; 

#endif


