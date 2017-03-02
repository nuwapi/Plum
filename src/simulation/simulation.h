/** The Simulation object contains the molecule vector as well as info about how
    many things are in the box and general MC params related to moving mols
    (beta, delta, etc).
*/

#ifndef SRC_SIMULATION_SIMULATION_H_
#define SRC_SIMULATION_SIMULATION_H_ 

#include <fstream>
#include <random> 
#include <set>
#include <vector>

#include "../force_field/force_field.h"
#include "../molecules/molecule.h"
#include "../utilities/constants.h"

using namespace std; 

class Simulation {
 private:
  //////////////////////////////////////
  // Parameters read from input file. //
  //////////////////////////////////////
  /** Input coordinate file. */
  string crd_name;
  /** Input topology file. */
  string top_name;
  /** Name of the simulation and the prefix for the output file names. */
  string run_name;
  /** The current step number. */
  long step;
  /** Total number of steps. */
  long steps;
  /** Number of steps of equilibration, which are not counted towards sampling.
    */
  long steps_eq;
  /** The frequency to calculate energy, force, density, Rg etc. */
  long sample_freq;
  /** The frequency to print energy, force etc. */
  long stat_out_freq;
  /** The frequency to print trajectory file. */
  long traj_out_freq;
  /** The number of dimensions to apply PBC, use 2 or 3. */
  int npbc;
  /** The length of displacement for random moves, in unit length. */
  double move_size;
  /** 1/kBT, should be set to 1. The energy unit of the simulation. */
  double beta;
  /** The number of phantom beads, each phantom bead is considered as one
      molecule, they don't have a size and they never move. Phantom beads
      should be used together with appropriate parameters in the input file. */
  int phantom;
  /** Coions are there to neutralize the surface charges represented by the
      phantom molecules. */
  int coion;
  /** Number of grafted chains, i.e. chains with monomer(s) that have symbol
      "L" (left, z=0) or "R" (right, z=box_l[2]). */
  int grafted;
  /** The counterions need for the grafted chains, the counterions themselves
      are not grafted. */
  int grafted_counterion;
  ifstream crd_in;
  ifstream top_in;
  /** Contains info on energy, pressure, Rg etc. */
  ofstream info_out;
  /** For polymer trajectories. */
  ofstream traj_p_out;
  /** For cation trajectories. */
  ofstream traj_c_out;
  /** For anion trajectories. */
  ofstream traj_a_out;
  /** For "neutral ion" trajectories. */
  ofstream traj_n_out;
  ofstream crd_last_out[2];
  ofstream top_last_out;
  ofstream rho_last_out[2];

  /////////////////////////////////////////////
  // Parameters read from top and crd files. //
  /////////////////////////////////////////////
  /** Total number of Beads. */
  int n_bead;
  /** Total number of Molecules. */
  int n_mol;
  /** Total number of chain Molecules. */
  int n_chain;
  /** Total number of beads in all chains. */
  int n_chain_b;
  /** Total number cations. */
  int n_cion;
  /** Total number fo anions. */
  int n_aion;
  /** Total number of "neutral ions". */
  int n_nion;
  /** The x,y,z length of the simulation box, unit length. */
  double box_l[3];
  /** The probabilities to choose one of the four MC translational moves. */
  double move_prob[kNoMoveType];
  /** Molecule vector. */
  vector<Molecule> mols;
  /** Force field for the simulation. */
  ForceField force_field;

  //////////////////////////////////////////
  // Variables for simulation statistics. //
  //////////////////////////////////////////
  bool calc_chem_pot;
  double chem_pot_cumu;
  int accepted[kNoMoveType];
  int attempted[kNoMoveType];
  int insertion_accepted;
  int deletion_accepted;
  int insertion_attempted;
  int deletion_attempted;
  int ff_avg_counter;
  int mol_avg_counter;
  double pair_e_cumu;
  double ewald_e_cumu;
  double ewald_e_real_cumu;
  double ewald_e_repl_cumu;
  double ewald_e_self_cumu;
  double bond_e_cumu;
  double ext_e_cumu;
  /** Running total density value. */
  double density_cumu;
  int density_z_bin;
  double density_z_res;
  double * density_z_cumu;
  double rg_tot_cumu;
  double rg_x_cumu;
  double rg_y_cumu;
  double rg_z_cumu;
  /** End to end distance. */
  double e_to_e_cumu;
  /** Chains adsorbed per unit area, 1/A2. */
  int adsorbed_chains;
  /** Beads adsorbed per unit area, 1/A2. */
  int adsorbed_beads;
  /** Adsportion percentage, %. */
  double adsorption_percent;
  int id_counter;
  set<int> id_list;

  //////////////////////////////
  // Random number generator. //
  //////////////////////////////
  mt19937 rand_gen;

 public:
  /** Constructor. */
  Simulation();
  /** Destructor. */
  ~Simulation();

  ///////////////////////////
  // Simulation functions. //
  ///////////////////////////
  /** Run this simulation. */
  void Run();
  /** Attempts to perform one of the available translational moves on a random
      molecule; Accepts/rejects and adjusts all coordinate & energy arrays. */
  void TranslationalMove(); 
  /** Attempts a grand-canonical Monte Carlo move. */
  void GCMove();

  ////////////////
  // Utilities. //
  ////////////////
  /** Sample energy, pressure, molecular conformation etc. */
  void Sample();
  void CalcRhoZ();
  /** Calculates average value of sqrt(Rg^2) for all chains. */
  double RadiusOfGyration();
  string RadiusOfGyrationXYZ();
  double EndToEndDistance();
  void DetectAdsorption();
  /** Always generate new ID, never reuse. Warning: could run out of id due to
      int out of bound error. */
  int GenBeadID();
  void UpdateMolCounts();

  ///////////////////////////////////
  // Read input crd and top files. //
  ///////////////////////////////////
  void ReadCrd(); 
  void ReadTop();
  void CoordinateObeyRigidBond(double);

  //////////////////////
  // Print functions. //
  //////////////////////
  void PrintStatHeader();
  void PrintStat();
  void PrintTraj();
  void PrintLastCrd();
  void PrintLastTop();
  void PrintLastRhoZ();
  void PrintEndToEndVector();

};
 
#endif


