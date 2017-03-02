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
  int steps_eq;
  /** The frequency to calculate energy, force, density, Rg etc. */
  int sample_freq;
  /** The frequency to print energy, force etc. */
  int stat_out_freq;
  /** The frequency to print trajectory file. */
  int traj_out_freq;
  /** The number of dimensions to apply PBC, use 2 or 3. */
  int npbc;
  /** The length of displacement for random moves, in unit length. */
  double move_size;
  /** 1/kBT, should be set to 1. The energy unit of the simulation. */
  double beta;
  /** The number of phantom beads, they don't have a size and they never
      move. Phantom beads should be used together with appropriate parameters
      in the input file. */
  int phantom;
  ifstream crd_in;
  ifstream top_in;
  /** Contains info on energy, pressure, Rg etc. */
  ofstream info_out;
  ofstream traj_out;
  ofstream crd_last_out;
  ofstream top_last_out;

  /////////////////////////////////////////////
  // Parameters read from top and crd files. //
  /////////////////////////////////////////////
  /** Total number of Beads. */
  int n_bead;
  /** Total number of Molecules. */
  int n_mol;
  /** Total number of chain Molecules. */
  int n_chain;
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
  /** Calculates average value of sqrt(Rg^2) for all chains. */
  double RadiusOfGyration();
  string RadiusOfGyrationXYZ();
  double EndToEndDistance();
  void DetectAdsorption();
  /** Always generate new ID, never reuse. Warning: could run out of id due to
      int out of bound error. */
  int GenBeadID();

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
  void PrintEndToEndVector();

};
 
#endif


