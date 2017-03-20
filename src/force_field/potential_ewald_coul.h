#ifndef SRC_FORCE_FIELD_POTENTIAL_EWALD_COUL_H_
#define SRC_FORCE_FIELD_POTENTIAL_EWALD_COUL_H_

#include <iostream>
#include <map>
#include <string>

#include "../molecules/bead.h"
#include "potential_ewald.h"

using namespace std; 

class PotentialEwaldCoul : public PotentialEwald {
 private:
  /** The dielectric constant, unitless. */
  double dielectric;
  /** Bjerrum length. Derived from the dielectric constant. */
  double lB;
  double kBT_derived;
  double T_derived;
  /** The "radial" cutoff of the number of unit cells to do real space sum.
      num_x+num_y+num_z < real_cutoff.\n
      real_cutoff has to be an integrate equal or bigger than 1. */
  double real_cutoff;
  /** Maximum number of cells to encompass. */
  int real_cell[3];
  /** Cut off after this number of k in the reciprocal space. */
  double repl_cutoff;
  /** Maximum number of cells to encompass. */
  int repl_cell[3];
  /** Proportional to the inverse of the STD of the Gaussians, unit A^-1. */
  double alpha;

  // For efficiency of calculating the reciprocal component.
  /** The total number of k number according to repl_cell[3]. */
  int repl_ceto[3];
  /** k vector values. */
  double * kx;
  double * ky;
  double * kz;
  /** "k2" for all vector k. */
  double * k2;
  /** "exp(-k2/(4*alpha))/k2" for all vector k. */
  double * ek2;
  double * kz_forP;
  double * k2_forP;
  double * ek2_forP;

 public: 
  // Initialization functions.
  PotentialEwaldCoul(string, double[3]);
  ~PotentialEwaldCoul();
  /** Read parameters from file, can read in the sigmas and epsilons from
      multiple chain types (symbols). */
  void ReadParameters();

  // Energy functions.
  /** Energy between two beads. */
  double PairEnergyReal(Bead&, Bead&, int);   
  double PairEnergyRepl(Bead&, Bead&, int);
  double PairEnergyRealForP(Bead&, Bead&, int);
  double PairEnergyReplForP(Bead&, Bead&, int);
  double SelfEnergy(Bead&);
  double DipoleE(vector<Molecule>&);
  double DipoleE(vector<Molecule>&, vector<Bead>&);
  double DipoleE(vector<Molecule>&, int, int);
  double DipoleEDiff(vector<Molecule>&, vector<Bead>&, Bead&, Bead&, int, int,
                     double, int);
  /** The forces along Z direction. */
  double PairForceZReal(Bead&, Bead&, int);
  double PairForceZRepl(Bead&, Bead&, int);
  double ForceZDipole(Bead&, double);
  /** Vector D multiplies the force, see Yethiraj's papers. */
  double PairDForceReal(Bead&, Bead&, Bead&, Bead&, int);
  double PairDForceRepl(Bead&, Bead&, Bead&, Bead&, int);

  double GetlB();

}; 

#endif


