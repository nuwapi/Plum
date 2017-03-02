//A simple bonded potential with only harmonic bonds between specified atoms
//no angle or dihedral potentials 
//currently the only PotentialBond supported for gc simulations 

#ifndef SPRING_H
#define SPRING_H

#include <stdlib.h> 
#include <iomanip> 
#include <iostream> 
#include <random> 

#include "../molecules/bead.h"
#include "../molecules/molecule.h"
#include "../utilities/misc.h"
#include "potential_bond.h"

class PotentialSpring: public PotentialBond {
 private:
  double m_kBond;
  double m_r0;

 public:
  PotentialSpring(int, string); 
  void ReadParameters();
  double MoleculeEnergy(Molecule&, double[], int);
  double RandomBondLen(double, mt19937&); 
  double EnergyDifference(vector < Molecule >&, double[], int, int);

};

#endif


