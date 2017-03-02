#ifndef SRC_MOLECULES_MOLECULE_H_
#define SRC_MOLECULES_MOLECULE_H_

#include <iostream>
#include <random>
#include <vector>

#include "bead.h"

/** A molecule contains a vector of monomer beads. For single particles, this
    vector just contains one bead. It also contains 2D arrays giving indices
    of beads involved in bonds/angles/dihedrals. The molecule is in charge of
    bead positions, so it has functions to do translational MC moves, scale 
    coordinates for pressure calculation, and print bead positions. */
class Molecule {
 private:
  int len; 
  int n_bond; 
  int n_angle; 
  int n_dihed; 

 public: 
  Molecule(); 
  Molecule(const Molecule&); 

  /** Storing the array of monomers/beads belonging to the molecule. */
  vector<Bead> bds;
  // Bonds are stored in a data structure of size nBonds x 2, where the two
  // entries are the indices of the two involved beads from the monomer vector.
  // for each mol, these include 0-nBead.
  vector<vector<int>> bonds;
  vector<vector<int>> angles; 
  vector<vector<int>> diheds; 

  // Add beads/topological features
  void AddBead(Bead);
  void AddBead(string, int, int, double, double, double, double); 
  void AddBond(int, int); 
  void AddAngle(int, int, int); 
  void AddDihed(int, int, int, int); 
  // Get info about molecule
  int Size(); 
  int NBond(); 
  int NAngle(); 
  int NDihed(); 

  ////////////////////////
  // Monte Carlo moves. //
  ////////////////////////
  /** Translate a random monomer or a single-bead molecule by a small random
      distance. */
  void BeadTranslate(double, double[3],  std::mt19937&); 
  /** Translate the entire molecule by a small random distance. */
  void COMTranslate(double delta, std::mt19937&);
  /** Pivot algorithm changed by Nuo from rotating around an existing bond to
      rotating around a point such that pivot itself is ergodic. [9/28/2016]\n
      Now it also can choose to rotate either the left or right part of the
      chain. [9/28/2016].\n
      The current pivot algorithm does not require rigid bond to use.  */
  void Pivot(double delta, std::mt19937&, double);
  /** Basically (Rachel's original) pivot for an internal section of the
      molecule. */
  void Crankshaft(double delta, std::mt19937&); 
  /** Unbiased reptation in both forward and backward direcitons. Currently,
      the routine has to be used with the "rigid bond" setting. */
  void RandomReptation(std::mt19937&, double);

}; 

#endif


