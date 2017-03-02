#ifndef SRC_MOLECULES_BEAD_H_
#define SRC_MOLECULES_BEAD_H_

#include <iomanip>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

/** Basic object describing a single particle contains position along with other
    attributes. */
class Bead {
 private:
  /** Current bead position. */
  double current_pos[3];
  /** Trial position after MC move. */
  double trial_pos[3];
  /** The charge of the bead. */
  double q;
  /** Unique identifier (not reused). */
  int id;
  /** The ID of the chain the bead belongs to. */
  int c_id;
  /** Bead symbol is used to decide what parameter to use for the bead. */
  string symbol;
  /** The type of this bead, ion, monomer etc. */
  int type;
  /** Record is the bead was moved. */
  bool moved; //did the bead move this time? 

 public: 
  /** Default constructor, should not be used. */
  Bead();
  /** Constructor that should be used. The arguments are symbol, id, cid,
      charge, x, y, z. */
  Bead(string, int, int, double, double, double, double);
  /** Copy constructor. */
  Bead(const Bead&);

  // Utilities.
  double Charge();
  void SetCharge(double);
  int ID();
  int ChainID();
  void SetID(int); 
  void SetChainID(int);
  string Symbol();
  void SetSymbol(string); 
  int Type();
  void SetType(int); 
  bool GetMoved();
  void SetMoved();
  void UnsetMoved();
  /** First arguement is a flag, 0 - current_pos, 1 - trial_pos.
      Second argument is index, 0 - x, 1 - y, 2 - z. */
  double GetCrd(int, int);
  /** First two arguments are the same as GetCrd, last one is the value
      to set. */ 
  void SetCrd(int, int, double);
  void SetAllCrd(double[3]);  //set xyz using a coord array
  /** Copy trial_pos to current_pos. */
  void UpdateCurrentPos();
  /** Copy current_pos to trial_pos. */
  void UpdateTrialPos(); 

  // Distance calculations.
  /** Calculate the minimum distance between two beads considering all of the
      periodic images. This function uses the trial coordinates! */
  double BBDist(Bead&, double[], int);
  double BBDistC(Bead&, double[], int);
  double BBDistVec(Bead&, double[], int, int);
  double BBDistVecWithRef(Bead&, Bead&, Bead&, double[], int, int);
  double BBDistScale(Bead&, double[], int, int, double);
  /** The "naive distance" directly returns the distance between two beads
      without shifting due to periodic boundary conditions. This is used to
      calculate Rg for a connected polymer that shouldn't be folded back into
      into the same unit cell. In fact, if the "fold back" case ever happens, a
      larger unit cell should be used. These functions use the trial
      coordinates! */
  double NaiveBBDistX(Bead&);  // Nuo added: 8/9/2016.
  double NaiveBBDistY(Bead&);  // Nuo added: 8/9/2016.
  double NaiveBBDistZ(Bead&);  // Nuo added: 8/9/2016.
  /** Returns a string that stores the distance to the z=0 and the distance to
      the z=box_l[2] wall. */
  string DistToWall(double[3]);      // Nuo added: 9/7/2016.

};
 
#endif


