#include "molecule.h"

#include <stdlib.h>
#include <cmath>
#include <iostream> 
#include <iomanip>
#include <random>
#include <vector>

#include <Eigen/Dense>
#include "bead.h"
#include "../utilities/constants.h"
#include "../utilities/misc.h"

using namespace std; 
using namespace Eigen;

Molecule::Molecule() {
  len = 0;
  n_bond = 0;
  n_angle = 0;
  n_dihed = 0;

}

// Copy constructor.
Molecule::Molecule(const Molecule& mol) {
  len = 0;  // These values will be incremented 
  n_bond = 0; // as the add functions are called. 
  n_angle = 0;
        n_dihed = 0;
  for (int i = 0; i < (int)mol.bonds.size(); i++) {
    AddBond(mol.bonds[i][0], mol.bonds[i][1]);
  }
  for (int i = 0; i < (int)mol.angles.size(); i++) {
    AddAngle(mol.angles[i][0], mol.angles[i][1], mol.angles[i][2]); 
  }
  for (int i = 0; i < (int)mol.diheds.size(); i++) {
    AddDihed(mol.diheds[i][0], mol.diheds[i][1], mol.diheds[i][2], mol.diheds[i][3]); 
  }
  for (int i = 0; i < (int)mol.bds.size(); i++) {
    AddBead(mol.bds[i]);
  }

}

// Adds a bead to the first empty slot in the bead array.
void Molecule::AddBead(Bead b) {
  bds.push_back(b); 
  len++; 

}

void Molecule::AddBead(string symbol, int id, int c_id, double charge, double x, double y,
                       double z) {
  bds.push_back(Bead(symbol, id, c_id, charge, x, y, z)); 
  bds[len].SetCrd(0, 0, x);
  bds[len].SetCrd(0, 1, y);
  bds[len].SetCrd(0, 2, z);
  len++;

}

void Molecule::AddBond(int ind1, int ind2){
  vector<int>bond = {ind1, ind2}; 
  bonds.push_back(bond); 
  n_bond++; 
}

int Molecule::NBond() {
  return n_bond; 

}

void Molecule::AddAngle(int ind1, int ind2, int ind3){
  vector < int > angle = {ind1, ind2, ind3}; 
  angles.push_back(angle); 
  n_angle++; 
}

void Molecule::AddDihed(int ind1, int ind2, int ind3, int ind4){
  vector < int > di = {ind1, ind2, ind3, ind4};
  diheds.push_back(di); 
  n_angle++;  
}

//returns number of beads in molecule
int Molecule::Size() {
  return len; 

}

int Molecule::NAngle() {
  return n_angle;

}

int Molecule::NDihed() {
  return n_dihed;

}

void Molecule::BeadTranslate(double move_size, double box_l[3],
                             mt19937& rand_gen) {
  double vec[3];

  // Routine for translating beads.
  randSphere(vec, rand_gen);
  double vec_len = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  if (vec_len > 0) {
    // 3 times the move size of chains.
    vec_len = 3*move_size / vec_len;
  }

  bds[0].SetMoved();

  for (int i = 0; i < 3; i++) {
    double old = bds[0].GetCrd(0, i);
    bds[0].SetCrd(1, i, old + vec_len * vec[i]);
  }

  /* 
  // Routine for dropping beads randomly in the box.
  vec[0] = (double)rand_gen()/rand_gen.max() * box_l[0];
  vec[1] = (double)rand_gen()/rand_gen.max() * box_l[1];
  vec[2] = (double)rand_gen()/rand_gen.max() * box_l[2];
  bds[0].SetMoved();
  for (int i = 0; i < 3; i++) {
    bds[0].SetCrd(1, i, vec[i]);
  }
  */

}

void Molecule::COMTranslate(double move_size, mt19937& rand_gen) {
  double vec[3] = {0.5*move_size * (double)rand_gen()/rand_gen.max(),
                   0.5*move_size * (double)rand_gen()/rand_gen.max(),
                   0.5*move_size * (double)rand_gen()/rand_gen.max()};
  if (rand_gen()%2 == 0)  vec[0] = -vec[0];
  if (rand_gen()%2 == 0)  vec[1] = -vec[1];
  if (rand_gen()%2 == 0)  vec[2] = -vec[2];

  for (int i = 0; i < len; i++) { 
    for (int j = 0; j < 3; j++) {
      double old = bds[i].GetCrd(0, j);
      bds[i].SetCrd(1, j, old + vec[j]); 
    }
    bds[i].SetMoved();
  }

}

void Molecule::Pivot(double move_size, mt19937& rand_gen, double rigid_bond) {
  // Choose the pivot bead.
  int pivot = floor(len * (double)rand_gen()/rand_gen.max());
  if (pivot == len)  pivot--;
  // This if for grafted polymers.
  if (bds[0].Symbol() == "R" || bds[0].Symbol() == "L")
    pivot = 0;

  // Determine the move size.
  double move_size_rand = move_size * (double)rand_gen()/rand_gen.max();

  // Forward rotation.
  for (int i = pivot+1; i < len; i++) {
    double vec[3];
    randSphere(vec, rand_gen);
    double x = bds[i].GetCrd(1, 0) + move_size_rand*vec[0];
    double y = bds[i].GetCrd(1, 1) + move_size_rand*vec[1];
    double z = bds[i].GetCrd(1, 2) + move_size_rand*vec[2];
    double distx = x - bds[i-1].GetCrd(1, 0);
    double disty = y - bds[i-1].GetCrd(1, 1);
    double distz = z - bds[i-1].GetCrd(1, 2);
    double norm = rigid_bond / sqrt(distx*distx + disty*disty + distz*distz);
    // This is the new center.
    double cx = norm*distx + bds[i-1].GetCrd(1, 0);
    double cy = norm*disty + bds[i-1].GetCrd(1, 1);
    double cz = norm*distz + bds[i-1].GetCrd(1, 2);
    // This is amount to translate for the rest of the chain.
    double mx = cx - bds[i].GetCrd(1, 0);
    double my = cy - bds[i].GetCrd(1, 1);
    double mz = cz - bds[i].GetCrd(1, 2);
    // Translate.
    for (int j = i; j < len; j++) {
      bds[j].SetCrd(1, 0, bds[j].GetCrd(1, 0)+mx);
      bds[j].SetCrd(1, 1, bds[j].GetCrd(1, 1)+my);
      bds[j].SetCrd(1, 2, bds[j].GetCrd(1, 2)+mz);
    }
  }

  // Backward rotation.
  for (int i = pivot-1; i >= 0; i--) {
    double vec[3];
    randSphere(vec, rand_gen);
    double x = bds[i].GetCrd(1, 0) + move_size_rand*vec[0];
    double y = bds[i].GetCrd(1, 1) + move_size_rand*vec[1];
    double z = bds[i].GetCrd(1, 2) + move_size_rand*vec[2];
    double distx = x - bds[i+1].GetCrd(1, 0);
    double disty = y - bds[i+1].GetCrd(1, 1);
    double distz = z - bds[i+1].GetCrd(1, 2);
    double norm = rigid_bond / sqrt(distx*distx + disty*disty + distz*distz);
    // This is the new center.
    double cx = norm*distx + bds[i+1].GetCrd(1, 0);
    double cy = norm*disty + bds[i+1].GetCrd(1, 1);
    double cz = norm*distz + bds[i+1].GetCrd(1, 2);
    // This is amount to translate for the rest of the chain.
    double mx = cx - bds[i].GetCrd(1, 0);
    double my = cy - bds[i].GetCrd(1, 1);
    double mz = cz - bds[i].GetCrd(1, 2);
    // Translate.
    for (int j = i; j >= 0; j--) {
      bds[j].SetCrd(1, 0, bds[j].GetCrd(1, 0)+mx);
      bds[j].SetCrd(1, 1, bds[j].GetCrd(1, 1)+my);
      bds[j].SetCrd(1, 2, bds[j].GetCrd(1, 2)+mz);
    }
  }

  // Finally, set moved.
  for (int i = 0; i < len; i++) {
    if (i != pivot)  bds[i].SetMoved();
  }

}

void Molecule::Crankshaft(double move_size, mt19937& rand_gen) {
  int first = floor((len-2) * (double)rand_gen()/rand_gen.max()); 
  int last = floor((len - first) * (double)rand_gen()/rand_gen.max()) + first; 
  Vector3d axis; 
  for (int i = 0; i < 3; i++) {
    axis(i) = bds[last].GetCrd(0,i) - bds[first].GetCrd(0,i); 
  }
  axis.normalize();
  double angle = move_size * ((double)rand_gen()/rand_gen.max() * 2 * M_PI - M_PI);
  AngleAxisd a = AngleAxisd(angle, axis); 
  Matrix3d rot;
  rot = a.toRotationMatrix(); 
  for (int i = first + 1; i < last; i++) {
    // Toggle moved flags for affected beads to true.
    bds[i].SetMoved();
    Vector3d v; // Put coords into eigen vector object, make the first bead the origin.
    for (int k = 0; k < 3; k++) {
      v(k) = bds[i].GetCrd(0,k) - bds[first].GetCrd(0,k);
    }
    // Rotate and put them back, restoring normal origin.
    v = rot * v;
    for (int k = 0; k < 3; k++) {
      bds[i].SetCrd(1, k, v(k) + bds[first].GetCrd(0,k));
    }
  }
  
}

void Molecule::RandomReptation(mt19937& rand_gen, double rigid_bond) {
  if (rigid_bond < 0) {
    cout << "Rigid bond is not properly defined. The current reptation move "
         << "only supports freely joined hard sphere chains. "
         << "Exiting! Program complete." << endl;
    exit(1);
  }

  // Forward direction.
  int direction = 1;
  int begin = 0;
  int end = len - 1;
  if (rand_gen() % 2 == 0) {
    // Backward direction.
    direction = -1;
    begin = len - 1;
    end = 0;
  }

  for (int i = begin; i != end; i+=direction) {
    for (int j = 0; j < 3; j++) {
      bds[i].SetCrd(1, j, bds[i+direction].GetCrd(0, j));
    }
  }

  double vec[3];
  randSphere(vec, rand_gen);
  for (int i = 0; i < 3; i++) {
    vec[i] *= (rigid_bond);
    bds[end].SetCrd(1, i, bds[end].GetCrd(0, i)+vec[i]);
  }


  for (int i = 0; i < len; i++) {
    bds[i].SetMoved();
  }

}


