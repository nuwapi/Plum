#include "bead.h"
#include <cmath> 
#include <iomanip> 
#include <string> 
#include <iostream> 
#include <fstream> 
#include <sstream>

using namespace std; 

Bead::Bead() {
  for(int i = 0; i < 3; i++) {
    trial_pos[i] = 0; 
    current_pos[i] = 0; 
  }
  q = 0; 
  symbol = "NONE";
  moved = false;

}

Bead::Bead(string symbol_in, int id_in, int c_id_in, double q_in, double x,
           double y, double z) {
  trial_pos[0] = x; 
  trial_pos[1] = y; 
  trial_pos[2] = z; 
  for(int i = 0; i < 3; i++) {
    current_pos[i] = trial_pos[i]; 
  }
  symbol = symbol_in;
  q = q_in;
  moved = false;  
  id = id_in;
  c_id = c_id_in;

}

Bead::Bead(const Bead& bead) {
  for(int i = 0; i < 3; i++) {
    trial_pos[i] = bead.trial_pos[i]; 
  }
  for(int i = 0; i < 3; i++){
    current_pos[i] = bead.current_pos[i]; 
  }  
  q = bead.q; 
  type = bead.type; 
  symbol = bead.symbol;
  moved = bead.moved;  
  id = bead.id;

}

string Bead::Symbol() {
  return symbol; 

}

int Bead::Type() {
  return type; 

}

double Bead::Charge() {
  return q;

}

int Bead::ID() {
  return id; 

}

int Bead::ChainID() {
  return c_id;

}

double Bead::GetCrd(int flag, int index) {
  if (flag == 0) {
    return current_pos[index]; 
  }
  else if (flag == 1) {
    return trial_pos[index]; 
  }
  else {
    cout << "Bead::GetCrd:\n  Invalid flag!" << endl;
    exit(1);
  }
  return 0; 

}

void Bead::SetCrd(int flag, int index, double value) {
  if (flag == 0) {
    current_pos[index] = value; 
  }
  else if (flag == 1) {
    trial_pos[index] = value; 
  }
  else {
    cout << "Bead::SetCrd:\n  Invalid flag!" << endl;
  }

}

void Bead::SetAllCrd(double pos[]) {
  for(int i = 0; i < 3; i++) {
    current_pos[i] = pos[i]; 
    trial_pos[i] = pos[i]; 
  }

}


void Bead::SetMoved() {
  moved = true;
}

void Bead::UnsetMoved() {
  moved = false;
}

bool Bead::GetMoved() {
  return moved; 
}

void Bead::SetType(int type_in) {
  type = type_in; 
}

void Bead::SetID(int id_in) {
  id = id_in;
}

void Bead::SetChainID(int c_id_in) {
  c_id = c_id_in;

}

void Bead::SetSymbol(string symbol_in) {
  symbol = symbol_in; 
}

void Bead::SetCharge(double q_in) {
  q = q_in;
}

void Bead::UpdateTrialPos() {
  for (int i = 0; i < 3; i++) {
    trial_pos[i] = current_pos[i]; 
  }

}

void Bead::UpdateCurrentPos() {
  for (int i = 0; i < 3; i++) {
    current_pos[i] = trial_pos[i]; 
  }

}

double Bead::BBDist(Bead& bead, double box_l[], int npbc) {
  double dist = 0;
  for (int i = 0; i < 3; i++) {
    double dist_i = trial_pos[i] - bead.trial_pos[i]; 
    if (i < npbc) {
      dist_i -= box_l[i] * round(dist_i / box_l[i]);
      if (abs(dist_i) > 0.5*box_l[i]) {
        dist_i = box_l[i] - abs(dist_i);
      }
    }
    dist += (dist_i * dist_i); 
  }
  return sqrt(dist); 

}

double Bead::BBDistC(Bead& bead, double box_l[], int npbc) {
  double dist = 0;
  for (int i = 0; i < 3; i++) {
    double dist_i = current_pos[i] - bead.current_pos[i];
    if (i < npbc) {
      dist_i -= box_l[i] * round(dist_i / box_l[i]);
      if (abs(dist_i) > 0.5*box_l[i]) {
        dist_i = box_l[i] - abs(dist_i);
      }
    }
    dist += (dist_i * dist_i);
  }
  return sqrt(dist);

}

double Bead::BBDistVec(Bead& bead, double box_l[], int npbc, int index) {
  double dist_i = trial_pos[index] - bead.trial_pos[index];
  if (index < npbc) {
    dist_i -= box_l[index] * round(dist_i / box_l[index]);
    /*
    if (dist_i < -0.5*box_l[index]) {
      dist_i += box_l[index];
    }
    else if (dist_i > 0.5*box_l[index]) {
      dist_i -= box_l[index];
    }
    */
  }

  return dist_i;

}

double Bead::BBDistVecWithRef(Bead& bead, Bead& ref1, Bead& ref2,
                              double box_l[], int npbc, int index) {
  double dist_i = trial_pos[index] - bead.trial_pos[index];
  double ref_i = ref1.trial_pos[index] - ref2.trial_pos[index];
  if (index < npbc) {
    dist_i -= box_l[index] * round(ref_i / box_l[index]);
  }

  return dist_i;

}

double Bead::BBDistScale(Bead& bead, double box_l[], int npbc, int index, double p_scaling) {
  double dist = 0;
  for (int i = 0; i < 3; i++) {
    double dist_i = trial_pos[i] - bead.trial_pos[i];
    if (i < npbc) {
      dist_i -= box_l[i] * round(dist_i / box_l[i]);
      if (abs(dist_i) > 0.5*box_l[i]) {
        dist_i = box_l[i] - abs(dist_i);
      }
    }
    if (i == index) {
      dist_i -= dist_i*p_scaling;
    }
    dist += (dist_i * dist_i);
  }
  return sqrt(dist);

}

double Bead::NaiveBBDistX(Bead& bead) {
  return abs(trial_pos[0] - bead.trial_pos[0]); 

}

double Bead::NaiveBBDistY(Bead& bead) {
  return abs(trial_pos[1] - bead.trial_pos[1]);

}

double Bead::NaiveBBDistZ(Bead& bead) {
  return abs(trial_pos[2] - bead.trial_pos[2]);

}

string Bead::DistToWall(double box_l[]) {
  double z_coord = trial_pos[2] - box_l[2] * floor(trial_pos[2]/box_l[2]);
  std::ostringstream foo;
  foo << z_coord << " " << box_l[2] - z_coord;
  return foo.str();

}


