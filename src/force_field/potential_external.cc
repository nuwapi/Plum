#include "potential_external.h"

#include "../utilities/constants.h"

using namespace std;

// Default constructor.
PotentialExternal::PotentialExternal(string potential_name) {
  name = potential_name;
  dE = 0;
  E_tot = 0;

}

void PotentialExternal::SetE(int flag, int key, double val) {
  if (flag == 0) {
    current_energy_map[key] = val; 
  }
  else if (flag == 1) {
    trial_energy_map[key] = val; 
  }
  else {
    cout << "  Invalid map requested!" << endl;
  }

}

void PotentialExternal::SetEBothMaps(int key, double val) {
  current_energy_map[key] = val;
  trial_energy_map[key] = val;

}

double PotentialExternal::GetE(int flag, int key) {
  if (flag == 0) {
    return current_energy_map[key]; 
  }
  else if (flag == 1) {
    return trial_energy_map[key]; 
  }
  else {
    cout << "  Invalid map requested!" << endl; 
  }
  return 1;

}

void PotentialExternal::EnergyInitialization(vector<Molecule>& mols,
                                             double box_l[],
                                             int npbc) {
  E_tot = 0;

  for (int i = 0; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      double en = BeadEnergy(mols[i].bds[j], box_l); 
      E_tot += en; 
      SetEBothMaps(mols[i].bds[j].ID(), en); 
    }
  }

}

void PotentialExternal::EnergyInitForLastMol(vector<Molecule>& mols, int chain_len,
                                             double bead_charge, double box_l[],
                                             int npbc) {
  int added = 1;
  // Assume all chains are the same and always go before their counterions!!!
  // Assume monovalent ions.
  if (bead_charge != 0)
    added += chain_len;

  for (int i = (int)mols.size()-added; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      double en = BeadEnergy(mols[i].bds[j], box_l);
      E_tot += en;
      SetEBothMaps(mols[i].bds[j].ID(), en);
    }
  }

}

void PotentialExternal::AdjustEnergyUponMolDeletion(vector<Molecule>& mols,
                                                    int delete_id) {
  // Assume monovalent counterion!
  int counterion = 0;
  for (int i = 0; i < mols[delete_id].Size(); i++) {
    counterion += (int)round(abs(mols[delete_id].bds[i].Charge()));
  }

  for (int i = delete_id; i <= delete_id+counterion; i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      E_tot -= current_energy_map[mols[i].bds[j].ID()]; 
      current_energy_map.erase(mols[i].bds[j].ID()); 
      trial_energy_map.erase(mols[i].bds[j].ID());
    }
  }

}

double PotentialExternal::EnergyDifference(vector<Molecule>& mols,
                                           int active_mol, double box_l[],
                                           int npbc) {
  dE = 0; 
  for (int i = 0; i < (int)mols[active_mol].Size(); i++) {
    if (mols[active_mol].bds[i].GetMoved()) {
      double eNew = BeadEnergy(mols[active_mol].bds[i], box_l); 
      if (eNew >= kVeryLargeEnergy) {
        dE = kVeryLargeEnergy;
        return dE;
      }
      SetE(1, mols[active_mol].bds[i].ID(), eNew); 
      dE += eNew - GetE(0, mols[active_mol].bds[i].ID()); 
    }
  }
  return dE; 

}

void PotentialExternal::FinalizeEnergyBothMaps(vector<Molecule>& mols,
                                               int active_mol, bool accepted) {
  if (accepted)  E_tot += dE; 
  for (int i = 0; i < (int) mols[active_mol].Size(); i++) {
    if (mols[active_mol].bds[i].GetMoved()) {
      if (accepted) {
        SetE(0, mols[active_mol].bds[i].ID(),
             GetE(1, mols[active_mol].bds[i].ID())); 
      }
      else {
        SetE(1, mols[active_mol].bds[i].ID(),
             GetE(0, mols[active_mol].bds[i].ID()));
      }
    }
  }

}

double PotentialExternal::GetTotalEnergy() {
  return E_tot; 

}

double PotentialExternal::CalcTrialTotalEnergy(vector < Molecule >& mols,
                                               double box_l[], int npbc) {
  double total_energy = 0;
  for (int i = 0; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      total_energy += BeadEnergy(mols[i].bds[j], box_l);
    }
  }

  return total_energy;

}

string PotentialExternal::PotentialName() {
  return name;

}


