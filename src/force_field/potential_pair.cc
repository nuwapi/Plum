#include "potential_pair.h"

#include <cmath>
#include <vector> 

#include "../utilities/constants.h"

using namespace std; 

// Default constructor.
PotentialPair::PotentialPair(string potential_name) {
  name = potential_name;
  dE = 0;
  E_tot = 0;

}

// Sets a pairwise energy value in either array.
void PotentialPair::SetE(int flag, int key1, int key2, double val) {
  if (flag == 0) {
    current_energy_map[make_pair(min(key1, key2), max(key1, key2))] = val; 
  }
  else if (flag == 1) {
    trial_energy_map[make_pair(min(key1, key2), max(key1, key2))] = val;
  }
  else {
    cout << "PotentialPair::SetE:\n  Invalid map requested!" << endl;
  }

}

// Puts the value into all pair slots in both arrays.
void PotentialPair::SetEBothMaps(int key1, int key2, double val) {
  current_energy_map[make_pair(min(key1, key2), max(key1, key2))] = val;
  trial_energy_map[make_pair(min(key1, key2), max(key1, key2))] = val;

}

// Gets a pairwise energy value from either array.
double PotentialPair::GetE(int flag, int key1, int key2) {
  if (flag == 0) {
    return current_energy_map[make_pair(min(key1, key2), max(key1, key2))]; 
  }
  else if (flag == 1) {
    return trial_energy_map[make_pair(min(key1, key2), max(key1, key2))];
  }
  else {
    cout << "PotentialPair::GetE:\n  Invalid flag." << endl;
    return 1; 
  }
  return 0;
 
}

void PotentialPair::EnergyInitialization(vector<Molecule>& mols, double box_l[],
                                         int npbc) {
  E_tot = 0;

  // First put in intramolecular pair energies.
  for (int i = 0; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size()-1; j++) {
      for (int k = j+1; k < mols[i].Size(); k++) {
        double energy;
        if (name == "HardSphere" && k == j+1) {
          //energy = PairEnergy(mols[i].bds[j], mols[i].bds[k], box_l, npbc);
          energy = 0;
        }
        else {
          energy = PairEnergy(mols[i].bds[j], mols[i].bds[k], box_l, npbc); 
        }
        E_tot += energy;
        SetEBothMaps(mols[i].bds[j].ID(), mols[i].bds[k].ID(), energy); 
        if (energy >= kVeryLargeEnergy) {
          cout << "  Overlap: "
               << mols[i].bds[j].BBDist(mols[i].bds[k], box_l, npbc) << ", "
               << i << " " << j << " | " << i << " " << k << endl;
        }
      }
    }
  }

  // Now intermolecular pair energies.
  for (int i = 0; i < (int)mols.size()-1; i++) {
    for (int j = i+1; j < (int)mols.size(); j++) {
      for (int k = 0; k < mols[i].Size(); k++) {
        for (int l = 0; l < mols[j].Size(); l++) {
          double energy;
          energy = PairEnergy(mols[i].bds[k], mols[j].bds[l], box_l, npbc);
          E_tot += energy;
          SetEBothMaps(mols[i].bds[k].ID(), mols[j].bds[l].ID(), energy); 
          if (energy >= kVeryLargeEnergy) {
            cout << "  Overlap: " << i << " " << k
                 << " | " << j << " " << l << endl;
          }
        }
      }
    }
  }

}

void PotentialPair::EnergyInitForLastMol(vector<Molecule>& mols, int chain_len,
                                         double bead_charge, double box_l[],
                                         int npbc) {
  double energy = 0;
  int added = 1;
  // Assume all chains are the same and always go before their counterions!!!
  // Assume monovalent ions.
  if (bead_charge != 0)
    added += chain_len;

  for (int i = (int)mols.size()-added; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      // With existing beads.
      for (int k = 0; k < (int)mols.size()-added; k++) {
        for (int l = 0; l < mols[k].Size(); l++) {
          int id1 = min(mols[k].bds[l].ID(), mols[i].bds[j].ID());
          int id2 = max(mols[k].bds[l].ID(), mols[i].bds[j].ID());
          if (name == "HardSphere") {
            energy = 0;
          }
          else {
            energy = PairEnergy(mols[i].bds[j], mols[k].bds[l], box_l, npbc);
          }
          E_tot += energy;
          SetEBothMaps(id1, id2, energy);
        }
      }

      // Within added beads.
      for (int k = (int)mols.size()-added; k <= i; k++) {
        for (int l = 0; l < mols[k].Size(); l++) {
          if (k < i || (k == i && l < j)) {
            int id1 = min(mols[k].bds[l].ID(), mols[i].bds[j].ID());
            int id2 = max(mols[k].bds[l].ID(), mols[i].bds[j].ID());
            if (name == "HardSphere") {
              energy = 0;
            }
            else {
              energy = PairEnergy(mols[i].bds[j], mols[k].bds[l], box_l, npbc);
            }
            E_tot += energy;
            SetEBothMaps(id1, id2, energy);
          }
        }
      }
    }
  }

}

double PotentialPair::EnergyDifference(vector<Molecule>& mols, int moved_mol,
                                       double box_l[], int npbc) {
  dE = 0;

  // For all pairs in the molecule that moved.
  for (int i = 0; i < mols[moved_mol].Size()-1; i++) {
    int id1 = mols[moved_mol].bds[i].ID();
    for (int j = i+1; j < mols[moved_mol].Size(); j++) {
      // Only when either i or j is moved.
      if (mols[moved_mol].bds[i].GetMoved() ||
          mols[moved_mol].bds[j].GetMoved()) {
        double new_e;
        int id2 = mols[moved_mol].bds[j].ID();
        if (name == "HardSphere" && j == i+1) {
          new_e = 0;
          //new_e = PairEnergy(mols[moved_mol].bds[i], mols[moved_mol].bds[j],
          //                   box_l, npbc);
        }
        else {
          new_e = PairEnergy(mols[moved_mol].bds[i], mols[moved_mol].bds[j],
                             box_l, npbc);
        }
        SetE(1, id1, id2, new_e);
        dE += (new_e - GetE(0, id1, id2));
      }
    }
  }

  // For all pairs with other molecules.
  for (int i = 0; i < (int)mols.size(); i++) {
    if (i != moved_mol) {
      for (int j = 0; j < mols[i].Size(); j++) {
        for (int k = 0; k < mols[moved_mol].Size(); k++) {
          // Only do it for moved beads.
          if (mols[moved_mol].bds[k].GetMoved()) {
            int id1 = min(mols[moved_mol].bds[k].ID(), mols[i].bds[j].ID());
            int id2 = max(mols[moved_mol].bds[k].ID(), mols[i].bds[j].ID());
            double new_e = PairEnergy(mols[moved_mol].bds[k], mols[i].bds[j],
                                      box_l, npbc);
            SetE(1, id1, id2, new_e);
            dE += (new_e - GetE(0, id1, id2));
          }
        }
      }
    }
  }

  return dE; 

}

void PotentialPair::FinalizeEnergyBothMaps(vector<Molecule>& mols,
                                           int moved_mol, bool accept) {
  // Adjust total energy variable.
  if (accept) {
     E_tot += dE;
  }

  // Within the moved molecule.
  for (int i = 0; i < mols[moved_mol].Size()-1; i++) {
    for (int j = i+1; j < mols[moved_mol].Size(); j++) {
      if (mols[moved_mol].bds[i].GetMoved() ||
          mols[moved_mol].bds[j].GetMoved()) {
        int id1 = mols[moved_mol].bds[i].ID();
        int id2 = mols[moved_mol].bds[j].ID();
        if (accept) {
          SetE(0, id1, id2, GetE(1, id1, id2));
        }
        else {
          SetE(1, id1, id2, GetE(0, id1, id2));
        }
      }
    }
  }

  // Between the moved molecule and the other molecules.
  for (int i = 0; i < (int)mols.size(); i++) {
    if (i != moved_mol) {
      for (int j = 0; j < mols[i].Size(); j++) {
        for (int k = 0; k < mols[moved_mol].Size(); k++) {
          if (mols[moved_mol].bds[k].GetMoved()) {
            int id1 = min(mols[moved_mol].bds[k].ID(), mols[i].bds[j].ID());
            int id2 = max(mols[moved_mol].bds[k].ID(), mols[i].bds[j].ID());
            if (accept) {
              SetE(0, id1, id2, GetE(1, id1, id2)); 
            }
            else {
              SetE(1, id1, id2, GetE(0, id1, id2));
            }
          }
        }
      }
    }
  }

}

void PotentialPair::AdjustEnergyUponMolDeletion(vector<Molecule>& mols,
                                                int delete_id) {
  // Assume monovalent counterion!
  int counterion = 0;
  for (int i = 0; i < mols[delete_id].Size(); i++) {
    counterion += (int)abs(round(mols[delete_id].bds[i].Charge()));
  }

  int gap = 0;
  if (name == "HardSphere")  gap = 1;

  ////////////
  // Chain. //
  ////////////
  for (int i = 0; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      for (int k = 0; k < mols[delete_id].Size(); k++) {
        if (i != delete_id || j > k+gap) {
          int id1 = min(mols[delete_id].bds[k].ID(), mols[i].bds[j].ID());
          int id2 = max(mols[delete_id].bds[k].ID(), mols[i].bds[j].ID());
          E_tot -= current_energy_map[make_pair(id1, id2)];
          current_energy_map.erase(make_pair(id1, id2));
          trial_energy_map.erase(make_pair(id1, id2));
        }
      }
    }
  }

  //////////////////
  // Counterions. //
  //////////////////
  for (int i = delete_id+1; i <= delete_id+counterion; i++) {
    for (int j = 0; j < (int)mols.size(); j++) {
      if (j < delete_id || j > i) {
        for (int k = 0; k < mols[j].Size(); k++) {
          int id1 = min(mols[i].bds[0].ID(), mols[j].bds[k].ID());
          int id2 = max(mols[i].bds[0].ID(), mols[j].bds[k].ID());
          E_tot -= current_energy_map[make_pair(id1, id2)];
          current_energy_map.erase(make_pair(id1, id2));
          trial_energy_map.erase(make_pair(id1, id2));
        }
      }
    }
  }

}

double PotentialPair::GetTotalEnergy() {
  return E_tot; 

}

string PotentialPair::PotentialName() {
  return name;

}

