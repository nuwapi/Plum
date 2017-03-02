#include "potential_bond.h"

using namespace std;

// Construct with initial no. of mols.
PotentialBond::PotentialBond(int n_mol, string potential_name) {
  name = potential_name;
  current_energy_array.resize(n_mol);
  trial_energy_array.resize(n_mol);
  dE = 0;
  E_tot = 0;

}

// calc all initial molecular energies, fill in both vectors. 
void PotentialBond::EnergyInitialization(vector <Molecule>& mols, double length[], int npbc) {
  E_tot = 0;

  for(int i = 0; i < (int)mols.size(); i++){
    double en = MoleculeEnergy(mols[i], length, npbc); 
    current_energy_array[i] = en; 
    trial_energy_array[i] = en; 
    E_tot += en; 
  }

}

// EnergyInitializations the last mol in the array (presumed to be new one) 
void PotentialBond::EnergyInitForLastMol(vector<Molecule>& mols, double length[], int npbc) {
  double en = MoleculeEnergy(mols[(int)mols.size()-1], length, npbc);
  current_energy_array.push_back(en); 
  trial_energy_array.push_back(en); 
  E_tot += en;  
}

void PotentialBond::AdjustEnergyUponMolDeletion(int delete_id){
  E_tot -= current_energy_array[delete_id]; 
  current_energy_array.erase(current_energy_array.begin() + delete_id); 
  trial_energy_array.erase(trial_energy_array.begin() + delete_id); 
}

void PotentialBond::SetE(int flag, int index, double val) {
  if (flag == 0) {
    current_energy_array[index] = val;
  }
  else if (flag == 1) {
    trial_energy_array[index] = val;
  }
  else {
    cout << "  Invalid array!" << endl;
  }

}

double PotentialBond::GetE(int flag, int index) {
  if (flag == 0) {
    return current_energy_array[index];
  }
  else if (flag == 1) {
    return trial_energy_array[index];
  }
  else {
    cout << "  Invalid array!" << endl;
    return 1;
  }

}
 
void PotentialBond::FinalizeEnergy(int active_mol, bool accepted) {
  if (accepted) {
    current_energy_array[active_mol] = trial_energy_array[active_mol]; 
    E_tot += dE; 
  }
  else {
    trial_energy_array[active_mol] = current_energy_array[active_mol]; 
  }

  // Energy difference has now been processed.
  dE = 0;

}

// Returns energy for all bonds in box, based on trial position.
double PotentialBond::CalcTrialTotalEnergy(vector<Molecule>& mols, double length[], int npbc) {
  double total_energy = 0;
  for (int i = 0; i < (int)mols.size(); i++) {
    total_energy += MoleculeEnergy(mols[i], length, npbc);
  }

  return total_energy;

}

double PotentialBond::GetTotalEnergy() {
  return E_tot; 

}

string PotentialBond::PotentialName() {
  return name;

}


