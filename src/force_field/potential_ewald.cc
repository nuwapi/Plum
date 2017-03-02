#include "potential_ewald.h"

#include <cmath>
#include <string>

using namespace std; 

PotentialEwald::PotentialEwald(string potential_name, double box_l_in[3]) {
  name = potential_name;
  E_tot = 0;
  current_real_E = 0;
  current_repl_E = 0;
  current_self_E = 0;
  dE = 0;
  for (int i = 0; i < 3; i++) {
    box_l[i] = box_l_in[i];
  }

  calc_pphi = false;
  trial_chain_e = new double [2];

}

PotentialEwald::~PotentialEwald() {
  delete [] trial_chain_e;

}

void PotentialEwald::SetERealRepl(int flag, int key1, int key2, double val_real,
                                  double val_repl) {
  if (flag == 0) {
    current_real_energy_map[make_pair(key1, key2)] = val_real;
    current_repl_energy_map[make_pair(key1, key2)] = val_repl;
  }
  else if (flag == 1) {
    trial_real_energy_map[make_pair(key1, key2)] = val_real;
    trial_repl_energy_map[make_pair(key1, key2)] = val_repl;
  }
  else {
    cout << "PotentialEwald::SetERealRepl:\n  Invalid map requested!" 
         << endl;
  }

}

void PotentialEwald::SetPPhiRealRepl(int flag, int key1, int key2,
                                     double val_real, double val_repl) {
  if (flag == 0) {
    current_real_pphi_map[make_pair(key1, key2)] = val_real;
    current_repl_pphi_map[make_pair(key1, key2)] = val_repl;
  }
  else if (flag == 1) {
    trial_real_pphi_map[make_pair(key1, key2)] = val_real;
    trial_repl_pphi_map[make_pair(key1, key2)] = val_repl;
  }
  else {
    cout << "PotentialEwald::SetERealRepl:\n  Invalid map requested!"
         << endl;
  }

}

void PotentialEwald::SetESelf(int flag, int key, double val_self) {
  if (flag == 0) {
    current_self_energy_map[key] = val_self;
  }
  else if (flag == 1) {
    trial_self_energy_map[key] = val_self;
  }
  else {
    cout << "PotentialEwald::SetESelf:\n  Invalid map requested!" << endl;
  }

}

void PotentialEwald::SetEBoth2DMaps(int key1, int key2, double val_real,
                                    double val_repl) {
  current_real_energy_map[make_pair(key1, key2)] = val_real; 
  current_repl_energy_map[make_pair(key1, key2)] = val_repl;
  trial_real_energy_map[make_pair(key1, key2)] = val_real;
  trial_repl_energy_map[make_pair(key1, key2)] = val_repl;

}

void PotentialEwald::SetPPhiBoth2DMaps(int key1, int key2, double val_real, 
                                       double val_repl) {
  current_real_pphi_map[make_pair(key1, key2)] = val_real;
  current_repl_pphi_map[make_pair(key1, key2)] = val_repl;
  trial_real_pphi_map[make_pair(key1, key2)] = val_real;
  trial_repl_pphi_map[make_pair(key1, key2)] = val_repl;

}

void PotentialEwald::SetEBothSelfMaps(int key, double val_self) {
  current_self_energy_map[key] = val_self;
  trial_self_energy_map[key] = val_self;

}

double PotentialEwald::GetERealRepl(int flag, int key1, int key2) {
  double ene_real, ene_repl;
  if (flag == 0) {
    ene_real = current_real_energy_map[make_pair(key1, key2)];
    ene_repl = current_repl_energy_map[make_pair(key1, key2)];
    return ene_real + ene_repl; 
  }
  else if (flag == 1) {
    ene_real = trial_real_energy_map[make_pair(key1, key2)];
    ene_repl = trial_repl_energy_map[make_pair(key1, key2)];
    return ene_real + ene_repl;
  }
  else {
    cout << "PotentialEwald::GetERealRepl:\n  Invalid flag." << endl;
    return 1; 
  }
  return 0;
 
}

double PotentialEwald::GetESelf(int flag, int key) {
 if (flag == 0) {
   return current_self_energy_map[key];
  }
  else if (flag == 1) {
   return trial_self_energy_map[key];
  }
  else {
    cout << "PotentialEwald::GetESelf:\n  Invalid flag." << endl;
    return 1; 
  }
  return 0;

}

double PotentialEwald::GetTotalEnergy() {
  return E_tot; 

}

void PotentialEwald::EnergyInitialization(vector<Molecule>& mols, int npbc) {
  double ene_real, ene_repl, ene_self;
  double pphi_real, pphi_repl;
  E_tot = 0;

  // Real and repl energies.
  for (int i = 0; i < (int)mols.size(); i++) {      // For every mol.
    for (int j = i; j < (int)mols.size(); j++) {    // For every mol.
      for (int k = 0; k < mols[i].Size(); k++) {    // For every idx in mol 1.
        for (int l = 0; l < mols[j].Size(); l++) {  // For every idx in mol 2.
          if ((j == i && l >= k) || (j > i)) {
            ene_real = PairEnergyReal(mols[i].bds[k], mols[j].bds[l], npbc);
            ene_repl = PairEnergyRepl(mols[i].bds[k], mols[j].bds[l], npbc); 
            if (j == i && l == k) {
              ene_real *= 0.5; 
              ene_repl *= 0.5;
            }
            E_tot += (ene_real + ene_repl); 
            SetEBoth2DMaps(mols[i].bds[k].ID(), mols[j].bds[l].ID(), ene_real,
                           ene_repl);

            if (calc_pphi) {
              pphi_real = PairDForceReal(mols[i].bds[k], mols[j].bds[l],
                                         mols[i].bds[0], mols[j].bds[0], npbc);
              pphi_repl = PairDForceRepl(mols[i].bds[k], mols[j].bds[l],
                                         mols[i].bds[0], mols[j].bds[0],  npbc);
              if (j == i && l == k) {
                pphi_real *= 0.5;
                pphi_repl *= 0.5;
              }
              SetPPhiBoth2DMaps(mols[i].bds[k].ID(), mols[j].bds[l].ID(),
                                ene_real, ene_repl);
            }
          }
        }
      }
    }
  }
  // Self energy.
  for (int i = 0; i < (int)mols.size(); i++) {  // For every molecule.
    for (int j = 0; j < mols[i].Size(); j++) {  // For every bead in mol.
      ene_self = SelfEnergy(mols[i].bds[j]);
      E_tot += ene_self;
      SetEBothSelfMaps(mols[i].bds[j].ID(), ene_self);
    }
  }

  // Only used when a confining potential is used.
  if (dipole_correction) {
    current_dipl_E = DipoleE(mols);
    trial_dipl_E = current_dipl_E;
    E_tot += current_dipl_E;
  }
   
}

double PotentialEwald::TrialChainEnergy(vector<Molecule>& mols,
                                        vector<Bead>& chain, int chain_len,
                                        int delete_id, int npbc) {
  double dE = 0;
  double ene_real, ene_repl, ene_self;

  // For insertion.
  int beads = 0;
  for (int i = 0; i < (int)mols.size(); i++)
    beads += mols[i].Size();
  int tot = chain_len * (2*(beads+chain_len) + 1);
  int sub = 2*(beads+chain_len) + 1;

  ///////////////////
  // If insertion. //
  ///////////////////
  if (delete_id == -1) {
    delete [] trial_chain_e;
    trial_chain_e = new double[tot];
    // Intramolecular.
    for (int i = 0; i < chain_len; i++) {
      for (int j = 0; j <= i; j++) {
        ene_real = PairEnergyReal(chain[i], chain[j], npbc);
        ene_repl = PairEnergyRepl(chain[i], chain[j], npbc);
        if (j == i) {
          ene_real *= 0.5;
          ene_repl *= 0.5;
        }
        trial_chain_e[i*sub + 2*(beads + j) + 0] = ene_real;
        trial_chain_e[i*sub + 2*(beads + j) + 1] = ene_repl;
        dE += ene_real + ene_repl;
      }
    }
    // Intermolecular.
    for (int i = 0; i < chain_len; i++) {
      int bead_counter = 0;
      for (int j = 0; j < (int)mols.size(); j++) {
        for (int k = 0; k < mols[j].Size(); k++) {
          ene_real = PairEnergyReal(chain[i], mols[j].bds[k], npbc);
          ene_repl = PairEnergyRepl(chain[i], mols[j].bds[k], npbc);
          trial_chain_e[i*sub + 2*bead_counter + 0] = ene_real;
          trial_chain_e[i*sub + 2*bead_counter + 1] = ene_repl;
          dE += ene_real + ene_repl;
          bead_counter++;
        }
      }
    }
    // Self.
    for (int i = 0; i < chain_len; i++) {
      ene_self = SelfEnergy(chain[i]);
      trial_chain_e[i*sub + 2*(beads + chain_len)] = ene_self;
      dE += ene_self;
    }
    // Only used when a confining potential is used.
    if (dipole_correction) {
      double dipl_E = DipoleE(mols, chain);
      dE += dipl_E - current_dipl_E;
    }
  }
  //////////////////
  // If deletion. //
  //////////////////
  else {
    int start = delete_id;
    int end = delete_id;
    // Meaning chain is charged.
    if (mols[delete_id].Size() < chain_len)
      end += chain_len/2;
    // Real/repl.
    for (int i = start; i <= end; i++) {
      for (int j = 0; j < mols[i].Size(); j++) {
        for (int k = 0; k < (int)mols.size(); k++) {
          for (int l = 0; l < mols[k].Size(); l++) {
            if ((k < start || k > end) || (k == i && l >= j) || (k > i)) {
              int id1 = min(mols[i].bds[j].ID(), mols[k].bds[l].ID());
              int id2 = max(mols[i].bds[j].ID(), mols[k].bds[l].ID());
              pair<int,int> indices = make_pair(id1, id2);
              ene_real = current_real_energy_map[indices];
              ene_repl = current_repl_energy_map[indices];
              dE += ene_real + ene_repl;
            }
          }
        }
      }
    }
    // Self.
    for (int i = start; i <= end; i++) {
      for (int j = 0; j < mols[i].Size(); j++) {
        ene_self = current_self_energy_map[mols[i].bds[j].ID()];
        dE += ene_self;
      }
    }
    // Dipole.
    if (dipole_correction) {
      double dipl_E = DipoleE(mols, delete_id, chain_len/2);
      dE += current_dipl_E - dipl_E;
    }
  }

  return dE;

}

void PotentialEwald::EnergyInitForLastMol(vector<Molecule>& mols, int chain_len,
                                          double bead_charge, int npbc) {
  double ene_real, ene_repl, ene_self;
  double pphi_real, pphi_repl;
  int added = 1;
  // Assume all chains are the same and always go before their counterions!!!
  // Assume monovalent ions.
  for (int i = 0; i < mols[0].Size(); i++)
    added += (int)abs(round(mols[0].bds[i].Charge()));
  int added_b = 0;
  for (int i = (int)mols.size()-added; i < (int)mols.size(); i++)
    added_b += mols[i].Size();
  int exist_b = 0;
  for (int i = 0; i < (int)mols.size()-added; i++)
    exist_b += mols[i].Size();
  int sub = 2*(exist_b + added_b) + 1;

  // Real/repl.
  int c_add = 0;
  for (int i = (int)mols.size()-added; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
     // With existing beads.
     int c_all = 0;
      for (int k = 0; k < (int)mols.size()-added; k++) {
        for (int l = 0; l < mols[k].Size(); l++) {
          int id1 = mols[k].bds[l].ID();
          int id2 = mols[i].bds[j].ID();
          ene_real = trial_chain_e[c_add*sub + 2*c_all + 0];
          ene_repl = trial_chain_e[c_add*sub + 2*c_all + 1];
          E_tot += (ene_real + ene_repl);
          SetEBoth2DMaps(id1, id2, ene_real, ene_repl);

          if (calc_pphi) {
            pphi_real = PairDForceReal(mols[k].bds[l], mols[i].bds[j],
                                       mols[k].bds[0], mols[i].bds[0], npbc);
            pphi_repl = PairDForceRepl(mols[k].bds[l], mols[i].bds[j],
                                       mols[k].bds[0], mols[i].bds[0], npbc);
            SetPPhiBoth2DMaps(id1, id2, pphi_real, pphi_repl);
          }
          c_all++;
        }
      }
      // Within added beads.
      c_all = 0;
      for (int k = (int)mols.size()-added; k <= i; k++) {
        for (int l = 0; l < mols[k].Size(); l++) {
          if (k < i || (k == i && l <= j)) {
            int id1 = mols[k].bds[l].ID();
            int id2 = mols[i].bds[j].ID();
            ene_real = trial_chain_e[c_add*sub + 2*(exist_b+c_all) + 0];
            ene_repl = trial_chain_e[c_add*sub + 2*(exist_b+c_all) + 1];
            E_tot += (ene_real + ene_repl);
            SetEBoth2DMaps(id1, id2, ene_real, ene_repl);

            if (calc_pphi) {
              pphi_real = PairDForceReal(mols[k].bds[l], mols[i].bds[j],
                                         mols[k].bds[0], mols[i].bds[0], npbc);
              pphi_repl = PairDForceRepl(mols[k].bds[l], mols[i].bds[j],
                                         mols[k].bds[0], mols[i].bds[0], npbc);
              if (i == k && j == l) {
                pphi_real *= 0.5;
                pphi_repl *= 0.5;
              }
              SetPPhiBoth2DMaps(id1, id2, pphi_real, pphi_repl);
            }
            c_all++;
          }
        }
      }

      c_add++;
    }
  }

  // Self.
  for (int i = (int)mols.size()-added; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      ene_self = SelfEnergy(mols[i].bds[j]);
      E_tot += ene_self;
      SetEBothSelfMaps(mols[i].bds[j].ID(), ene_self);
    }
  }

  // Dipole.
  if (dipole_correction) {
    current_dipl_E = DipoleE(mols);
    E_tot += current_dipl_E - trial_dipl_E;
    trial_dipl_E = current_dipl_E;
  }

}

double PotentialEwald::EnergyDifference(vector<Molecule>& mols, int active_mol,
                                        int npbc) {
  dE = 0;
  double new_ene_real, new_ene_repl, new_ene_self;
  double new_pphi_real, new_pphi_repl;

  // For all pairs in the active molecule that moved.
  for (int i = 0; i < mols[active_mol].Size(); i++) {
    for (int j = i; j < mols[active_mol].Size(); j++) {
      // Only when either i or j is moved.
      if (mols[active_mol].bds[i].GetMoved() ||
          mols[active_mol].bds[j].GetMoved()) {
        new_ene_real = PairEnergyReal(mols[active_mol].bds[i],
                                      mols[active_mol].bds[j], npbc);
        new_ene_repl = PairEnergyRepl(mols[active_mol].bds[i],
                                      mols[active_mol].bds[j], npbc);
        if (j == i) {
          new_ene_real *= 0.5;
          new_ene_repl *= 0.5;
        }

        SetERealRepl(1, mols[active_mol].bds[i].ID(),
                        mols[active_mol].bds[j].ID(),
                        new_ene_real, new_ene_repl);
        dE += (new_ene_real + new_ene_repl -
                       GetERealRepl(0, mols[active_mol].bds[i].ID(),
                                       mols[active_mol].bds[j].ID()));

        if (calc_pphi) {
          new_pphi_real = PairDForceReal(mols[active_mol].bds[i],
                                         mols[active_mol].bds[j], 
                                         mols[active_mol].bds[0],
                                         mols[active_mol].bds[0], npbc);
          new_pphi_repl = PairDForceRepl(mols[active_mol].bds[i], 
                                         mols[active_mol].bds[j],
                                         mols[active_mol].bds[0],
                                         mols[active_mol].bds[0], npbc);
          if (j == i) {
            new_pphi_real *= 0.5;
            new_pphi_repl *= 0.5;
          }
          SetPPhiRealRepl(1, mols[active_mol].bds[i].ID(),
                             mols[active_mol].bds[j].ID(),
                             new_pphi_real, new_pphi_repl);
        }

      }
    }
  }
  // For all pairs with other molecules.
  for (int i = 0; i < (int)mols.size(); i++) {
    if (i != active_mol) {
      for (int j = 0; j < mols[i].Size(); j++) {
        for (int k = 0; k < mols[active_mol].Size(); k++) {
          // Only do it for moved beads.
          if (mols[active_mol].bds[k].GetMoved()) {
            new_ene_real = PairEnergyReal(mols[active_mol].bds[k],
                                          mols[i         ].bds[j],
                                          npbc);
            new_ene_repl = PairEnergyRepl(mols[active_mol].bds[k],
                                          mols[i         ].bds[j],
                                          npbc);
            int id1 = min(mols[active_mol].bds[k].ID(),
                          mols[i         ].bds[j].ID());
            int id2 = max(mols[active_mol].bds[k].ID(),
                          mols[i         ].bds[j].ID());
            SetERealRepl(1, id1, id2, new_ene_real, new_ene_repl);
            dE += (new_ene_real + new_ene_repl -
                           GetERealRepl(0, id1, id2));

            if (calc_pphi) {
              new_pphi_real = PairDForceReal(mols[active_mol].bds[k],
                                             mols[i         ].bds[j],
                                             mols[active_mol].bds[0],
                                             mols[i         ].bds[0], npbc);
              new_pphi_repl = PairDForceRepl(mols[active_mol].bds[k],
                                             mols[i         ].bds[j],
                                             mols[active_mol].bds[0],
                                             mols[i         ].bds[0], npbc);
              SetPPhiRealRepl(1, id1, id2, new_pphi_real, new_pphi_repl);
            }
          }
        }
      }
    }
  }
  // Self-energy.
  for (int i = 0; i < mols[active_mol].Size(); i++) {
    // Only when i is moved.
    if (mols[active_mol].bds[i].GetMoved()) {
      new_ene_self = SelfEnergy(mols[active_mol].bds[i]);
      SetESelf(1, mols[active_mol].bds[i].ID(), new_ene_self);
      dE += (new_ene_self -
                     GetESelf(0, mols[active_mol].bds[i].ID()));
    }
  }


  // Only used when a confining potential is used.
  if (dipole_correction) {
    trial_dipl_E = DipoleE(mols);
    dE += trial_dipl_E - current_dipl_E;
  }
 
  return dE;

}

void PotentialEwald::FinalizeEnergyBothMaps(vector<Molecule>& mols,
                                            int active_mol, bool accepted) {
  // Adjust total energy variable.
  if (accepted) {
     E_tot += dE;
  }

  // Within the active molecule.
  for (int i = 0; i < mols[active_mol].Size(); i++) {
    for (int j = i; j < mols[active_mol].Size(); j++) {
      if (mols[active_mol].bds[i].GetMoved() ||
          mols[active_mol].bds[j].GetMoved()) {
        pair<int,int> indices = make_pair(mols[active_mol].bds[i].ID(),
                                          mols[active_mol].bds[j].ID());
        if (accepted) {
          // Update the current energy map if accepted.
          current_real_energy_map[indices] = trial_real_energy_map[indices];
          current_repl_energy_map[indices] = trial_repl_energy_map[indices];
          current_real_pphi_map[indices] = trial_real_pphi_map[indices];
          current_repl_pphi_map[indices] = trial_repl_pphi_map[indices];
        }
        else {
          // Cover up the trial energy map if rejected.
          trial_real_energy_map[indices] = current_real_energy_map[indices];
          trial_repl_energy_map[indices] = current_repl_energy_map[indices];
          trial_real_pphi_map[indices] = current_real_pphi_map[indices];
          trial_repl_pphi_map[indices] = current_repl_pphi_map[indices];
        }
      }
    }
  }
  // Between the active molecule and the other molecules.
  for (int i = 0; i < (int)mols.size(); i++) {
    if(i != active_mol) {
      for (int j = 0; j < mols[i].Size(); j++) {
        for (int k = 0; k < mols[active_mol].Size(); k++) {
          if (mols[active_mol].bds[k].GetMoved()) {
            int id1 = min(mols[active_mol].bds[k].ID(),
                          mols[i         ].bds[j].ID());
            int id2 = max(mols[active_mol].bds[k].ID(),
                          mols[i         ].bds[j].ID());
            pair<int,int> indices = make_pair(id1, id2);
            if (accepted) {
              // Update the current energy map if accepted.
              current_real_energy_map[indices] = trial_real_energy_map[indices];
              current_repl_energy_map[indices] = trial_repl_energy_map[indices];
              current_real_pphi_map[indices] = trial_real_pphi_map[indices];
              current_repl_pphi_map[indices] = trial_repl_pphi_map[indices];
            }
            else {
              // Cover up the trial energy map if rejected.
              trial_real_energy_map[indices] = current_real_energy_map[indices];
              trial_repl_energy_map[indices] = current_repl_energy_map[indices];
              trial_real_pphi_map[indices] = current_real_pphi_map[indices];
              trial_repl_pphi_map[indices] = current_repl_pphi_map[indices];
            }
          }
        }
      }
    }
  }
  // Self-energy.
  for (int i = 0; i < mols[active_mol].Size(); i++) {
    if (mols[active_mol].bds[i].GetMoved()) {
      int index = mols[active_mol].bds[i].ID();
      if (accepted) {
        // Update the current energy map if accepted.
        current_self_energy_map[index] = trial_self_energy_map[index];
      }
      else {
        // Cover up the trial energy map if rejected.
        trial_self_energy_map[index] = current_self_energy_map[index];
      }
    } 
  }


  // Only used when a confining potential is used.
  if (dipole_correction) {
    if (accepted) {
      current_dipl_E = trial_dipl_E;
    }
    else {
      trial_dipl_E = current_dipl_E;
    }
  }

}

void PotentialEwald::AdjustEnergyUponMolDeletion(vector<Molecule>& mols,
                                                 int delete_id) {
  // Assume monovalent counterion!
  int counterion = 0;
  for (int i = 0; i < mols[delete_id].Size(); i++) {
    counterion += (int)abs(round(mols[delete_id].bds[i].Charge()));
  }

  //////////////////////////
  // Chain with the rest. //
  //////////////////////////
  for (int i = 0; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      for (int k = 0; k < mols[delete_id].Size(); k++) {
        if (i != delete_id || (i == delete_id && j >= k)) {
          int id1 = min(mols[delete_id].bds[k].ID(), mols[i].bds[j].ID());
          int id2 = max(mols[delete_id].bds[k].ID(), mols[i].bds[j].ID());
          pair<int,int> indices = make_pair(id1, id2);
          E_tot -= current_real_energy_map[indices];
          E_tot -= current_repl_energy_map[indices];
          current_real_energy_map.erase(indices);
          trial_real_energy_map.erase(indices);
          current_repl_energy_map.erase(indices);
          trial_repl_energy_map.erase(indices);
          current_real_pphi_map.erase(indices);
          trial_real_pphi_map.erase(indices);
          current_repl_pphi_map.erase(indices);
          trial_repl_pphi_map.erase(indices);
        }
      }
    }
  }

  //////////////////
  // Counterions. //
  //////////////////
  for (int i = delete_id+1; i <= delete_id+counterion; i++) {
    for (int j = 0; j < (int)mols.size(); j++) {
      if (j < delete_id || j >= i) {
        for (int k = 0; k < mols[j].Size(); k++) {
          int id1 = min(mols[i].bds[0].ID(), mols[j].bds[k].ID());
          int id2 = max(mols[i].bds[0].ID(), mols[j].bds[k].ID());
          pair<int,int> indices = make_pair(id1, id2);
          E_tot -= current_real_energy_map[indices];
          E_tot -= current_repl_energy_map[indices];
          current_real_energy_map.erase(indices);
          trial_real_energy_map.erase(indices);
          current_repl_energy_map.erase(indices);
          trial_repl_energy_map.erase(indices);
          current_real_pphi_map.erase(indices);
          trial_real_pphi_map.erase(indices);
          current_repl_pphi_map.erase(indices);
          trial_repl_pphi_map.erase(indices);
        }
      }
    }
  }


  ///////////
  // Self. //
  ///////////
  for (int i = delete_id; i <= delete_id+counterion; i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      int index = mols[i].bds[j].ID();
      E_tot -= current_self_energy_map[index];
      current_self_energy_map.erase(index);
      trial_self_energy_map.erase(index);
    }
  }

  ////////////////////////
  // Dipole correction. //
  ////////////////////////
  if (dipole_correction) {
    double dipl_E = DipoleE(mols, delete_id, counterion);
    E_tot += dipl_E - current_dipl_E;
    current_dipl_E = dipl_E;
    trial_dipl_E = dipl_E;
  }

}

void PotentialEwald::UpdateEnergyComponents(vector<Molecule>& mols, int npbc) {
  current_real_E = 0;
  current_repl_E = 0;
  current_self_E = 0;

  // Real/repl energies.
  for (int i = 0; i < (int)mols.size(); i++) {
    for (int j = i; j < (int)mols.size(); j++) {
      for (int k = 0; k < mols[i].Size(); k++) {
        for (int l = 0; l < mols[j].Size(); l++) {
          if ((j == i && l >= k) || (j > i)) {
            pair<int,int> indices = make_pair(mols[i].bds[k].ID(),
                                              mols[j].bds[l].ID());
            current_real_E += current_real_energy_map[indices];
            current_repl_E += current_repl_energy_map[indices];
          }
        }
      }
    }
  }
  // Self-energy.
  for (int i = 0; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      int index = mols[i].bds[j].ID();
      current_self_E += current_self_energy_map[index];
    }
  }

}

double PotentialEwald::PUPV(vector<Molecule>& mols, double vol, int npbc) {
  double pUpV = 0;

  for (int i = 0; i < (int)mols.size(); i++) {      // For every mol.
    int i_s = mols[i].Size();
    for (int j = i; j < (int)mols.size(); j++) {    // For every mol.
      int j_s = mols[j].Size();
      for (int k = 0; k < i_s; k++) {    // For every idx in mol 1.
        for (int l = 0; l < j_s; l++) {  // For every idx in mol 2.
          if ((j == i && l >= k) || (j > i)) {
            int id1 = mols[i].bds[k].ID();
            int id2 = mols[j].bds[l].ID(); 
            pair<int,int> indices = make_pair(id1, id2);
            double phi_r = current_real_energy_map[indices];
            double phi_k = current_repl_energy_map[indices];
            double ddphi_r = current_real_pphi_map[indices];
            double ddphi_k = current_repl_pphi_map[indices];
            if (j == i && l == k) {
              ddphi_r = 0;
              ddphi_k = 0;
            }
            pUpV += - phi_r - phi_k - ddphi_r - ddphi_k;
          }
        }
      }
    }
  }

  // Add self energy.
  for (int i = 0; i < (int)mols.size(); i++) {
    int i_s = mols[i].Size();
    for (int j = 0; j < i_s; j++) {
      int id = mols[i].bds[j].ID();
      double phi_s = current_self_energy_map[id];
      pUpV += - phi_s;
    }
  }

  return pUpV / (3.0 * vol);

}

double PotentialEwald::RDotF(vector<Molecule>& mols, double vol, int npbc) {
  double r_dot_f_zz = 0;
  
  for (int i = 0; i < (int)mols.size(); i++) {
    int i_s = mols[i].Size();
    for (int j = i; j < (int)mols.size(); j++) {
      int j_s = mols[j].Size();
      for (int k = 0; k < i_s; k++) {
        for (int l = 0; l < j_s; l++) {
          if ((j == i && l >= k) || (j > i)) {
            int id1 = mols[i].bds[k].ID();
            int id2 = mols[j].bds[l].ID();
            pair<int,int> indices = make_pair(id1, id2);
            double pphi_r = current_real_pphi_map[indices];
            double pphi_k = current_repl_pphi_map[indices];
            if (j == i && l == k) {
              pphi_r = 0;
              pphi_k = 0;
            } 
            r_dot_f_zz += - pphi_r - pphi_k;
          } 
        } 
      } 
    } 
  } 

  return r_dot_f_zz;

}

double PotentialEwald::GetRealEnergy() {
  return current_real_E;

}

double PotentialEwald::GetReplEnergy() {
  return current_repl_E;

}

double PotentialEwald::GetSelfEnergy() {
  return current_self_E;

}

string PotentialEwald::PotentialName() {
  return name;

}


