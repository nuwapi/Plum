#include "force_field.h"

#include "../utilities/constants.h"

double ForceField::BeadsEnergy(Bead& bead1, Bead& bead2, vector<Molecule>& mols,
                               int current_len, int delete_id) {
  double energy   = 0;
  double pair_e   = 0;  double ewald_e  = 0;
  double pair_e1  = 0;  double pair_e2  = 0;
  double ewald_r1 = 0;  double ewald_r2 = 0;
  double ewald_k1 = 0;  double ewald_k2 = 0;
  int counterion = 0;
  if (gc_bead_charge != 0)  counterion = gc_chain_len;

  /////////////////////////////////////////////////////////////////
  // 1. Pairwise interactions for the bead with other molecules. //
  /////////////////////////////////////////////////////////////////
  for (int i = 0; i < (int)mols.size(); i++) {
    if (delete_id == -1 || (i < delete_id || i > delete_id + counterion)) {
      for (int j = 0; j < mols[i].Size(); j++) {
        if (use_pair_pot) {
          pair_e1 = pair_pot->PairEnergy(bead1, mols[i].bds[j], box_l, npbc);
          if (gc_bead_charge != 0)
            pair_e2 = pair_pot->PairEnergy(bead2, mols[i].bds[j], box_l, npbc);
        }
        pair_e += pair_e1 + pair_e2;

        if (pair_e >= kVeryLargeEnergy)  break;

        if (use_ewald_pot) {
          ewald_r1 = ewald_pot->PairEnergyReal(bead1, mols[i].bds[j], npbc);
          ewald_k1 = ewald_pot->PairEnergyRepl(bead1, mols[i].bds[j], npbc);
          if (gc_bead_charge != 0) {
            ewald_r2 = ewald_pot->PairEnergyReal(bead2, mols[i].bds[j], npbc);
            ewald_k2 = ewald_pot->PairEnergyRepl(bead2, mols[i].bds[j], npbc);
          }
        }
        ewald_e += (ewald_r1 + ewald_k1 + ewald_r2 + ewald_k2);
      }
    }
    if (pair_e >= kVeryLargeEnergy)  break;
  }

  ////////////////////////////////////////////////////////
  // 2. Loop over the beads in the current trial chain. //
  ////////////////////////////////////////////////////////
  for (int i = 0; i < current_len; i++) {
    if (use_pair_pot) {
      pair_e1 = 0;
      if (i < current_len - 1)
        pair_e1 += pair_pot->PairEnergy(bead1, cbmc_chain[i], box_l, npbc);
      if (gc_bead_charge != 0) {
        pair_e1 += pair_pot->PairEnergy(bead1, cbmc_chain[i+gc_chain_len],
                                        box_l, npbc);
        pair_e2  = pair_pot->PairEnergy(bead2, cbmc_chain[i], box_l, npbc);
        pair_e2 += pair_pot->PairEnergy(bead2, cbmc_chain[i+gc_chain_len], 
                                        box_l, npbc);
      }
    }
    pair_e += pair_e1 + pair_e2;

    if (pair_e >= kVeryLargeEnergy)  break;

    if (use_ewald_pot) {
      ewald_r1  = ewald_pot->PairEnergyReal(bead1, cbmc_chain[i], npbc);
      ewald_k1  = ewald_pot->PairEnergyRepl(bead1, cbmc_chain[i], npbc);
      if (gc_bead_charge != 0) {
        ewald_r1 += ewald_pot->PairEnergyReal(bead1, cbmc_chain[i+gc_chain_len],
                                              npbc);
        ewald_k1 += ewald_pot->PairEnergyRepl(bead1, cbmc_chain[i+gc_chain_len], 
                                              npbc);
        ewald_r2  = ewald_pot->PairEnergyReal(bead2, cbmc_chain[i], npbc);
        ewald_k2  = ewald_pot->PairEnergyRepl(bead2, cbmc_chain[i], npbc);
        ewald_r2 += ewald_pot->PairEnergyReal(bead2, cbmc_chain[i+gc_chain_len],
                                              npbc);
        ewald_k2 += ewald_pot->PairEnergyRepl(bead2, cbmc_chain[i+gc_chain_len],
                                              npbc);
     }
    }
    ewald_e += (ewald_r1 + ewald_k1 + ewald_r2 + ewald_k2);

  }

  //////////////////////////////////////////////////////
  // 3. Inter bead pair & Ewald self/dipole energies. //
  //////////////////////////////////////////////////////
  if (use_pair_pot && gc_bead_charge != 0) {
    pair_e += pair_pot->PairEnergy(bead1, bead2, box_l, npbc);
  }
  if (use_ewald_pot && pair_e < kVeryLargeEnergy) {
    ewald_e += ewald_pot->SelfEnergy(bead1);
    ewald_e += ewald_pot->SelfEnergy(bead2);
    ewald_e += ewald_pot->PairEnergyReal(bead1, bead2, npbc);
    ewald_e += ewald_pot->PairEnergyRepl(bead1, bead2, npbc);
    if (ewald_pot->UseDipoleCorrection())
      ewald_e += ewald_pot->DipoleEDiff(mols, cbmc_chain, bead1, bead2,
                                        current_len, gc_chain_len,
                                        gc_bead_charge,  delete_id);
  }

  /////////////////////////////////////////////////
  // 4. Interaction with the external hard wall. //
  /////////////////////////////////////////////////
  if (use_ext_pot) {
    pair_e += ext_pot->BeadEnergy(bead1, box_l);
    if (gc_bead_charge != 0)
      pair_e += ext_pot->BeadEnergy(bead2, box_l);
  }

  ////////////////
  // 5. Return. //
  ////////////////
  energy = pair_e + ewald_e;
  if (pair_e >= kVeryLargeEnergy)
    return kVeryLargeEnergy;
  return energy;
 
}

/*
 * Input
 * end_bead   : the previously chosen bead in the trial chain.
 * mols       : molecule array of the existing beads.
 * current_len: current length of the trial chain.
 * rand_gen   : random number generator.
 * delete_id  : used for chain deletion, the ID of the chain to be deleted.
 */
double ForceField::CBMCFGenTrialBeads(Bead& end_bead, vector<Molecule>& mols,
                                      int current_len, mt19937& rand_gen,
                                      int delete_id) {
  double Wi = 0;

  for (int i = 0; i < cbmc_no_of_trials; i++) {
    /////////////////////////////////////
    // Decide indices and bond length. //
    /////////////////////////////////////
    int c_index = i;  // Chain trial bead index.
    int i_index = i;  // Ion trial bead index.
    if (gc_bead_charge != 0)  i_index = i + cbmc_no_of_trials;
    double bond_len = 0;
    if (use_bond_pot) {
      bond_len = bond_pot->RandomBondLen(beta, rand_gen); 
    }
    else if (use_bond_rigid) {
      bond_len = rigid_bond;
    }

    //////////////////////////////
    // Decide bead coordinates. //
    //////////////////////////////
    double bead_coord[3];
    // Chain particle.
    randSphere(bead_coord, rand_gen);
    for (int j = 0; j < 3; j++) {
      bead_coord[j] *= bond_len;
      bead_coord[j] += end_bead.GetCrd(0, j);
    }
    cbmc_trial_beads[c_index].SetAllCrd(bead_coord);
    // Chain counterion.
    if (gc_bead_charge != 0) {
      bead_coord[0] = (double)rand_gen()/rand_gen.max() * box_l[0];
      bead_coord[1] = (double)rand_gen()/rand_gen.max() * box_l[1];
      bead_coord[2] = (double)rand_gen()/rand_gen.max() * box_l[2];
      cbmc_trial_beads[i_index].SetAllCrd(bead_coord);
    }

    //////////////////////////////
    // Calculating bead energy. //
    //////////////////////////////
    double bead_energy = BeadsEnergy(cbmc_trial_beads[c_index],
                                     cbmc_trial_beads[i_index],
                                     mols, current_len, delete_id);
    cbmc_trial_weights[i] = exp(-beta * bead_energy);    
    Wi += cbmc_trial_weights[i];
  }

  return Wi;

}

bool ForceField::CBMCFChainInsertion(vector<Molecule>& mols,
                                     mt19937& rand_gen) {
  bool accept = false;
  double weight = 1.0;

  ///////////////////
  // Grow the 1st. //
  ///////////////////
  double xyz[3];
  int i_index = 0;  // Ion index.
  int initial_beads = 1;
  if (gc_bead_charge != 0) {
    i_index = gc_chain_len;
    initial_beads = 2;
  }
  for (int i = 0; i < initial_beads; i++) {
    for (int j = 0; j < 3; j++)
      xyz[j] = (double)rand_gen() / rand_gen.max() * box_l[j];
    cbmc_chain[i*gc_chain_len].SetAllCrd(xyz);
  }
  weight *= exp(-beta * BeadsEnergy(cbmc_chain[0], cbmc_chain[i_index], mols, 0,
                                    -1));
  if (weight <= 0)  return false;


  ////////////////////
  // Grow the rest. //
  ////////////////////
  for (int i = 1; i < gc_chain_len; i++) {
    // Generate beads calculate their energies.
    double Wi = CBMCFGenTrialBeads(cbmc_chain[i-1], mols, i, rand_gen, -1);
    weight *= Wi/cbmc_no_of_trials;
    if (weight <= 0)  return accept;

    // Choose bead(s).
    double rand_num = (double)rand_gen() / rand_gen.max() * Wi;
    int current_bead = 0;
    double cumulate_weight = cbmc_trial_weights[0];
    while (cumulate_weight < rand_num) {
      current_bead++;
      cumulate_weight += cbmc_trial_weights[current_bead];
    }

    // Assign coordinates.
    for (int j = 0; j < 3; j++)
      xyz[j] = cbmc_trial_beads[current_bead].GetCrd(0, j);
    cbmc_chain[i].SetAllCrd(xyz);
    if (gc_bead_charge != 0) {
      current_bead += cbmc_no_of_trials;
      for (int j = 0; j < 3; j++) 
        xyz[j] = cbmc_trial_beads[current_bead].GetCrd(0, j);
      cbmc_chain[i + gc_chain_len].SetAllCrd(xyz);
    }
  }

  //////////////////////////////////
  // Decide acceptance/rejection. //
  //////////////////////////////////
  // Calculate the prefactor used in the acceptance rule.
  UpdateMolCounts(mols);
  int spc1;
  if      (gc_chain_len   >  1)  spc1 = n_chain;
  else if (gc_bead_charge >= 0)  spc1 = n_cion - coion;
  else                           spc1 = n_aion - coion;
  int spc2;
  if      (gc_bead_charge == 0)  spc2 = 0;
  else if (gc_bead_charge >  0)  spc2 = n_aion;
  else                           spc2 = n_cion;
  int spc1_add = 1;
  int spc2_add = gc_chain_len;
  if      (gc_bead_charge == 0)  spc2_add = 0;
  double factorial = 1;
  double m1 = gc_chain_len;
  double m2 = 1;
  double vol_over_lam = 1;
  for (int i = spc1 + 1; i <= spc1 + spc1_add; i++)
    factorial *= i;
  for (int i = spc2 + 1; i <= spc2 + spc2_add; i++)
    factorial *= i;
  vol_over_lam *= pow(vol/(gc_deBroglie_prefactor/pow(m1, 1.5)), spc1_add);
  vol_over_lam *= pow(vol/(gc_deBroglie_prefactor/pow(m2, 1.5)), spc2_add);
  double C = vol_over_lam * (1.0/(double)factorial);

  // Decide acceptance or rejection.
  double rand_num = (double)rand_gen() / rand_gen.max();

  if (rand_num < (exp(beta*chem_pot) * weight) * C) {
    accept = true;
    // Initialize Ewald energy.
    if (use_ewald_pot)
      ewald_pot->TrialChainEnergy(mols, cbmc_chain, gc_chain_len*2, -1, npbc);

    mols.push_back(Molecule());
    for (int i = 0; i < gc_chain_len; i++) {
      mols[(int)mols.size() - 1].AddBead(cbmc_chain[i]); 
    }
    // Add bonds to the new chain ASSUMING LINEAR MOLECULE!
    for (int i = 0; i < mols[(int)mols.size() - 1].Size()-1; i++) {
      mols[(int)mols.size() - 1].AddBond(i, i+1); 
    }
    // Add counterions.
    if (gc_bead_charge != 0) {
      for (int i = gc_chain_len; i < gc_chain_len*2; i++) {
        mols.push_back(Molecule());
        mols[(int)mols.size() - 1].AddBead(cbmc_chain[i]);
      }
    }
  }

  // At this point, the new beads are added with the correct coordinates, but
  // they still need the proper IDs and energy initialization, also Simulation
  // nParticle fields need to be updated.
  return accept; 

}

int ForceField::CBMCFChainDeletion(vector<Molecule>& mols, mt19937& rand_gen) {
  UpdateMolCounts(mols);
  int spc1;
  if      (gc_chain_len   >  1)  spc1 = n_chain;
  else if (gc_bead_charge >= 0)  spc1 = n_cion - coion;
  else                           spc1 = n_aion - coion;

  int delete_id;
  delete_id = floor((double)rand_gen()/rand_gen.max() * spc1);
  if (delete_id >= spc1)  delete_id = spc1 - 1;
  if (gc_bead_charge != 0)
    delete_id = phantom + coion + delete_id * (1 + gc_chain_len);
  else
    delete_id = phantom + coion + delete_id;

  int second_bead = 0;
  if (gc_bead_charge != 0)  second_bead = 1;
  
  ////////////////////////
  // "Grow" first bead. //
  ////////////////////////
  double weight = 1.0;
  // Weight for the first bead in chain.
  int second_mol = delete_id;
  if (gc_bead_charge != 0)  second_mol = delete_id + 1;
  weight *= exp(-beta*BeadsEnergy(mols[delete_id].bds[0],
                                  mols[second_mol].bds[0],
                                  mols, 0, delete_id));
  // Replace the xyz bead by the existing bead on the chain to be deleted.
  double xyz[3] = {mols[delete_id].bds[0].GetCrd(0, 0),
                   mols[delete_id].bds[0].GetCrd(0, 1),
                   mols[delete_id].bds[0].GetCrd(0, 2)};
  cbmc_chain[0].SetAllCrd(xyz);
  if (gc_bead_charge != 0) {
    xyz[0] = mols[second_mol].bds[0].GetCrd(0, 0);
    xyz[1] = mols[second_mol].bds[0].GetCrd(0, 1);
    xyz[2] = mols[second_mol].bds[0].GetCrd(0, 2);
    cbmc_chain[gc_chain_len].SetAllCrd(xyz);
  }

  ////////////////////
  // Grow the rest. //
  ////////////////////
  for (int i = 1; i < gc_chain_len; i++) {
    xyz[0] = mols[delete_id].bds[i].GetCrd(0, 0);
    xyz[1] = mols[delete_id].bds[i].GetCrd(0, 1);
    xyz[2] = mols[delete_id].bds[i].GetCrd(0, 2);
    cbmc_chain[i].SetAllCrd(xyz);
    if (gc_bead_charge != 0) {
      xyz[0] = mols[delete_id+i+1].bds[0].GetCrd(0, 0);
      xyz[1] = mols[delete_id+i+1].bds[0].GetCrd(0, 1);
      xyz[2] = mols[delete_id+i+1].bds[0].GetCrd(0, 2);
      cbmc_chain[i+gc_chain_len].SetAllCrd(xyz);
    }

    CBMCFGenTrialBeads(cbmc_chain[i-1], mols, i, rand_gen, delete_id);
    cbmc_trial_weights[0] = exp(-beta * BeadsEnergy(cbmc_chain[i],
                            cbmc_chain[i+second_bead*gc_chain_len], mols, i,
                            delete_id));
    double Wi = 0;
    for (int j = 0; j < cbmc_no_of_trials; j++) {
      Wi += cbmc_trial_weights[j];
    }
    weight *= Wi/cbmc_no_of_trials;
  }

  //////////////////////////////////
  // Decide acceptance/rejection. //
  //////////////////////////////////
  int spc2;
  if      (gc_bead_charge == 0)  spc2 = 0;
  else if (gc_bead_charge >  0)  spc2 = n_aion;
  else                           spc2 = n_cion;
  int spc1_del = 1;
  int spc2_del = gc_chain_len;
  if      (gc_bead_charge == 0)  spc2_del = 0;
  double factorial = 1;
  double m1 = gc_chain_len;
  double m2 = 1;
  double lam_over_vol = 1;
  for (int i = spc1; i > spc1 - spc1_del; i--)
    factorial *= i;
  for (int i = spc2; i > spc2 - spc2_del; i--)
    factorial *= i;
  lam_over_vol *= pow((gc_deBroglie_prefactor/pow(m1, 1.5))/vol, spc1_del);
  lam_over_vol *= pow((gc_deBroglie_prefactor/pow(m2, 1.5))/vol, spc2_del);
  double C = lam_over_vol * (double)factorial;
  double rand_num = (double)rand_gen() / rand_gen.max();

  if (rand_num < C / (exp(beta*chem_pot) * weight)) {
    if (use_pair_pot) {
      pair_pot->AdjustEnergyUponMolDeletion(mols, delete_id); 
    }
    if (use_ewald_pot) {
      ewald_pot->AdjustEnergyUponMolDeletion(mols, delete_id);
    }
    if (use_bond_pot) {
      bond_pot->AdjustEnergyUponMolDeletion(delete_id); 
    }
    if (use_ext_pot) {
      ext_pot->AdjustEnergyUponMolDeletion(mols, delete_id); 
    }
    return delete_id; 
  }

  // Signal unsuccessful deletion.
  return -1;

}

double ForceField::CalcChemicalPotentialF(vector<Molecule>& mols,
                                          mt19937& rand_gen) {
  double total = 0;

  int spc1;
  if      (gc_chain_len   >  1)  spc1 = n_chain;
  else if (gc_bead_charge >= 0)  spc1 = n_cion - coion;
  else                           spc1 = n_aion - coion;
  int spc2;
  if      (gc_bead_charge == 0)  spc2 = 0;
  else if (gc_bead_charge >  0)  spc2 = n_aion;
  else                           spc2 = n_cion;
  int spc1_add = 1;
  int spc2_add = gc_chain_len;
  if      (gc_bead_charge == 0)  spc2_add =0;
  double factorial = 1;
  double m1 = gc_chain_len;
  double m2 = 1;
  double vol_over_lam = 1;
  for (int i = spc1 + 1; i <= spc1 + spc1_add; i++)
    factorial *= i;
  for (int i = spc2 + 1; i <= spc2 + spc2_add; i++)
    factorial *= i;
  vol_over_lam *= pow(vol/(gc_deBroglie_prefactor/pow(m1, 1.5)), spc1_add);
  vol_over_lam *= pow(vol/(gc_deBroglie_prefactor/pow(m2, 1.5)), spc2_add);
  double C = vol_over_lam * (1.0/(double)factorial);

  int toadd = gc_chain_len;
  if (gc_bead_charge != 0)  toadd *= 2;

  for (int c = 0; c < mu_tot_ins; c++) {
    double weight = 1.0;

    // First beads.
    double xyz[3];
    int i_index = 0;
    int initial_beads = 1;
    if (gc_bead_charge != 0) {
      i_index = gc_chain_len;
      initial_beads = 2;
    }
    for (int i = 0; i < initial_beads; i++) {
      for (int j = 0; j < 3; j++)
        xyz[j] = (double)rand_gen() / rand_gen.max() * box_l[j];
      cbmc_chain[i*gc_chain_len].SetAllCrd(xyz);
    }
    weight *= exp(-beta * BeadsEnergy(cbmc_chain[0], cbmc_chain[i_index], mols, 0,
                                      -1));

    // The rest.
    for (int i = 1; i < gc_chain_len; i++) {
      // Generate beads calculate their energies.
      double Wi = CBMCFGenTrialBeads(cbmc_chain[i-1], mols, i, rand_gen, -1);
      weight *= Wi/cbmc_no_of_trials;
      if (weight <= 0)  break;
  
      // Choose bead(s).
      double rand_num = (double)rand_gen() / rand_gen.max() * Wi;
      int current_bead = 0;
      double cumulate_weight = cbmc_trial_weights[0];
      while (cumulate_weight < rand_num) {
        current_bead++;
        cumulate_weight += cbmc_trial_weights[current_bead];
      }
  
      // Assign coordinates.
      for (int j = 0; j < 3; j++)
        xyz[j] = cbmc_trial_beads[current_bead].GetCrd(0, j);
      cbmc_chain[i].SetAllCrd(xyz);
      if (gc_bead_charge != 0) {
        current_bead += cbmc_no_of_trials;
        for (int j = 0; j < 3; j++)
          xyz[j] = cbmc_trial_beads[current_bead].GetCrd(0, j);
        cbmc_chain[i + gc_chain_len].SetAllCrd(xyz);
      }
    }

    total += weight * C;
  }
  
  return total/(double)mu_tot_ins;

}


