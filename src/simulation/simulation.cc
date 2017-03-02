#include "simulation.h" 

#include <stdlib.h>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set> 
#include <sstream>
#include <string>
#include <vector> 

#include "../force_field/force_field.h"
#include "../molecules/bead.h"
#include "../molecules/molecule.h"
#include "../utilities/constants.h"

using namespace std; 

Simulation::Simulation() {
  cout << "  General simulation parameters." << endl;
  string flag;

  /////////////////////////////////
  // Read in and initialization. //
  /////////////////////////////////
  cin >> flag >> crd_name;
  cin >> flag >> top_name;
  cin >> flag >> run_name;
  cin >> flag >> steps;
  cin >> flag >> steps_eq;
  cin >> flag >> sample_freq;
  cin >> flag >> stat_out_freq;
  cin >> flag >> traj_out_freq;
  cin >> flag >> npbc;
  cin >> flag >> beta;
  cin >> flag >> move_size;
  cin >> flag >> calc_chem_pot;
  cin >> flag >> phantom;
  cin >> flag >> coion;
  cin >> flag >> grafted;
  cin >> flag >> grafted_counterion;
  cin >> flag >> move_prob[0];
  cin >> flag >> move_prob[1];
  cin >> flag >> move_prob[2];
  cin >> flag >> move_prob[3];
  cin >> flag >> move_prob[4];

  cout << setw(35) << "Input coordinate file       : " << crd_name      << endl;
  cout << setw(35) << "Input topology file         : " << top_name      << endl;
  cout << setw(35) << "Prefix of output file names : " << run_name      << endl;
  cout << setw(35) << "Total simulation steps      : " << steps         << endl;
  cout << setw(35) << "Equilibrium steps           : " << steps_eq      << endl;
  cout << setw(35) << "Sampling frequency          : " << sample_freq   << endl;
  cout << setw(35) << "Statistics output frequency : " << stat_out_freq << endl;
  cout << setw(35) << "Trajectory output frequency : " << traj_out_freq << endl;
  cout << setw(35) << "Number of dimension of PBC  : " << npbc          << endl;
  cout << setw(35) << "MC move size (ul)           : " << move_size     << endl;
  cout << setw(35) << "beta (1/kBT)                : " << beta          << endl;
  cout << setw(35) << "Calculate chemical pot?     : " << YesOrNo(calc_chem_pot)
                                                       << endl;
  cout << setw(35) << "Number of phantom molecules : " << phantom       << endl;
  cout << setw(35) << "Number of coions            : " << coion         << endl;
  cout << setw(35) << "p(bead translation move)    : " << move_prob[0]  << endl;
  cout << setw(35) << "p(pol COM translation move) : " << move_prob[1]  << endl;
  cout << setw(35) << "p(pol pivot move)           : " << move_prob[2]  << endl;
  cout << setw(35) << "p(pol crankshaft move)      : " << move_prob[3]  << endl;
  cout << setw(35) << "p(pol random reptation move): " << move_prob[4]  << endl;

  // Set up variables for simulation statistics.
  step = 0;
  id_counter = 0;
  mol_avg_counter = 0;
  ff_avg_counter = 0;
  chem_pot_cumu = 0;
  pair_e_cumu = 0;
  ewald_e_cumu = 0;
  ewald_e_real_cumu = 0;
  ewald_e_repl_cumu = 0;
  ewald_e_self_cumu = 0;
  bond_e_cumu = 0;
  ext_e_cumu = 0;
  rg_tot_cumu = 0;
  rg_x_cumu = 0;
  rg_y_cumu = 0;
  rg_z_cumu = 0;
  e_to_e_cumu = 0;
  adsorbed_chains = 0;
  adsorbed_beads = 0;
  adsorption_percent = 0;
  for (int i = 0; i < kNoMoveType; i++) {
    accepted[i] = 0;
    attempted[i] = 0;
  }
  insertion_accepted = 0;
  deletion_accepted = 0;
  insertion_attempted = 0;
  deletion_attempted = 0;

  // Read bead connectivity info and initialize bead/molecule numbers.
  ReadTop();
  // Read coordinates.
  ReadCrd();
  // Set up ForceField object.
  force_field.Initialize(beta, npbc, box_l, mols, phantom, coion, grafted,
                         grafted_counterion);

  // Densities.
  density_cumu = 0;
  density_z_res = 0.02;
  density_z_bin = (int)ceil(box_l[2]/density_z_res);
  density_z_cumu = new double[3*density_z_bin];
  for (int i = 0; i < 3*density_z_bin; i++)
    density_z_cumu[i] = 0;

  // Open output files.
  string info_out_name = run_name + "_stat.dat";
  string traj_p_out_name = run_name + "_traj_p.xyz";
  string traj_c_out_name = run_name + "_traj_c.xyz";
  string traj_a_out_name = run_name + "_traj_a.xyz";
  string traj_n_out_name = run_name + "_traj_n.xyz";
  info_out.open(info_out_name.c_str());
  traj_p_out.open(traj_p_out_name.c_str());
  traj_c_out.open(traj_c_out_name.c_str());
  traj_a_out.open(traj_a_out_name.c_str());
  traj_n_out.open(traj_n_out_name.c_str());
  // Print header for output statistics file.
  PrintStatHeader();

  // Initialize random number generator using time.
  unsigned s = chrono::system_clock::now().time_since_epoch().count();
  rand_gen.seed(s);
  cout << "  The random number seed is: " << s << endl;

  ///////////////////
  // A few checks. //
  ///////////////////
  if (beta <= 0) {
    cout << "  " << beta << " is not an acceptable beta value. Exiting! "
         << "Program complete." << endl;
    exit(1);
  }
  else if (beta != 1.0) {
    cout << "  Note: It is highly recommended that beta be set to 1 " << endl;
    cout << "        such that the energy unit of the simulation is " << endl;
    cout << "        kBT!" << endl;
  }
  if (npbc > 3 || npbc < 2) {
    cout << "  " << npbc << " is an invalid number of periodic "
         << "dimensions! Exiting! Program complete." << endl;
    exit(1);
  }
  // Check move probabilities.
  double tot_move_prob = 0;
  for (int i = 0; i < kNoMoveType; i++) {
    tot_move_prob += move_prob[i];
  }
  if (abs(tot_move_prob-1.0) > kMedSmallNumber) {
    cout << "  MC move probabilities do not sum up to 1.0. Exiting! Program "
         << "complete." << endl;
    exit(1);
  }
  if (grafted + n_chain == 0 && 1-move_prob[0] > 0) {
    cout << "  Warning: There is no chain in the system but the chain" << endl;
    cout << "           MC move(s) is requested!" << endl;
  }
  if (n_mol - grafted - n_chain == 0 && move_prob[0] > 0) {
    cout << "  Warning: There is no single bead in the system but the" << endl; 
    cout << "           bead MC move is requested!" << endl;
  }
  // If crankshaft move and rigid bond are used at the same time.
  if (move_prob[3] > 0 && force_field.UseBondRigid()) {
    cout << "  The current crankshaft MC move and rigid bond cannot be used at "
         << "the same time. Exiting! Program complete." << endl;
    //exit(1);
  }
  // If reptation and rigid bond are not used at the same time.
  if (move_prob[4] > 0 && (!force_field.UseBondRigid())) {
    cout << "  The current reptation MC move and rigid bond have to be used at "
         << "the same time. Exiting! Program complete." << endl;
    exit(1);
  }
  if (calc_chem_pot && force_field.UseGC()) {
    cout << "  ERROR: Chemical potential cannot be calculate when GCMC" << endl;
    cout << "         is used. Exiting. Program complete." << endl;
    exit(1);
  }
  // Note for phantom beads.
  if (phantom > 0) {
    cout << "  Note: [!!!] "<< phantom << " phantom beads are used. Check" << endl;
    cout << "        if they are used with the appropriate parameters!" << endl;
    cout << "        Their bead type should have < 0 radius and their" << endl;
    cout << "        coordinates should appear at the BEGINING of the" << endl;
    cout << "        input coordiante file." << endl;
  }

}

Simulation::~Simulation() {
  PrintLastCrd();
  PrintLastTop();
  info_out.close();
  traj_p_out.close();
  traj_c_out.close();
  traj_a_out.close();
  traj_n_out.close();
  delete [] density_z_cumu;
  cout << "\n  Simulation complete." << endl;

}

void Simulation::Run() {
  int gc_freq = force_field.GCFrequency();
  UpdateMolCounts();

  for (step = 1; step <= steps; step++) {
    int rand_num = rand_gen();
    // Randomly choose whether to do a GC move.
    if (force_field.UseGC() && (rand_num % gc_freq == 0)) {
      GCMove();
    }
    // If not doing a GC move, do a translational move.
    else {
      TranslationalMove();
    }

    Sample();
    PrintStat();
    PrintTraj();
    PrintLastCrd();
    PrintLastTop();
    PrintLastRhoZ();
  }

}

void Simulation::TranslationalMove() {
  if (n_mol - phantom > 0) {
    // Pre-choose a chain and an ion.
    int chain_id = -1;
    int ion_id = -1;
    int mol_id = -1;
    int which_chain = (int)floor((double)rand_gen()/rand_gen.max() *
                                 (grafted + n_chain));
    int which_ion   = (int)floor((double)rand_gen()/rand_gen.max() *
                                 (n_cion + n_aion + n_nion));
    if (which_chain == grafted + n_chain)       which_chain--;
    if (which_ion == n_cion + n_aion + n_nion)  which_ion--;
    int chain_choosing_counter = -1;
    int ion_choosing_counter = -1;
    for (int i = phantom; i < n_mol; i++) {
      if (mols[i].Size() > 1)  chain_choosing_counter++;
      else                     ion_choosing_counter++;
      if (chain_choosing_counter == which_chain) {
        chain_id = i;
        which_chain = -2;
      }
      if (ion_choosing_counter == which_ion) {
        ion_id = i;
        which_ion = -2;
      }
    }

    // Decide move type.
    int move_type = 0;
    double rand_num = (double)rand_gen() / rand_gen.max(); 
    double current_move_type = move_prob[0];
    while (current_move_type < rand_num) {
      move_type++;
      current_move_type += move_prob[move_type];
    }
    if (move_type >= kNoMoveType)  move_type = kNoMoveType - 1;
  
    // If choose to move a single bead and there is one.
    if (move_type == 0 && n_cion + n_aion + n_nion > 0) {
      mol_id = ion_id;
      mols[mol_id].BeadTranslate(move_size, box_l, rand_gen);
    }
    // If choose to move a chain and there is one.
    else if (grafted + n_chain > 0) {
      mol_id = chain_id;
      switch (move_type) {
        case 1:
          mols[mol_id].COMTranslate(move_size, rand_gen);
          break;
        case 2:
          mols[mol_id].Pivot(move_size, rand_gen, force_field.RigidBondLen());
          break;
        case 3:
          mols[mol_id].Crankshaft(move_size, rand_gen);
          break;
        case 4:
          mols[mol_id].RandomReptation(rand_gen, force_field.RigidBondLen());
          break;
      }
    }

    // If a move is actually attempted.
    if ((move_type == 0 && n_cion + n_aion + n_nion > 0) ||
        (move_type > 0 && grafted + n_chain > 0)) {
      double dE = 0;
      bool accept = false;
      attempted[move_type]++;

      // Assign kVeryLargeEnergy if the polymer leaves the z surfaces.
      // dE = force_field.EnsureGrafting(mols, mol_id);
      // Only calculate energy difference if grafting is ensured.
      if (dE < kVeryLargeEnergy)
        dE = force_field.EnergyDifference(mols, mol_id);

      if (dE >= kVeryLargeEnergy) {
        accept = false;
      }
      else {
        accept = ((double)rand_gen() / rand_gen.max() < exp(-beta*dE));
      }

      force_field.FinalizeEnergies(mols, accept, mol_id);

      // Make trial positions current positions.
      if (accept) {
        for (int i = 0; i < mols[mol_id].Size(); i++) {
          mols[mol_id].bds[i].UpdateCurrentPos();   
        }
        accepted[move_type]++; 
      }
      // Set trial positions to the current positions.
      else {
        for (int i = 0; i < mols[mol_id].Size(); i++) {
          mols[mol_id].bds[i].UpdateTrialPos();
        }
      }
      for (int i = 0; i < mols[mol_id].Size(); i++) {
        mols[mol_id].bds[i].UnsetMoved(); 
      }
    }
  }

}

void Simulation::GCMove() {
  // Choose whether to attempt insertion or deletion.
  bool insert = rand_gen() % 2;

  // Attempt insertion.
  if (insert) {
    if (step > steps_eq)  insertion_attempted++;
    bool accept = force_field.CBMCFChainInsertion(mols, rand_gen);
    if (accept) {
      if (step > steps_eq)  insertion_accepted++;
      for (int i = n_mol; i < (int)mols.size(); i++) {
        for (int j = 0; j < mols[i].Size(); j++) {
          mols[i].bds[j].SetID(GenBeadID());
          mols[i].bds[j].SetChainID(i);
        }
      }
      force_field.EnergyInitForAddedMolecule(mols);
    }
  }
  // Attempt deletion.
  else {
    if (step > steps_eq)  deletion_attempted++;
    int deleted_chain;

    // If there is no mobile pairs, there is nothing to delete.
    if (n_mol - phantom - coion <= 0) {
      deleted_chain = -1;
    }
    else {
      deleted_chain = force_field.CBMCFChainDeletion(mols, rand_gen);
      if (deleted_chain >= phantom + coion) {
        if (step > steps_eq)  deletion_accepted++;

        int deleted_ions = 0;
        for (int i = 0; i < mols[deleted_chain].Size(); i++)
          deleted_ions += (int)round(abs(mols[deleted_chain].bds[i].Charge()));

        for (int i = deleted_chain + deleted_ions; i >= deleted_chain; i--) {
          mols.erase(mols.begin()+i);
        }
      }
    }
  }

  UpdateMolCounts();
  force_field.UpdateMolCounts(mols);

}

double Simulation::RadiusOfGyration() {
  double rgAvg = 0;

  if ((int)mols.size() > 0) {
    for (int i = 0; i < (int)mols.size(); i++) {
      // If this is a chain molecule.
      if (mols[i].Size() > 1) {
        double rgMol = 0;
        for (int j = 0; j < mols[i].Size()-1; j++) {
          for (int k = j+1; k < mols[i].Size(); k++) {
            // Setting pbc_dimesion to 0 to prevent folding of a polymer into
            // the same unit cell.
            double r = mols[i].bds[j].BBDist(mols[i].bds[k],
                                                       box_l, 0);
            rgMol += (r*r);
          }
        }
        rgMol /= (mols[i].Size() * mols[i].Size());
        rgAvg += rgMol;
      }
    }
    if (grafted + n_chain > 0) {
      rgAvg = rgAvg / (grafted + n_chain);
    }
  }

  return rgAvg;

}

// Nuo added: 08/09/2016.
string Simulation::RadiusOfGyrationXYZ() {
  string rgAvgXYZ = "";
  double rgAvgX = 0;
  double rgAvgY = 0;
  double rgAvgZ = 0;

  if ((int)mols.size() > 0) {
    for (int i = 0; i < (int)mols.size(); i++) {
      // If this is a chain molecule.
      if (mols[i].Size() > 1) {
        double rgMolX = 0;
        double rgMolY = 0;
        double rgMolZ = 0;
        for (int j = 0; j < mols[i].Size()-1; j++) {
          for (int k = j+1; k < mols[i].Size(); k++) {
            double x = mols[i].bds[j].NaiveBBDistX(mols[i].bds[k]);
            double y = mols[i].bds[j].NaiveBBDistY(mols[i].bds[k]);
            double z = mols[i].bds[j].NaiveBBDistZ(mols[i].bds[k]);
            rgMolX += (x*x); 
            rgMolY += (y*y);
            rgMolZ += (z*z);
          }
        }
        rgMolX /= (mols[i].Size() * mols[i].Size()); 
        rgMolY /= (mols[i].Size() * mols[i].Size());
        rgMolZ /= (mols[i].Size() * mols[i].Size());
        rgAvgX += rgMolX;
        rgAvgY += rgMolY;
        rgAvgZ += rgMolZ; 
      }
    }

    if (grafted + n_chain > 0) {
      rgAvgX = rgAvgX / (grafted + n_chain);
      rgAvgY = rgAvgY / (grafted + n_chain);
      rgAvgZ = rgAvgZ / (grafted + n_chain);
    }
  }

  std::ostringstream foo;
  foo << rgAvgX << " " << rgAvgY << " " << rgAvgZ;
  return foo.str();

}

// Nuo added: 09/27/2016.
double Simulation::EndToEndDistance() {
  double avg_ete = 0;
  if ((int)mols.size() > 0) {
    for (int i = 0; i < (int)mols.size(); i++) {
      // If this is a chain molecule.
      if (mols[i].Size() > 1) {
        double dist = mols[i].bds[0].BBDist(
                      mols[i].bds[mols[i].Size()-1], box_l, 0);
        avg_ete += dist*dist;
      }
    }
    if (grafted + n_chain > 0) {
      avg_ete /= (grafted + n_chain);
    }
  }

  return avg_ete;

}

void Simulation::PrintEndToEndVector() {
  if ((int)mols.size() > 0) {
    cout << "X";
    for (int i = 0; i < (int)mols.size(); i++) {
      if (mols[i].Size() > 1) {
        double dimex = mols[i].bds[0               ].GetCrd(0, 0) -
                       mols[i].bds[mols[i].Size()-1].GetCrd(0, 0);
        double dimey = mols[i].bds[0               ].GetCrd(0, 1) -
                       mols[i].bds[mols[i].Size()-1].GetCrd(0, 1);
        double dimez = mols[i].bds[0               ].GetCrd(0, 2) -
                       mols[i].bds[mols[i].Size()-1].GetCrd(0, 2);
        double dist = sqrt(dimex*dimex+dimey*dimey+dimez*dimez);
        cout << dist << " " << dimex << " " << dimey << " " << dimez << " ";
      }
    }
    cout << endl;
  }

}

// Nuo added: 09/07/2016.
void Simulation::DetectAdsorption() {
  int current_adsorbed_chain = 0;

  if ((int)mols.size() > 0) {
    // For every molecule.
    for (int i = 0; i < (int)mols.size(); i++) {
      // If this is a chain molecule.
      if (mols[i].Size() > 1) {
        bool chain_adsorbed_left = false;
        bool chain_adsorbed_right = false;
        // For every bead in the chain.
        for (int j = 0; j < mols[i].Size(); j++) {
          istringstream iss(mols[i].bds[j].DistToWall(box_l));
          double dist_to_left, dist_to_right;
          iss >> dist_to_left;
          iss >> dist_to_right;

          if (dist_to_left < kAdsorbCutoff) {
            chain_adsorbed_left = true;
            adsorbed_beads++;
          }
          if (dist_to_right < kAdsorbCutoff) {
            chain_adsorbed_right = true;
            adsorbed_beads++;
          }

        }
        if (chain_adsorbed_left == true) {
          adsorbed_chains++;
          current_adsorbed_chain++;
        }
        if (chain_adsorbed_right == true) {
          adsorbed_chains++;
        }
      }
    }
  }

  adsorption_percent += (double)current_adsorbed_chain/(double)mols.size();

}

int Simulation::GenBeadID() {
  pair<set<int>::iterator, bool> result; 
  bool inserted = false;

  while (!inserted) {
    result = id_list.insert(id_counter);
    inserted = result.second;
    id_counter++;
  }

  // At end of loop, id_counter equals to the next id to be inserted in the
  // future, so minus 1 for the just inserted id.
  return id_counter - 1;

}

void Simulation::ReadCrd() {
  crd_in.open(crd_name.c_str());

  int mol_id;
  double x, y, z;
  double charge;
  string symbol;
  string line;

  getline(crd_in, line);
  getline(crd_in, line); 
  while (crd_in >> mol_id >> symbol >> x >> y >> z >> charge) {
    if (mol_id < 0 || mol_id > n_mol) {
      cout << "  Molecule ID " << mol_id << " is out of bounds!" << endl;
      cout << "  Ending program. Program complete. :-(" << endl;
      exit(1);
    }
    // If molecule ID is valid, create bead ID and add bead to that molecule.
    else {
      int id = GenBeadID();
      mols[mol_id].AddBead(Bead(symbol, id, mol_id, charge, x, y, z));
    }
  }

  /*
  // Count number of chain molecules.
  n_chain = 0;
  for (int i = phantom; i < (int)mols.size(); i++) {
    if (mols[i].Size() > 1) {
      n_chain++; 
    }
  }
  */

  crd_in.close();

}

void Simulation::ReadTop() {
  top_in.open(top_name.c_str());

  string flag;
  top_in >> flag >> n_bead;
  top_in >> flag >> n_mol;
  top_in >> flag >> box_l[0];
  top_in >> flag >> box_l[1];
  top_in >> flag >> box_l[2];

  int ind1, ind2, ind3, ind4, mol_id;
  top_in >> flag >> mol_id; 
  // Populate mol array with initial number of empty molecules.
  for(int i = 0; i < n_mol; i++){
    mols.push_back(Molecule()); 
  }
  while (mol_id != -1) {
    top_in >> ind1 >> ind2; 
    mols[mol_id].AddBond(ind1, ind2); 
    top_in >> mol_id; 
  }
  top_in >> flag >> mol_id; 
  while (mol_id != -1) {
    top_in >> ind1 >> ind2 >> ind3; 
    mols[mol_id].AddAngle(ind1, ind2, ind3); 
    top_in >> mol_id; 
  }
  top_in >> flag >> mol_id; 
  while (mol_id != -1) {
    top_in >> ind1 >> ind2 >> ind3 >> ind4; 
    mols[mol_id].AddDihed(ind1, ind2, ind3, ind4); 
    top_in >> mol_id; 
  }

  top_in.close();

}

void Simulation::CoordinateObeyRigidBond(double rigid_bond) {
  // Maintain rigid bond length.
  for (int i = phantom; i < n_mol; i++) {
    for (int j = 0; j< mols[i].Size()-1; j++) {
      double x1 = mols[i].bds[j].GetCrd(0, 0);
      double y1 = mols[i].bds[j].GetCrd(0, 1);
      double z1 = mols[i].bds[j].GetCrd(0, 2);
      double x2 = mols[i].bds[j+1].GetCrd(0, 0);
      double y2 = mols[i].bds[j+1].GetCrd(0, 1);
      double z2 = mols[i].bds[j+1].GetCrd(0, 2);
      double dx = x2-x1;
      double dy = y2-y1;
      double dz = z2-z1;
      double dist = sqrt(dx*dx + dy*dy + dz*dz);
      double incr = 0;
      double cdist = dist;
      while (cdist <= rigid_bond) {
        x2 = x1 + (rigid_bond+incr) * dx / dist;
        y2 = y1 + (rigid_bond+incr) * dy / dist;
        z2 = z1 + (rigid_bond+incr) * dz / dist;
        dx = x2-x1;
        dy = y2-y1;
        dz = z2-z1;
        cdist = sqrt(dx*dx + dy*dy + dz*dz);
        incr += kSmallNumber;
      }
      mols[i].bds[j+1].SetCrd(0, 0, x2);
      mols[i].bds[j+1].SetCrd(1, 0, x2);
      mols[i].bds[j+1].SetCrd(0, 1, y2);
      mols[i].bds[j+1].SetCrd(1, 1, y2);
      mols[i].bds[j+1].SetCrd(0, 2, z2);
      mols[i].bds[j+1].SetCrd(1, 2, z2);
    }
  }

}

void Simulation::PrintStatHeader() {
  info_out << "#Step";
  if (force_field.UsePairPot()) {
    info_out << " " << "<PairEne>";
  }
  if (force_field.UseEwaldPot()) {
    info_out << " " << "<EwaldEne>";
  }
  if (force_field.UseBondPot()) {
    info_out << " " << "<BondEne>";
  }
  if (force_field.UseExtPot()) {
    info_out << " " << "<ExtPotEne>";
  }
  info_out << " " << "<TotalEnergy>";
  if (force_field.UseGC()) {
    info_out << " " << "NoOfMol";
    info_out << " " << "<Density>";
  }
  info_out << " " << "<Pxx> <Pyy> <Pzz> <Phxx> <Phyy> <Phzz> <Pexx> <Peyy> <Pezz>";
  if (calc_chem_pot) {
    info_out << " " << "mu";
  }
  info_out << " " << "<|Rg|> <|Rgx|> <|Rgy|> <|Rgz|> <ete>";
  info_out << " " << "beadtrans comtrans pivot crankshaft reptate";
  if (force_field.UseGC()) {
    info_out << " " << "insert delete";
  }
  info_out << endl;

}

void Simulation::Sample() {
  if (step > steps_eq && step % sample_freq == 0) {
    ff_avg_counter++;
    //////////////
    // Energies //
    //////////////
    if (force_field.UsePairPot()) {
      pair_e_cumu += force_field.TotPairEnergy();
    }
    if (force_field.UseEwaldPot()) {
      ewald_e_cumu += force_field.TotEwaldEnergy();
      //double * temp = new double[3];
      //force_field.GetEwaldEnergyComponents(mols, temp);
      //ewald_e_real_cumu += temp[0];
      //ewald_e_repl_cumu += temp[1];
      //ewald_e_self_cumu += temp[2];
      //delete [] temp;
    }
    if (force_field.UseBondPot()) {
      bond_e_cumu += force_field.TotBondEnergy();
    }
    if (force_field.UseExtPot()) {
      ext_e_cumu += force_field.TotExtEnergy();
    }
    /////////////
    // Density //
    /////////////
    density_cumu += (double)mols.size() / (box_l[0]*box_l[1]*box_l[2]);
    CalcRhoZ();
    //////////////
    // Pressure //
    //////////////

    if (force_field.UseExtPot() && step % (sample_freq*10) == 0) {
      force_field.CalcPressureVolScalingHSELSlit(mols);
      //force_field.CalcPressureVirialHSELSlit(mols, density_cumu/ff_avg_counter);
      //force_field.CalcPressureForceELSlit(mols);
    }
    else if (step % (sample_freq*10) == 0) {
      force_field.CalcPressureVolScalingHSELSlit(mols);
      //force_field.CalcPressureVirialHSEL(mols, density_cumu/ff_avg_counter);
      //force_field.CalcPressureVirialEL(mols);
    }
    /////////////////////////
    // Chemical potential //
    ////////////////////////
    if (calc_chem_pot) {
      chem_pot_cumu += force_field.CalcChemicalPotentialF(mols, rand_gen);
    }
  }

  if (step > steps_eq && step % sample_freq == 0 && (int)mols.size() > 0) {
    mol_avg_counter++;
    ////////////////////////
    // Radius of gyration //
    ////////////////////////
    istringstream iss(RadiusOfGyrationXYZ());
    double rg_x_cumu_in, rg_y_cumu_in, rg_z_cumu_in;
    iss >> rg_x_cumu_in;
    iss >> rg_y_cumu_in;
    iss >> rg_z_cumu_in;
    rg_tot_cumu += RadiusOfGyration();
    rg_x_cumu += rg_x_cumu_in;
    rg_y_cumu += rg_y_cumu_in;
    rg_z_cumu += rg_z_cumu_in;
    /////////////////////////
    // End to end distance //
    /////////////////////////
    e_to_e_cumu += EndToEndDistance();
    ////////////////
    // Adsorption //
    ////////////////
    /*
    if (force_field.UseExtPot()) {
      DetectAdsorption();
    }
    */
  }

}

void Simulation::PrintStat() {
  if (step % stat_out_freq == 0) {
    info_out << step;

    double tot_energy = 0;
    if (force_field.UsePairPot()) {
      info_out << " " << setprecision(7) << pair_e_cumu / ff_avg_counter;
      tot_energy += pair_e_cumu / ff_avg_counter;
    }
    if (force_field.UseEwaldPot()) {
      info_out << " " << setprecision(7) << ewald_e_cumu / ff_avg_counter;
      //info_out << " " << setprecision(12) << ewald_e_real_cumu / ff_avg_counter;
      //info_out << " " << setprecision(12) << ewald_e_repl_cumu / ff_avg_counter;
      //info_out << " " << setprecision(12) << ewald_e_self_cumu / ff_avg_counter;
      tot_energy += ewald_e_cumu / ff_avg_counter;
    }
    if (force_field.UseBondPot()) {
      info_out << " " << setprecision(7) << bond_e_cumu / ff_avg_counter;
      tot_energy += bond_e_cumu / ff_avg_counter;
    }
    if (force_field.UseExtPot()) {
      info_out << " " << setprecision(7) << ext_e_cumu / ff_avg_counter;
      tot_energy += ext_e_cumu / ff_avg_counter;
    }
    info_out << " " << tot_energy;
    if (force_field.UseGC()) {
      info_out << " " << (int)mols.size()
               << " " << setprecision(4) << density_cumu/ff_avg_counter;
    }
    info_out << " " << force_field.GetPressure();
    if (calc_chem_pot) {
      if (chem_pot_cumu/ff_avg_counter > 0)
        info_out << " " << -(1.0/beta)*log(chem_pot_cumu/ff_avg_counter);
      else
        info_out << " INF";
    }

    info_out << " " << setprecision(4) << sqrt(rg_tot_cumu/mol_avg_counter)
             << " " << setprecision(4) << sqrt(rg_x_cumu/mol_avg_counter)
             << " " << setprecision(4) << sqrt(rg_y_cumu/mol_avg_counter)
             << " " << setprecision(4) << sqrt(rg_z_cumu/mol_avg_counter)
             << " " << setprecision(4) << sqrt(e_to_e_cumu/mol_avg_counter);
    /*
    if (force_field.UseExtPot()) {
      // To convert to count per square of average radius of gyration.
      // 2 represents two walls.
      double area_convert = (rg_tot_cumu/mol_avg_counter) * (rg_tot_cumu/mol_avg_counter)
                            / (2*box_l[0]*box_l[1]);
      info_out << " " << setprecision(4)
               << area_convert*(1.0*adsorbed_chains)/mol_avg_counter
               << " " << setprecision(4)
               << area_convert*(1.0*adsorbed_beads)/mol_avg_counter
               << " " << setprecision(4) << adsorption_percent/mol_avg_counter;
    }
    */
    info_out << " " << accepted[0] / (double)attempted[0];
    info_out << " " << accepted[1] / (double)attempted[1];
    info_out << " " << accepted[2] / (double)attempted[2];
    info_out << " " << accepted[3] / (double)attempted[3];
    info_out << " " << accepted[4] / (double)attempted[4];
    if (force_field.UseGC()) {
      info_out << " " << insertion_accepted / (double)insertion_attempted;
      info_out << " " << deletion_accepted / (double)deletion_attempted;
    }
    info_out << endl;
  }

}

void Simulation::PrintTraj() {
  if (step > steps_eq && step % traj_out_freq == 0) {
    if (grafted + n_chain > 0) {
      traj_p_out << n_chain_b << endl;
      traj_p_out << "STEP: " << step << endl;
    }
    if (n_cion > 0) {
      traj_c_out << n_cion << endl;
      traj_c_out << "STEP: " << step << endl;
    }
    if (n_aion > 0) {
      traj_a_out << n_aion << endl;
      traj_a_out << "STEP: " << step << endl;
    }
    if (n_nion > 0) {
      traj_n_out << n_nion << endl;
      traj_n_out << "STEP: " << step << endl;
    }

    for (int i = phantom; i < n_mol; i++) {
      double centerx = 0;
      double centery = 0;
      double centerz = 0;
      for (int j = 0; j < mols[i].Size(); j++) {
        centerx += mols[i].bds[j].GetCrd(0, 0);
        centery += mols[i].bds[j].GetCrd(0, 1);
        centerz += mols[i].bds[j].GetCrd(0, 2);
      }
      centerx /= (double)mols[i].Size();
      centery /= (double)mols[i].Size();
      centerz /= (double)mols[i].Size();
      double dispx = floor(centerx/box_l[0])*box_l[0];
      double dispy = floor(centery/box_l[1])*box_l[1];
      double dispz = floor(centerz/box_l[2])*box_l[2];

      if (mols[i].Size() > 1) {
        for (int j = 0; j < mols[i].Size(); j++) {
          traj_p_out << mols[i].bds[j].Symbol() << " ";
          traj_p_out << mols[i].bds[j].GetCrd(0, 0) - dispx << " ";
          traj_p_out << mols[i].bds[j].GetCrd(0, 1) - dispy << " ";
          traj_p_out << mols[i].bds[j].GetCrd(0, 2) - dispz << endl;
        }
      }
      else if (mols[i].bds[0].Charge() > 0) {
        traj_c_out << mols[i].bds[0].Symbol() << " ";
        traj_c_out << mols[i].bds[0].GetCrd(0, 0) - dispx << " ";
        traj_c_out << mols[i].bds[0].GetCrd(0, 1) - dispy << " ";
        traj_c_out << mols[i].bds[0].GetCrd(0, 2) - dispz << endl;
      }
      else if (mols[i].bds[0].Charge() < 0) {
        traj_a_out << mols[i].bds[0].Symbol() << " ";
        traj_a_out << mols[i].bds[0].GetCrd(0, 0) - dispx << " ";
        traj_a_out << mols[i].bds[0].GetCrd(0, 1) - dispy << " ";
        traj_a_out << mols[i].bds[0].GetCrd(0, 2) - dispz << endl;
      }
      else {
        traj_n_out << mols[i].bds[0].Symbol() << " ";
        traj_n_out << mols[i].bds[0].GetCrd(0, 0) - dispx << " ";
        traj_n_out << mols[i].bds[0].GetCrd(0, 1) - dispy << " ";
        traj_n_out << mols[i].bds[0].GetCrd(0, 2) - dispz << endl;
      }
    }
  }

}

void Simulation::PrintLastCrd() {
  string bck[2] = {"", "2"};
  if (step > steps_eq && step % stat_out_freq == 0) {
    for (int f = 0; f < 2; f++) {
      crd_last_out[f].open(run_name + "_lastcrd.dat" + bck[f]);
  
      crd_last_out[f] << n_bead << endl;
      crd_last_out[f] << " " << endl; 
      for (int i = 0; i < n_mol; i++) {
        double centerx = 0;
        double centery = 0;
        double centerz = 0;
        for (int j = 0; j < mols[i].Size(); j++) {
          centerx += mols[i].bds[j].GetCrd(0, 0);
          centery += mols[i].bds[j].GetCrd(0, 1);
          centerz += mols[i].bds[j].GetCrd(0, 2);
        }
        centerx /= (double)mols[i].Size();
        centery /= (double)mols[i].Size();
        centerz /= (double)mols[i].Size();
        double dispx = floor(centerx/box_l[0])*box_l[0];
        double dispy = floor(centery/box_l[1])*box_l[1];
        double dispz = floor(centerz/box_l[2])*box_l[2];
        if (i < phantom) {
          dispx = dispy = dispz = 0;
        }
  
        for (int j = 0; j < mols[i].Size(); j++) {
          crd_last_out[f] << i << " ";
          crd_last_out[f] << mols[i].bds[j].Symbol() << " ";
          crd_last_out[f] << setprecision(18)
                       << mols[i].bds[j].GetCrd(0, 0) - dispx << " ";
          crd_last_out[f] << setprecision(18)
                       << mols[i].bds[j].GetCrd(0, 1) - dispy << " ";
          crd_last_out[f] << setprecision(18)
                       << mols[i].bds[j].GetCrd(0, 2) - dispz << " ";
          crd_last_out[f] << setprecision(4) << mols[i].bds[j].Charge() << endl;
        }
      }
  
      crd_last_out[f].close();
    }
  }

}

void Simulation::PrintLastTop() {
  if (step > steps_eq && step % stat_out_freq == 0) {
    top_last_out.open(run_name + "_lasttop.dat");

    top_last_out << "TotNoOfBeads: " << n_bead      << endl; 
    top_last_out << "TotNoOfMolec: " << n_mol       << endl; 
    top_last_out << "Box_Length_X: " << setprecision(18) << box_l[0] <<endl; 
    top_last_out << "Box_Length_Y: " << setprecision(18) << box_l[1] <<endl; 
    top_last_out << "Box_Length_Z: " << setprecision(18) << box_l[2] <<endl; 
    top_last_out << "bonds:" << endl; 
    /*
    for (int i = 0; i < n_mol; i++) {
      for (int j = 0; j < (int)mols[i].bonds.size(); j++) {
        top_last_out << i << " " << mols[i].bonds[j][0] << " "
                                 << mols[i].bonds[j][1] << endl; 
      }
    }
    */
    top_last_out << "-1" << endl; 
    top_last_out << "angles:" << endl; 
    /*
    for (int i = 0; i < n_mol; i++) {
      for(int j = 0; j < (int)mols[i].angles.size(); j++) {
        top_last_out << i << " " << mols[i].angles[j][0] << " "
                                 << mols[i].angles[j][1] << " "
                                 << mols[i].angles[j][2] << endl; 
      }
    } 
    */ 
    top_last_out << "-1" << endl; 
    top_last_out << "diheds" << endl; 
    /*
    for (int i = 0; i < n_mol; i++) {
      for(int j = 0; j < (int)mols[i].diheds.size(); j++) {
        top_last_out << i << " " << mols[i].diheds[j][0] << " "
                                 << mols[i].diheds[j][1] << " "
                                 << mols[i].diheds[j][2] << " "
                                 << mols[i].diheds[j][3] << endl;
      }
    }
    */
    top_last_out << "-1" << endl; 

    top_last_out.close();
  }

}

void Simulation::CalcRhoZ() {
  for (int i = phantom; i < n_mol; i++) {
    int id;
    if (mols[i].Size() > 1)                 id = 0;
    else if (mols[i].bds[0].Charge() >= 0)  id = 1;
    else                                    id = 2;
 
    for (int j = 0; j < mols[i].Size(); j++) {
      double z = mols[i].bds[j].GetCrd(0,2);
      z -= box_l[2] * floor(z / box_l[2]);
      int binz = (int)floor(z / density_z_res);
      if (binz >= 0 && binz < density_z_bin)
        density_z_cumu[3*binz + id] += 1.0/(box_l[0]*box_l[1]*density_z_res);
    }
  }

}

void Simulation::PrintLastRhoZ() {
  string bck[2] = {"", "2"};
  if (step > steps_eq && step % stat_out_freq == 0) {
    for (int i = 0; i < 2; i++) {
      rho_last_out[i].open(run_name + "_lastrhoZ.dat" + bck[i]);

      for (int bin = 0; bin < density_z_bin; bin++) {
        double pos = (bin + 0.5) * density_z_res;
        rho_last_out[i] << pos << '\t'
                        << density_z_cumu[3*bin + 0]/(double)ff_avg_counter << '\t'
                        << density_z_cumu[3*bin + 1]/(double)ff_avg_counter << '\t'
                        << density_z_cumu[3*bin + 2]/(double)ff_avg_counter << endl;
      }
      rho_last_out[i].close();
    }
  }

}

void Simulation::UpdateMolCounts() {
  n_mol = (int)mols.size();
  n_chain = 0;
  n_chain_b = 0;
  n_cion = 0;
  n_aion = 0;
  n_nion = 0;
  n_bead = phantom;  // Each phantom bead is one molecule.

  for (int i = phantom; i < n_mol; i++) {
    n_bead += mols[i].Size();
    if (mols[i].Size() > 1) {
      n_chain++;
      n_chain_b += mols[i].Size();
    }
    else if (mols[i].bds[0].Charge() > 0) {
      n_cion++;
    }
    else if (mols[i].bds[0].Charge() < 0) {
      n_aion++;
    }
    else {
      n_nion++;
    }
  }

  n_chain -= grafted;

}


