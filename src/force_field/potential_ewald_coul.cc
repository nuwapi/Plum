#include "potential_ewald_coul.h"

#include <cmath>
#include <iomanip>
#include <vector>

#include "../molecules/molecule.h"
#include "../utilities/constants.h"
#include "../utilities/misc.h"

PotentialEwaldCoul::PotentialEwaldCoul(string potential_name, double box_l[3])
                                     : PotentialEwald(potential_name, box_l) {
  ReadParameters(); 

}

PotentialEwaldCoul::~PotentialEwaldCoul() {
  delete [] kx;
  delete [] ky;
  delete [] kz;
  delete [] k2;
  delete [] ek2;
  delete [] kz_forP;
  delete [] k2_forP;
  delete [] ek2_forP;

}

void PotentialEwaldCoul::ReadParameters() {
  cout << setw(35) << "[EP] Ewald potential type   : " << "Coul" << endl;

  string flag;
  cin >> flag >> lB
      >> flag >> dielectric
      >> flag >> alpha
      >> flag >> dipole_correction;
  // e2 * kKe / (dielectric * lB), in unit J*A/(unit length)
  kBT_derived = kKe / (dielectric * lB);
  // In unit K*A/(unit length)
  T_derived = kBT_derived / kKB;
  // Determining box size for Ewald, which may not be the same as the rest
  // of the simulation.
  if (dipole_correction) {
    box_l[2] = box_l[2]*5;  // Why 5 times? See Frenkel and Smit.
  }
  box_vol = box_l[0]*box_l[1]*box_l[2];
  if (dipole_correction) {
    box_vol_forP = box_l[0]*box_l[1]*(box_l[2]+5*kDz);
  }
  else {
    box_vol_forP = box_l[0]*box_l[1]*(box_l[2]+kDz);
  }

  // Automatically decide real space cutoff. Using unit charge as a criteria.
  real_cutoff = 1;
  while (0.5*lB*1*1*erfc(sqrt(alpha)*real_cutoff)/real_cutoff
          > kEwaldCutoff) {
    // Add one unit length at a time.
    real_cutoff += 1;
  }
  for (int i = 0; i < 3; i++) {
    real_cell[i] = (int)ceil(real_cutoff / box_l[i]);
  }
  // Automaticall decide reciprocal space cutoff.
  repl_cutoff = alpha;
  while (lB*1*1/(2*kPi*box_vol) * (4*kPi*kPi) * exp(-repl_cutoff/(4*alpha))
         / repl_cutoff > kEwaldCutoff) {
    repl_cutoff += alpha;
  }
  for (int i = 0; i < 3; i++) {
    repl_cell[i] = (int)ceil(sqrt(repl_cutoff)*box_l[i]/(2*kPi));
  }

  // Setting up the variables to improve the repl calculation efficiency.
  repl_ceto[0] = 2*repl_cell[0] + 1;
  repl_ceto[1] = 2*repl_cell[1] + 1;
  repl_ceto[2] = 2*repl_cell[2] + 1;
  kx = new double[repl_ceto[0]];
  ky = new double[repl_ceto[1]];
  kz = new double[repl_ceto[2]];
  k2 = new double[repl_ceto[0]*repl_ceto[1]*repl_ceto[2]];
  ek2 = new double[repl_ceto[0]*repl_ceto[1]*repl_ceto[2]];
  kz_forP = new double[repl_ceto[2]];
  k2_forP = new double[repl_ceto[0]*repl_ceto[1]*repl_ceto[2]];
  ek2_forP = new double[repl_ceto[0]*repl_ceto[1]*repl_ceto[2]];
  for (int lx = -repl_cell[0]; lx <= repl_cell[0]; lx++) {
    kx[lx+repl_cell[0]] = lx * 2*kPi / box_l[0];
  }
  for (int ly = -repl_cell[1]; ly <= repl_cell[1]; ly++) {
    ky[ly+repl_cell[1]] = ly * 2*kPi / box_l[1];
  }
  for (int lz = -repl_cell[2]; lz <= repl_cell[2]; lz++) {
    kz[lz+repl_cell[2]] = lz * 2*kPi / box_l[2];
    kz_forP[lz+repl_cell[2]] = lz * 2*kPi / (box_l[2] + kDz);
  }
  for (int lx = -repl_cell[0]; lx <= repl_cell[0]; lx++) {
    for (int ly = -repl_cell[1]; ly <= repl_cell[1]; ly++) {
      for (int lz = -repl_cell[2]; lz <= repl_cell[2]; lz++) {
        int idx = lx+repl_cell[0];
        int idy = ly+repl_cell[1];
        int idz = lz+repl_cell[2];
        int index = repl_ceto[1]*repl_ceto[2]*idx + repl_ceto[2]*idy + idz;
        k2[index] = kx[idx]*kx[idx] + ky[idy]*ky[idy] + kz[idz]*kz[idz];
        k2_forP[index] = kx[idx]*kx[idx] + ky[idy]*ky[idy] + kz_forP[idz]*kz_forP[idz];
        ek2[index] = exp(-k2[index]/(4*alpha)) / k2[index];
        ek2_forP[index] = exp(-k2_forP[index]/(4*alpha)) / k2_forP[index];
      }
    }
  }
  ////////////////////////

  cout << setw(35) << "[EP] Bjerrum length (ul)    : " << lB          << endl;
  cout << setw(35) << "[EP] Dielectric constant    : " << dielectric  << endl;
  cout << setw(35) << "[EP] Derived kBT in J*A/ul  : " << kBT_derived << endl;
  cout << setw(35) << "[EP] Derived T K*A/ul       : " << T_derived   << endl;
  cout << setw(35) << "[EP] Ewald alpha (ul^-2)    : " << alpha       << endl;
  cout << setw(35) << "[EP] Real space cutoff (ul) : " << real_cutoff << endl;
  cout << setw(35) << "[EP] (cell sum number)      : "
       << ceil((real_cell[0]+real_cell[1]+real_cell[2])/3.0) << endl;
  cout << setw(35) << "[EP] k space cutoff (ul^-2) : " << repl_cutoff << endl;
  cout << setw(35) << "[EP] (k sum number)         : "
       << ceil((repl_cell[0]+repl_cell[1]+repl_cell[2])/3.0) << endl;
  cout << setw(35) << "[EP] Use dipole correction? : " 
       << YesOrNo(dipole_correction) << endl;

}

double PotentialEwaldCoul::PairEnergyReal(Bead& bead1, Bead& bead2, int npbc) {
  double energy = 0;
  double dist[3];
  GetDistVector(bead1, bead2, box_l, npbc, dist);
  double q1 = bead1.Charge();
  double q2 = bead2.Charge();
  double prefactor = lB*q1*q2;

  for (int i = -real_cell[0]; i <= real_cell[0]; i++) {
    for (int j = -real_cell[1]; j <= real_cell[1]; j++) {
      for (int k = -real_cell[2]; k <= real_cell[2]; k++) {
        double r_vec[3];
        r_vec[0] = dist[0] + i*box_l[0];
        r_vec[1] = dist[1] + j*box_l[1];
        r_vec[2] = dist[2] + k*box_l[2];
        double r = sqrt(r_vec[0]*r_vec[0]+r_vec[1]*r_vec[1]+r_vec[2]*r_vec[2]);
        if (r > 0 && r <= real_cutoff) {
          energy += prefactor * erfc(sqrt(alpha)*r)/r;
        }
      }
    }
  }

  return energy;

}

// PairEnergyReal for scaled volume.
double PotentialEwaldCoul::PairEnergyRealForP(Bead& bead1, Bead& bead2, int npbc) {
  double energy = 0;
  double dist[3];
  double box_l_scaled[3] = {box_l[0], box_l[1], box_l[2]+kDz};
  GetDistVector(bead1, bead2, box_l_scaled, npbc, dist);
  double q1 = bead1.Charge();
  double q2 = bead2.Charge();
  double prefactor = lB*q1*q2;

  for (int i = -real_cell[0]; i <= real_cell[0]; i++) {
    for (int j = -real_cell[1]; j <= real_cell[1]; j++) {
      for (int k = -real_cell[2]; k <= real_cell[2]; k++) {
        double r_vec[3];
        r_vec[0] = dist[0] + i*box_l[0];
        r_vec[1] = dist[1] + j*box_l[1];
        if (dipole_correction)
          r_vec[2] = dist[2] + k*(box_l[2]+kDz*5);
        else
          r_vec[2] = dist[2] + k*(box_l[2]+kDz);
        double r = sqrt(r_vec[0]*r_vec[0]+r_vec[1]*r_vec[1]+r_vec[2]*r_vec[2]);
        if (r > 0 && r <= real_cutoff) {
      
          energy += prefactor * erfc(sqrt(alpha)*r)/r;
        }
      }
    }
  }

/*
  if (q1 == q2 && q1 == 1)
    cout << "test!!!" << energy << " ";
  energy = 0;
  GetDistVectorC(bead1, bead2, box_l, npbc, dist);

  for (int i = -real_cell[0]; i <= real_cell[0]; i++) {
    for (int j = -real_cell[1]; j <= real_cell[1]; j++) {
      for (int k = -real_cell[2]; k <= real_cell[2]; k++) {
        double r_vec[3];
        r_vec[0] = dist[0] + i*box_l[0];
        r_vec[1] = dist[1] + j*box_l[1];
        if (dipole_correction)
          r_vec[2] = dist[2] + k*(box_l[2]);
        else
          r_vec[2] = dist[2] + k*(box_l[2]);
        double r = sqrt(r_vec[0]*r_vec[0]+r_vec[1]*r_vec[1]+r_vec[2]*r_vec[2]);
        if (r > 0 && r <= real_cutoff) {
          energy += prefactor * erfc(sqrt(alpha)*r)/r;
        }
      }
    }
  }
  if (q1 == q2 && q1 == 1)
    cout << energy << endl;;
*/

  return energy;

}

double PotentialEwaldCoul::PairEnergyRepl(Bead& bead1, Bead& bead2, int npbc) {
  double energy = 0;
  double r[3];
  GetDistVector(bead1, bead2, box_l, npbc, r);
  double q1 = bead1.Charge();
  double q2 = bead2.Charge();
  double prefactor = lB * q1*q2/(kPi*box_vol) * (4*kPi*kPi);

  for (int lx = 0; lx < repl_ceto[0]; lx++) {
    for (int ly = 0; ly < repl_ceto[1]; ly++) {
      for (int lz = 0; lz < repl_ceto[2]; lz++) {
        int idx = repl_ceto[1]*repl_ceto[2]*lx + repl_ceto[2]*ly + lz;
        if (k2[idx] > 0 && k2[idx] <= repl_cutoff) {
          energy += prefactor * ek2[idx]
                    * cos(kx[lx]*r[0] + ky[ly]*r[1] + kz[lz]*r[2]);
        }
      }
    }
  }

  return energy;

}

// PairEnergyRepl for volume scaling.
double PotentialEwaldCoul::PairEnergyReplForP(Bead& bead1, Bead& bead2, int npbc) {
  double energy = 0;
  double r[3];
  double box_l_scaled[3] = {box_l[0], box_l[1], box_l[2]+kDz};
  GetDistVector(bead1, bead2, box_l_scaled, npbc, r);
  double q1 = bead1.Charge();
  double q2 = bead2.Charge();
  double prefactor = lB * q1*q2/(kPi*box_vol_forP) * (4*kPi*kPi);

  for (int lx = 0; lx < repl_ceto[0]; lx++) {
    for (int ly = 0; ly < repl_ceto[1]; ly++) {
      for (int lz = 0; lz < repl_ceto[2]; lz++) {
        int idx = repl_ceto[1]*repl_ceto[2]*lx + repl_ceto[2]*ly + lz;
        if (k2_forP[idx] > 0 && k2_forP[idx] <= repl_cutoff) {
          energy += prefactor * ek2_forP[idx]
                    * cos(kx[lx]*r[0] + ky[ly]*r[1] + kz_forP[lz]*r[2]);
        }
      }
    }
  }

  return energy;

}
double PotentialEwaldCoul::SelfEnergy(Bead& bead) {
  double q = bead.Charge();
  return -lB*sqrt(alpha/kPi)*q*q;

}

// [[[Note]]]: This function should only be used on the point charges on wall z=0!!!
double PotentialEwaldCoul::PairForceZReal(Bead& bead1, Bead& bead2, int npbc) {
  double force_z = 0;
  double dist[3];
  GetDistVector(bead2, bead1, box_l, npbc, dist);
  double q1 = bead1.Charge();
  double q2 = bead2.Charge();
  double prefactor = lB*q1*q2;
  double prefactor2 = 2*sqrt(alpha/kPi);

  for (int i = -real_cell[0]; i <= real_cell[0]; i++) {
    for (int j = -real_cell[1]; j <= real_cell[1]; j++) {
      for (int k = -real_cell[2]; k <= real_cell[2]; k++) {
        double r_vec[3];
        r_vec[0] = dist[0] + i*box_l[0];
        r_vec[1] = dist[1] + j*box_l[1];
        r_vec[2] = dist[2] + k*box_l[2];
        double r = sqrt(r_vec[0]*r_vec[0]+r_vec[1]*r_vec[1]+r_vec[2]*r_vec[2]);
        if (r > 0 && r <= real_cutoff) {
          force_z += prefactor * ((prefactor2*exp(-alpha*r*r))
                     + (erfc(sqrt(alpha)*r)/r)) * r_vec[2]/(r*r);
        }
      }
    }
  }

  return force_z;

}

// [[[Note]]]: This function should only be used on the point charges on wall z=0!!!
double PotentialEwaldCoul::PairForceZRepl(Bead& bead1, Bead& bead2, int npbc) {
  double force_z = 0;
  double r[3];
  GetDistVector(bead2, bead1, box_l, npbc, r);
  double q1 = bead1.Charge();
  double q2 = bead2.Charge();
  double prefactor = lB * 4*kPi * q1*q2/box_vol;

  for (int lx = 0; lx < repl_ceto[0]; lx++) {
    for (int ly = 0; ly < repl_ceto[1]; ly++) {
      for (int lz = 0; lz < repl_ceto[2]; lz++) {
        int idx = repl_ceto[1]*repl_ceto[2]*lx + repl_ceto[2]*ly + lz;
        if (k2[idx] > 0 && k2[idx] <= repl_cutoff) {
          force_z += prefactor * kz[lz]*ek2[idx]  
                     * sin(kx[lx]*r[0] + ky[ly]*r[1] + kz[lz]*r[2]);
        }
      }
    }
  }

  return force_z;

}

// [[[Note]]]: This function should only be used on the point charges on wall z=0!!!
double PotentialEwaldCoul::ForceZDipole(Bead& bead, double dipole_z) {
  return -lB * 4*kPi*bead.Charge()*dipole_z/box_vol;

}

// For pressure calculation.
double PotentialEwaldCoul::PairDForceReal(Bead& bead1, Bead& bead2,
                                          Bead& bead1_ref, Bead& bead2_ref,
                                          int npbc) {
  // d dot force.
  double dforce = 0;
  double dist[3];
  double dist1r[3];
  double dist2r[3];
  double d[3];
  GetDistVector(bead2, bead1, box_l, npbc, dist);
  GetDistVector(bead1_ref, bead1, box_l, 0, dist1r);
  GetDistVector(bead2_ref, bead2, box_l, 0, dist2r);
  for (int i = 0; i < 3; i++) {
    d[i] = dist1r[i] - dist2r[i];
  }

  double q1 = bead1.Charge();
  double q2 = bead2.Charge();
  double prefactor = lB*q1*q2;

  for (int i = -real_cell[0]; i <= real_cell[0]; i++) {
    for (int j = -real_cell[1]; j <= real_cell[1]; j++) {
      for (int k = -real_cell[2]; k <= real_cell[2]; k++) {
        double r_vec[3];
        r_vec[0] = dist[0] + i*box_l[0];
        r_vec[1] = dist[1] + j*box_l[1];
        r_vec[2] = dist[2] + k*box_l[2];
        double r = sqrt(r_vec[0]*r_vec[0] +
                        r_vec[1]*r_vec[1] +
                        r_vec[2]*r_vec[2]);
        if (r > 0 && r <= real_cutoff) {
          // d * r over r^3.
          //double dror3 = (d[0]*r_vec[0] + d[1]*r_vec[1] + d[2]*r_vec[2])/(r*r*r);
          double dror3 = (dist[0]*r_vec[0]+dist[1]*r_vec[1]+dist[2]*r_vec[2])/(r*r*r);
          dforce += -prefactor * dror3 * (erfc(sqrt(alpha)*r) +
                                          2*sqrt(alpha/kPi)*r*exp(-alpha*r*r));
        }
      }
    }
  }

  return dforce;

}

double PotentialEwaldCoul::PairDForceRepl(Bead& bead1, Bead& bead2,
                                          Bead& bead1_ref, Bead& bead2_ref,
                                          int npbc) {
  double dforce = 0;
  double r[3];
  double dist1r[3];
  double dist2r[3];
  double d[3];
  GetDistVector(bead2, bead1, box_l, npbc, r);
  GetDistVector(bead1_ref, bead1, box_l, 0, dist1r);
  GetDistVector(bead2_ref, bead2, box_l, 0, dist2r);
  for (int i = 0; i < 3; i++) {
    d[i] = dist1r[i] - dist2r[i];
  }
  double q1 = bead1.Charge();
  double q2 = bead2.Charge();
  double prefactor = lB * q1*q2/(kPi*box_vol) * (4*kPi*kPi);

  for (int lx = 0; lx < repl_ceto[0]; lx++) {
    for (int ly = 0; ly < repl_ceto[1]; ly++) {
      for (int lz = 0; lz < repl_ceto[2]; lz++) {
        int idx = repl_ceto[1]*repl_ceto[2]*lx + repl_ceto[2]*ly + lz;
        if (k2[idx] > 0 && k2[idx] <= repl_cutoff) {
          //double dk = d[0]*kx[lx] + d[1]*ky[ly] + d[2]*kz[lz];
          double dk = r[0]*kz[lx] + r[1]*kz[ly] +  r[2]*kz[lz];
          dforce += -dk * prefactor * ek2[idx]
                    * sin(kx[lx]*r[0] + ky[ly]*r[1] + kz[lz]*r[2]);
        }
      }
    }
  }

  return dforce;

}


double PotentialEwaldCoul::DipoleE(vector<Molecule>& mols) {
  double Mz = 0;
  for (int i = 0; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      // If the simulation is done correctly, the z coordinates should all be
      // within the unit cell.
      double z = mols[i].bds[j].GetCrd(0, 2);
      double c = mols[i].bds[j].Charge();
      Mz += c*z;
    }
  }
  return lB * 2*kPi/box_vol * Mz * Mz;

}

double PotentialEwaldCoul::DipoleE(vector<Molecule>& mols,
                                   vector<Bead>& chain) {
  double Mz = 0;
  for (int i = 0; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      // If the simulation is done correctly, the z coordinates should all be
      // within the unit cell.
      double z = mols[i].bds[j].GetCrd(0, 2);
      double c = mols[i].bds[j].Charge();
      Mz += c*z;
    }
  }
  for (int i = 0; i < (int)chain.size(); i++) {
    double z = chain[i].GetCrd(0, 2);
    double c = chain[i].Charge();
    Mz += c*z;
  }

  return lB * 2*kPi/box_vol * Mz * Mz;

}

// The dipole energy without the deleted molecule.
// Here, the delete_id should be >= 0.
double PotentialEwaldCoul::DipoleE(vector<Molecule>& mols, int delete_id,
                                   int counterion) {
  double Mz = 0;
  for (int i = 0; i < (int)mols.size(); i++) {
    if (i < delete_id || i > delete_id+counterion) {
      for (int j = 0; j < mols[i].Size(); j++) {
        // If the simulation is done correctly, the z coordinates should all be
        // within the unit cell.
        double z = mols[i].bds[j].GetCrd(0, 2);
        double c = mols[i].bds[j].Charge();
        Mz += c*z;
      }
    }
  }
  return lB * 2*kPi/box_vol * Mz * Mz;

}

// The dipole energy difference before and after the +/- bead pair is added.
double PotentialEwaldCoul::DipoleEDiff(vector<Molecule>& mols,
                                       vector<Bead>& cbmc_chain,
                                       Bead& bead1, Bead& bead2, 
                                       int current_len, int gc_chain_len,
                                       double gc_bead_charge, int delete_id) {
  double Mz_n = 0;
  double Mz_o = 0;
  int counterion = 0;
  if (gc_bead_charge != 0)  counterion = gc_chain_len;

  // Mols.
  if (delete_id >= 0) {
    for (int i = 0; i < (int)mols.size(); i++) {
      if (i < delete_id || i > delete_id+counterion) {
        for (int j = 0; j < mols[i].Size(); j++) {
          double z = mols[i].bds[j].GetCrd(0, 2);
          double c = mols[i].bds[j].Charge();
          Mz_o += c*z;
        }
      }
    }
  }
  else {
    for (int i = 0; i < (int)mols.size(); i++) {
      for (int j = 0; j < mols[i].Size(); j++) {
        double z = mols[i].bds[j].GetCrd(0, 2);
        double c = mols[i].bds[j].Charge();
        Mz_o += c*z;
      }
    }
  }
  // CBMC chain.
  for (int i = 0; i < current_len; i++) {
    double z = cbmc_chain[i].GetCrd(0, 2);
    double c = cbmc_chain[i].Charge();
    Mz_o += c*z;
  }
  if (gc_bead_charge != 0) {
    for (int i = gc_chain_len; i < gc_chain_len+current_len; i++) {
      double z = cbmc_chain[i].GetCrd(0, 2);
      double c = cbmc_chain[i].Charge();
      Mz_o += c*z;
    }
  }
  // These beads.
  Mz_n = Mz_o;
  double z = bead1.GetCrd(0, 2);
  double c = bead1.Charge();
  Mz_n += c*z;
  if (gc_bead_charge != 0) {
    z = bead2.GetCrd(0, 2);
    c = bead2.Charge();
    Mz_n += c*z;
  }

  return lB * 2*kPi/box_vol * (Mz_n*Mz_n - Mz_o*Mz_o);

}

double PotentialEwaldCoul::GetlB() {
  return lB;

}


