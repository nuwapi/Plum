#include "force_field.h"

#include "../utilities/constants.h"

void ForceField::InitPressureVirialHSELSlit() {
  if (!use_bond_rigid) {
    cout << "  Rigid bond is needed for pressure calculation!" << endl;
    cout << "  Exiting! Program complete." << endl;
    exit(1);
  }
  if (box_l[0] != box_l[1]) {
    cout << "  It is required that the x and y dimensions of the" << endl;
    cout << "  box are equal for slit geometry! Exiting! Program" << endl;
    cout << "  complete." << endl;
    exit(1);
  }

  vp_z = 0;
  vp_slit_margin = rigid_bond/1000;
  vp_slit_vol = box_l[2] * kPi * (box_l[0]/2.0) * (box_l[1]/2.0);
  vp_slit_Nm = chain_len;
  double divideby = 1.0;  // Decides bin resolution.
  vp_slit_z1_bin = floor(box_l[2]/(rigid_bond/divideby));
  vp_slit_z2_bin = vp_slit_z1_bin;
  vp_slit_a_bin = floor(0.5*box_l[0]/(rigid_bond/divideby));
  vp_slit_z1_res = box_l[2] / vp_slit_z1_bin;
  vp_slit_z2_res = vp_slit_z1_res;
  vp_slit_a_res = box_l[0] / vp_slit_a_bin;
  vp_d1 = vp_slit_Nm + 2;
  vp_d2 = vp_slit_Nm + 2;
  vp_d3 = vp_slit_z1_bin;
  vp_d4 = vp_slit_z2_bin;
  vp_d5 = vp_slit_a_bin;
  vp_d6 = 3;

  int array_size = vp_d1 * vp_d2 * vp_d3 * vp_d4 * vp_d5 * vp_d6;
  vp_slit_rho = new double[1];  // currently not in use.
  vp_slit_el = new double[array_size];
  array_size = vp_d1 * vp_d2 * vp_d3 * vp_d4 * vp_d6;
  vp_slit_hs = new double[array_size];
  lB = 0;
  if (use_ewald_pot)  ewald_pot->GetlB();

}

void ForceField::CalcPressureVirialHSELSlit(vector<Molecule>& mols,
                                            double rho) {
  if (use_gc)  UpdateMolCounts(mols);
  int npbc_h = 2;  // The name means "npbc here (in this function)".

  vp_z++;

  for (int i = 0; i < n_mol; i++) {
    int i_len = mols[i].Size();
    for (int j = i+1; j < n_mol; j++) {
      int j_len = mols[j].Size();
      for (int k = 0; k < i_len; k++) {
        for (int l = 0; l < j_len; l++) {
          double z1 = mols[i].bds[k].GetCrd(0, 2);
          double z2 = mols[j].bds[l].GetCrd(0, 2);
          double vx  = mols[i].bds[k].BBDistVec(mols[j].bds[l], box_l, npbc_h, 0);
          double vy  = mols[i].bds[k].BBDistVec(mols[j].bds[l], box_l, npbc_h, 1);
          double a = sqrt(vx*vx + vy*vy);
          int binz1 = floor(z1/vp_slit_z1_res);
          int binz2 = floor(z2/vp_slit_z2_res);
          int bina = floor(a/vp_slit_a_res);
          // The next four IFs need to be integrated into the fifth IF because
          // no z should be out of bound, the fifty IF is not needed.
          if (binz1 < 0)                binz1 = 0;
          if (binz1 >= vp_slit_z1_bin)  binz1 = vp_slit_z1_bin - 1;
          if (binz2 < 0)                binz2 = 0;
          if (binz2 >= vp_slit_z2_bin)  binz2 = vp_slit_z2_bin - 1;
          if (binz1 >= 0 && binz1 < vp_slit_z1_bin &&
              binz2 >= 0 && binz2 < vp_slit_z2_bin &&
              a < vp_slit_a_bin) {
            int s1, s2;
            double c1 = mols[i].bds[k].Charge();
            double c2 = mols[j].bds[l].Charge();
            if (i_len > 1)     s1 = k;
            else if (c1 >= 0)  s1 = chain_len;
            else               s1 = chain_len + 1;
            if (j_len > 1)     s2 = l;
            else if (c2 >= 0)  s2 = chain_len;
            else               s2 = chain_len + 1;
  
            double vz  = mols[i].bds[k].BBDistVec(mols[j].bds[l], box_l,npbc_h,2);
            double vcx = mols[i].bds[0].BBDistVecWithRef(mols[j].bds[0],
                                mols[i].bds[k], mols[j].bds[l], box_l, npbc_h, 0);
            double vcy = mols[i].bds[0].BBDistVecWithRef(mols[j].bds[0],
                                mols[i].bds[k], mols[j].bds[l], box_l, npbc_h, 1);
            double vcz = mols[i].bds[0].BBDistVecWithRef(mols[j].bds[0],
                                mols[i].bds[k], mols[j].bds[l], box_l, npbc_h, 2);
            double vlen = sqrt(vx*vx + vy*vy + vz*vz);
            double hs_xx = vcx * vx / beta;
            double hs_yy = vcy * vy / beta;
            double hs_zz = vcz * vz / beta;
            double el_xx = vcx * vx / pow(vlen, 3) * lB * c1 * c2;
            double el_yy = vcy * vy / pow(vlen, 3) * lB * c1 * c2;
            double el_zz = vcz * vz / pow(vlen, 3) * lB * c1 * c2;

            if (vlen < rigid_bond + vp_slit_margin) {
              int ind_x_hs = (((s1*vp_d2 + s2)*vp_d3 + binz1)*vp_d4 + binz2)*vp_d6 + 0;
              int ind_y_hs = (((s1*vp_d2 + s2)*vp_d3 + binz1)*vp_d4 + binz2)*vp_d6 + 1;
              int ind_z_hs = (((s1*vp_d2 + s2)*vp_d3 + binz1)*vp_d4 + binz2)*vp_d6 + 2;
              vp_slit_hs[ind_x_hs] += hs_xx;
              vp_slit_hs[ind_y_hs] += hs_yy;
              vp_slit_hs[ind_z_hs] += hs_zz;
            }
            if (use_ewald_pot) {
              int ind_x_el =
              ((((s1*vp_d2 + s2)*vp_d3 + binz1)*vp_d4 + binz2)*vp_d5 + bina)*vp_d6 + 0;
              int ind_y_el =
              ((((s1*vp_d2 + s2)*vp_d3 + binz1)*vp_d4 + binz2)*vp_d5 + bina)*vp_d6 + 1;
              int ind_z_el =
              ((((s1*vp_d2 + s2)*vp_d3 + binz1)*vp_d4 + binz2)*vp_d5 + bina)*vp_d6 + 2;
              vp_slit_el[ind_x_el] += el_xx;
              vp_slit_el[ind_y_el] += el_yy;
              vp_slit_el[ind_z_el] += el_zz;
            }
          }
        }
      }
    }
  }

  p_tensor_hs[0] = p_tensor_hs[1] = p_tensor_hs[2] = 0;
  p_tensor_el[0] = p_tensor_el[1] = p_tensor_el[2] = 0;
  p_tensor[0]    = p_tensor[1]    = p_tensor[2]    = 0;

  for (int d1 = 0; d1 < vp_d1; d1++) {
    for (int d2 = 0; d2 < vp_d2; d2++) {
      for (int d3 = 0; d3 < vp_d3; d3++) {
        for (int d4 = 0; d4 < vp_d4; d4++) {
          int ind_x_hs = (((d1*vp_d2 + d2)*vp_d3 + d3)*vp_d4 + d4)*vp_d6 + 0;
          int ind_y_hs = (((d1*vp_d2 + d2)*vp_d3 + d3)*vp_d4 + d4)*vp_d6 + 1;
          int ind_z_hs = (((d1*vp_d2 + d2)*vp_d3 + d3)*vp_d4 + d4)*vp_d6 + 2;
          p_tensor_hs[0] += vp_slit_hs[ind_x_hs];
          p_tensor_hs[1] += vp_slit_hs[ind_y_hs];
          p_tensor_hs[2] += vp_slit_hs[ind_z_hs]; 

          if (use_ewald_pot) {
            for (int d5 = 0; d5 < vp_d5; d5++) {
              int ind_x_el =
              ((((d1*vp_d2 + d2)*vp_d3 + d3)*vp_d4 + d4)*vp_d5 + d5)*vp_d6 + 0;
              int ind_y_el =
              ((((d1*vp_d2 + d2)*vp_d3 + d3)*vp_d4 + d4)*vp_d5 + d5)*vp_d6 + 1;
              int ind_z_el =
              ((((d1*vp_d2 + d2)*vp_d3 + d3)*vp_d4 + d4)*vp_d5 + d5)*vp_d6 + 2;
              p_tensor_el[0] += vp_slit_el[ind_x_el];
              p_tensor_el[1] += vp_slit_el[ind_y_el];
              p_tensor_el[2] += vp_slit_el[ind_z_el];
            }
          }

        }
      }
    }
  }

  p_tensor_hs[0] *= 1.0/(vp_slit_vol*vp_z);
  p_tensor_hs[1] *= 1.0/(vp_slit_vol*vp_z);
  p_tensor_hs[2] *= 1.0/(vp_slit_vol*vp_z);
  p_tensor_el[0] *= 1.0/(vp_slit_vol*vp_z);
  p_tensor_el[1] *= 1.0/(vp_slit_vol*vp_z);
  p_tensor_el[2] *= 1.0/(vp_slit_vol*vp_z);

  p_tensor[0] = rho/beta + p_tensor_hs[0] + p_tensor_el[0];
  p_tensor[1] = rho/beta + p_tensor_hs[1] + p_tensor_el[1];
  p_tensor[2] = rho/beta + p_tensor_hs[2] + p_tensor_el[2];

}

void ForceField::CalcPressureForceELSlit(vector<Molecule>& mols) {
  vp_z++;

  /*
  // Assume that the two plates carry the same number of point charges.
  // Then the number of point charges on plate z=0 is:
  int phantom_z0 = phantom/2;
  // The number of point charges to use in the pressure sampling. Anywhere
  // between 1 and phantom_z0.
  int phantom_use = phantom_z0;
  double area = box_l[0]*box_l[1] / (double)phantom_z0;

  // Calculating dipole moment.
  double dipole_z = 0;
  for (int i = 0; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      dipole_z += mols[i].bds[j].Charge() * mols[i].bds[j].GetCrd(0,2);
    }
  }

  for (int i = 0; i < phantom_use; i++) {
    for (int j = phantom_z0; j < (int)mols.size(); j++) {
      for (int k = 0; k < mols[j].Size(); k++) {
        double p_real = 0;
        double p_repl = 0;
        if (use_ewald_pot) {
          p_real = ewald_pot->PairForceZReal(mols[i].bds[0], mols[j].bds[k], npbc) /
                   area;
          p_repl = ewald_pot->PairForceZRepl(mols[i].bds[0], mols[j].bds[k], npbc) /
                   area;
        }
        // Plate 2.
        if (j < phantom)
          p_tensor_el[3] += p_real + p_repl;
        // Polymer.
        else if (mols[j].Size() > 1)
          p_tensor_el[0] += p_real + p_repl;
        // Cation.
        else if (mols[j].bds[0].Charge() > 0)
          p_tensor_el[1] += p_real + p_repl;
        // Anion.
        else if (mols[j].bds[0].Charge() < 0)
          p_tensor_el[2] += p_real + p_repl;
      }
    }
    if (use_ewald_pot) {
      // Dipole correction contribution.
      p_tensor_el[4] += ewald_pot->ForceZDipole(mols[i].bds[0], dipole_z) / area;
    }
  }

  // !!! Alternative mapping:
  >// p_tensor[0] is polymer elec contribution.
  // p_tensor[1] is cation elec contribution.
  // p_tensor[2] is anion elec contribution.
  // p_tensor[3] is plate 2 elec contribution.
  // p_tensor[4] is dipole correction contribution.
  // p_tensor[5] is total.
  p_tensor[0] = p_tensor_el[0] / (vp_z*phantom_use);
  p_tensor[1] = p_tensor_el[1] / (vp_z*phantom_use);
  p_tensor[2] = p_tensor_el[2] / (vp_z*phantom_use);
  p_tensor[3] = p_tensor_el[3] / (vp_z*phantom_use);
  p_tensor[4] = p_tensor_el[4] / (vp_z*phantom_use);
  p_tensor[5] = p_tensor[0] + p_tensor[1] + p_tensor[2] + p_tensor[3] +
                p_tensor[4];
  */

  double total_force = 0;
  double area = box_l[0] * box_l[1];
  for (int i = 0; i < (int)mols.size(); i++) {
    for (int j = 0; j < mols[i].Size(); j++) {
      for (int k = i; k < (int)mols.size(); k++) {
        for (int l = 0; l < mols[k].Size(); l++) {
          if (k > i || (k == i && l >= j)) {
            total_force += ewald_pot->PairForceZReal(mols[i].bds[j], mols[k].bds[l], npbc);
            total_force += ewald_pot->PairForceZRepl(mols[i].bds[j], mols[k].bds[l], npbc);
          }
        }
      }
    }
  }
  p_tensor_el[0] += total_force / area;
  p_tensor[0] = p_tensor_el[0] / vp_z;

}

/*
 * ID codes: 0 - cation, 1 - anion, 2 - polymer, 3 - surface.
 * p_tensor[0]: average Yethiraj total total.
 * p_tensor[1]: average de Miguel.
 * p_tensor[2]: cumulative Yethiraj ideal gas part.
 * p_tensor[3]: cumulative de Miguel pressure.
 * p_tensor[4]: cumulative Yethiraj dipole part.
 * p_tensor_el[i*4+j]: the cumulative Yethiraj electrostatic part of the pair ij.
 * p_tensor_hs[i*4+j]: the cumulative Yethiraj LJ part of the pair ij.
 */
void ForceField::CalcPressureVolScalingHSELSlit(vector<Molecule>& mols) {
  /////////////////////////////////////////
  // 1. Preparation.
  /////////////////////////////////////////
  vp_z++;
  double dU = 0;
  double oldE = 0;
  double newE = 0;
  double box_l_scaled[3] = {box_l[0], box_l[1], box_l[2]+kDz};
  double dz = kDz;
  //if (use_ewald_pot && use_ext_pot)  dz = kDz*kDiCorrection;

  double * d_com_z;
  if (n_mol != 0)
    d_com_z = new double[n_mol];
  else
    d_com_z = new double[1];

  /////////////////////////////////////////
  // 2. Calculate displacements.
  /////////////////////////////////////////
  for (int i = 0; i < n_mol; i++) {
    // Calculate the center of mass position along z axis.
    double com_z = 0;
    for (int j = 0; j < mols[i].Size(); j++) {
      com_z += mols[i].bds[j].GetCrd(0, 2);
    }
    com_z /= mols[i].Size();
    // Find out how much the com should move.
    d_com_z[i] = kDz * (com_z / box_l[2]);
    // For grafted polymers.
    if (mols[i].bds[0].Symbol() == "R")
      d_com_z[i] = kDz;
    else if (mols[i].bds[0].Symbol() == "L")
      d_com_z[i] = 0;
  }
  /*
  // If a counter ion is condensed on a polymer, move it the same amount as the
  // polymer is moved.
  for (int i = phantom; i < n_mol; i++) {
    if (mols[i].Size() == 1 && mols[i].bds[0].Charge() > 0) {
      for (int j = phantom; j < n_mol; j++) {
        if (mols[j].Size() > 1) {
          for (int k = 0; k < mols[j].Size(); k++) {
            double r = mols[i].bds[0].BBDistC(mols[j].bds[k], box_l, npbc);
            if (r < 3.0)  d_com_z[i] = d_com_z[j];
          }
        }
      }
    }
  }
  */

  /////////////////////////////////////////
  // 3. Calculate all dU components.
  /////////////////////////////////////////
  for (int i = 0; i < n_mol; i++) {
    int id_i;
    if (i < phantom)                        id_i = 3;  // Surface.
    else if (mols[i].Size() > 1)            id_i = 2;  // Polymer.
    else if (mols[i].bds[0].Charge() >= 0)  id_i = 0;  // Cation.
    else                                    id_i = 1;  // Anion.

    // For every atom in the molecule, calculate its energy with all atoms that
    // are not in the same molecule.
    for (int j = 0; j < mols[i].Size(); j++) {
      // Update atom trial position.
      mols[i].bds[j].SetCrd(1, 2, mols[i].bds[j].GetCrd(0, 2) + d_com_z[i]);

      // For all other molecules.
      for (int k = i + 1; k < n_mol; k++) {
        for (int l = 0; l < mols[k].Size(); l++) {
          if (i >= phantom || (i < phantom/2 && k >= phantom/2)) {
            mols[k].bds[l].SetCrd(1, 2, mols[k].bds[l].GetCrd(0, 2) + d_com_z[k]);

            int id_k;
            if (k < phantom)                        id_k = 3;  // Surface.
            else if (mols[k].Size() > 1)            id_k = 2;  // Polymer.
            else if (mols[k].bds[0].Charge() >= 0)  id_k = 0;  // Cation.
            else                                    id_k = 1;  // Anion.

            int index1 = (id_i < id_k)? id_i : id_k;  // The smaller index.
            int index2 = (id_i > id_k)? id_i : id_k;  // The bigger index.
            int index = index1 * 4 + index2;
            // Reusing variables.
            index1 = min(mols[i].bds[j].ID(), mols[k].bds[l].ID());
            index2 = max(mols[i].bds[j].ID(), mols[k].bds[l].ID());

            if (use_ewald_pot) {
              oldE = ewald_pot->GetERealRepl(0, index1, index2);
              newE = ewald_pot->PairEnergyRealForP(mols[i].bds[j], mols[k].bds[l], npbc);
              newE += ewald_pot->PairEnergyReplForP(mols[i].bds[j], mols[k].bds[l], npbc);
              dU += newE - oldE;
              p_tensor_el[index] += newE - oldE;
/*
if (id_i == 0 && id_k == 0) {
double nowd = mols[i].bds[j].BBDist(mols[k].bds[l], box_l_scaled, npbc);
double orid = mols[i].bds[j].BBDistC(mols[k].bds[l], box_l_scaled, npbc);
double ko = ewald_pot->GetERepl(0, index1, index2);
double kn = ewald_pot->PairEnergyReplForP(mols[i].bds[j], mols[k].bds[l], npbc);
ewald_pot->PairEnergyRealA(mols[i].bds[j], mols[k].bds[l], npbc);
cout.precision(17);
double dO[3];
double dN[3];
GetDistVectorC(mols[i].bds[j], mols[k].bds[l], box_l_scaled, npbc, dO);
GetDistVector(mols[i].bds[j], mols[k].bds[l], box_l_scaled, npbc, dN);
cout << nowd - orid << " "
     << dO[0] << " " << dO[1] << " " << dO[2] << " "
     << dN[0] << " " << dN[1] << " " << dN[2] << " "
     << kn-ko << endl;
}
*/
            }
            if (use_pair_pot && i >= phantom && k >= phantom) {
              oldE = pair_pot->GetE(0, index1, index2);
              newE = pair_pot->PairEnergy(mols[i].bds[j], mols[k].bds[l], box_l_scaled, npbc);
              dU += newE - oldE;
              p_tensor_hs[index] += newE - oldE;
            }

            mols[k].bds[l].SetCrd(1, 2, mols[k].bds[l].GetCrd(0, 2));
          }
        }
      }

      // Calculate surface-particle LJ interaction.
      if (use_ext_pot && i >= phantom) {
        int index1 = (id_i < 3)? id_i : 3;  // The smaller index.
        int index2 = (id_i > 3)? id_i : 3;  // The bigger index.
        int index = index1 * 4 + index2;
        oldE = ext_pot->GetE(0, mols[i].bds[j].ID());
        newE = ext_pot->BeadEnergy(mols[i].bds[j], box_l_scaled);
        dU += newE - oldE;
        p_tensor_hs[index] += newE - oldE;
      }

      // Recover atom trial position
      mols[i].bds[j].SetCrd(1, 2, mols[i].bds[j].GetCrd(0, 2));
    }
  }

  // Calculate the ideal part of the pressure.
  p_tensor[2] += (n_mol-phantom)/vol;

  // Calculate the dipole part of pressure.
  if (use_ewald_pot && ewald_pot->UseDipoleCorrection()) {
    double Mz_old = 0;
    double Mz_new = 0;
    double lB = ewald_pot->GetlB();
    for (int i = 0; i < (int)mols.size(); i++) {
      for (int j = 0; j < mols[i].Size(); j++) {
        double z = mols[i].bds[j].GetCrd(0, 2);
        double c = mols[i].bds[j].Charge();
        Mz_old += c*z;
        Mz_new += c*(z+d_com_z[i]);
      }
    }
    double d_di = (lB*2*kPi)*(Mz_new*Mz_new/(box_l[0]*box_l[1]*(box_l[2]+dz))
                             -Mz_old*Mz_old/vol);
    dU += d_di;
    p_tensor[4] += d_di;
  }

  // de Miguel.
  p_tensor[3] += pow(1.0+dz/box_l[2], n_mol-phantom)*exp(-beta*dU);

  /////////////////////////////////////////
  // 4. Compile results and calculate pressure.
  /////////////////////////////////////////
  // Yethiraj.
  p_tensor[0] = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = i; j < 4; j++) {
      int index = i * 4 + j;
      p_tensor[0] += p_tensor_el[index] + p_tensor_hs[index];
    }
  }
  p_tensor[0] += p_tensor[4];
  p_tensor[0] /= (box_l[0]*box_l[1]*dz);
  p_tensor[0] = (beta*p_tensor[2] - p_tensor[0])/vp_z;

  // Calculate Yethiraj components of the average.
  for (int i = 0; i < 4; i++) {
    for (int j = i; j < 4; j++) {
      int index = i * 4 + j;
      p_tensor2[index] = -p_tensor_hs[index]/(box_l[0]*box_l[1]*dz*vp_z);
      p_tensor3[index] = -p_tensor_el[index]/(box_l[0]*box_l[1]*dz*vp_z);

    }
  }
  p_tensor2[18] = -p_tensor[4]/(box_l[0]*box_l[1]*dz*vp_z);
  p_tensor2[19] = beta*p_tensor[2]/vp_z; 

  // de Miguel.
  p_tensor[1] = beta/(box_l[0]*box_l[1]*dz)*log(p_tensor[3]/vp_z);

  delete [] d_com_z;
//cout << "a\n" << endl;

}


