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
  double divideby = 5.0;  // Decides bin resolution.
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

}

void ForceField::CalcPressureVirialHSELSlit(vector<Molecule>& mols,
                                            double rho) {
  if (use_gc)  UpdateMolCounts(mols);

  vp_z++;

  for (int i = 0; i < n_mol; i++) {
    int i_len = mols[i].Size();
    for (int j = i+1; j < n_mol; j++) {
      int j_len = mols[j].Size();
      for (int k = 0; k < i_len; k++) {
        for (int l = 0; l < j_len; l++) {
          double z1 = mols[i].bds[k].GetCrd(0, 2);
          double z2 = mols[j].bds[l].GetCrd(0, 2);
          double vx  = mols[i].bds[k].BBDistVec(mols[j].bds[l], box_l, npbc, 0);
          double vy  = mols[i].bds[k].BBDistVec(mols[j].bds[l], box_l, npbc, 1);
          double a = sqrt(vx*vx + vy*vy);
          int binz1 = floor(z1/vp_slit_z1_res);
          int binz2 = floor(z2/vp_slit_z2_res);
          int bina = floor(a/vp_slit_a_res);
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
  
            double vz  = mols[i].bds[k].BBDistVec(mols[j].bds[l], box_l,npbc,2);
            double vcx = mols[i].bds[0].BBDistVecWithRef(mols[j].bds[0],
                                mols[i].bds[k], mols[j].bds[l], box_l, npbc, 0);
            double vcy = mols[i].bds[0].BBDistVecWithRef(mols[j].bds[0],
                                mols[i].bds[k], mols[j].bds[l], box_l, npbc, 1);
            double vcz = mols[i].bds[0].BBDistVecWithRef(mols[j].bds[0],
                                mols[i].bds[k], mols[j].bds[l], box_l, npbc, 2);
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


