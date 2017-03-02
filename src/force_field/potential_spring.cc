#include "potential_spring.h" 

using namespace std; 

PotentialSpring::PotentialSpring(int n_mol, string potential_name)
                            : PotentialBond(n_mol, potential_name) {
  ReadParameters(); 

}

void PotentialSpring::ReadParameters() {
  cout << setw(35) << "Bond potential type         : "
       << "Harmonic spring" << endl;

  string flag; 
  cin >> flag >> m_kBond >> flag >> m_r0; 
  cout << setw(37) << "Bond strength in eV/ul^2    = " << m_kBond << endl;
  cout << setw(37) << "Resting bond length in ul   = " << m_r0 << endl;

}

double PotentialSpring::MoleculeEnergy(Molecule& mol, double box_l[],
                                       int npbc) {
  double bond_energy = 0; 
  for (int i = 0; i < mol.NBond(); i++) {
    // Indices of the two beads forming the bond.
    int ind1 = mol.bonds[i][0];
    int ind2 = mol.bonds[i][1];
    // Current bond length.
    // Setting npbc to 0 meaning that we do not wrap bonds.
    double r = mol.bds[ind1].BBDist(mol.bds[ind2], box_l, 0);
    bond_energy += 0.5 * m_kBond * (r - m_r0) * (r - m_r0); 
  }
  return bond_energy;

}

// Returns a bond length from appropriate distribution for harmonic bond.
double PotentialSpring::RandomBondLen(double beta, mt19937& ranGen) {
  double len = 0; 
  double sigma = sqrt(1/(beta * m_kBond));
  double a = (m_r0 + 3*sigma) * (m_r0 + 3*sigma);
  bool ready = false;
  while (!ready) {
    len = gasdev(m_r0, sigma, ranGen); 
    ready = ((double)ranGen()/ranGen.max() < (len * len/a)); 
  }
  return len; 

}

double PotentialSpring::EnergyDifference(vector<Molecule>& mols,
                                         double box_l[], int npbc,
                                         int mol_id) {
  // If a the bead in a chain is translated.
  if (mols[mol_id].Size() > 1) {
    double trial_energy = MoleculeEnergy(mols[mol_id], box_l, npbc); 
    SetE(1, mol_id, trial_energy);
    dE = (trial_energy - GetE(0, mol_id)); 
    return dE;
  }
  // If it is a single bead molecule that is translated.
  else {
    return 0; 
  }

  return 0;

}


