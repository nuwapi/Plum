#include "potential_truncated_lj_wall.h" 

#include "../utilities/constants.h"

using namespace std; 

PotentialTruncatedLJWall::PotentialTruncatedLJWall(string potential_name)
                                      : PotentialExternal(potential_name) {
  ReadParameters(); 

}

void PotentialTruncatedLJWall::ReadParameters() {
  cout << setw(35) << "External potential type     : "
       << "Truncated LJ wall" << endl;

  string flag; 
  string symbol; 
  double sigma, epsilon; 
  cin >> flag >> m_cut 
      >> flag >> m_sigWall 
      >> flag >> m_epWall;
  cout << setw(37) << "LJ cutoff                   = " << m_cut << endl;
  cout << setw(37) << "Sigma for wall              = " << m_sigWall << endl;
  cout << setw(37) << "Epsilon for wall            = " << m_epWall << endl;

  while (true) { 
     cin >> flag >> symbol; 
    if (symbol == "end") { 
      break; 
    } 
    else { 
      cin >> flag >> sigma 
          >> flag >> epsilon;  
      sigmas[symbol] = sigma;  
      epsilons[symbol] = epsilon;  
      cout << setw(31) << "Sigma for bead type     " << setw(3) << symbol
           << " = " << sigma << endl;
      cout << setw(31) << "Epsilon for bead type   " << setw(3) << symbol
           << " = " << epsilon << endl;
    } 
  }

}

double PotentialTruncatedLJWall::BeadEnergy(Bead& bead, double box_l[]) {
  double energy = 0;
  // Get the z coordinate from the trial coordinate array.
  double z = bead.GetCrd(1, 2);
  double sigma = sigmas[bead.Symbol()];
  double epsilon = epsilons[bead.Symbol()];

  if (epsilon == 0) {
    energy = 0;
  }
  // If the bead falls outside of the box in the z-direction.
  else if (z <= 0 || z >= box_l[2]) {
    return kVeryLargeEnergy;
  }
  else {
    // Use WCA.
    if (m_cut < 0) {
      // 1.20093695518 = 3^(1/6)
      if (z < sigma) {
        double r3 = pow((sigma/z), 3);
        // 2.59807621135 = 3*3^0.5/2
        energy += 2.59807621135 * epsilon * (r3*r3 - r3);
      }
      if (box_l[2] - z < sigma) {
        double r3 = pow((sigma/(box_l[2] - z)), 3);
        energy += 2.59807621135 * epsilon * (r3*r3 - r3);
      }
    }
    // Use wall L-J.
    else {
      if (z < m_cut) {
        double r3 = pow((sigma/z), 3);
        energy +=  2.59807621135 * epsilon * (r3*r3 - r3);
      }
      if (box_l[2] - z < m_cut) {
        double r3 = pow((sigma/(box_l[2] - z)), 3);
        energy +=  2.59807621135 * epsilon * (r3*r3 - r3);
      }
    }
  }

  return energy;

}

// Calculate force PER UNIT AREA due to molecule-wall interaction.
double PotentialTruncatedLJWall::CalculateForce(vector < Molecule >& mols, double box_l[]){
  double force = 0; 
  for(int i = 0; i < (int)mols.size(); i++){
    for(int j = 0; j < mols[i].Size(); j++){
      double zCoord = mols[i].bds[j].GetCrd(1,2); 
      if(zCoord < 0 || zCoord > box_l[2]){
        cout<<"bead outside the box for force calculation! endng program!"<<endl; 
        exit(1); 
      }
      if(zCoord < m_cut){
        double sig = sigmas[mols[i].bds[j].Symbol()];
        double ep = epsilons[mols[i].bds[j].Symbol()];
        force +=  4 * ep * (3*pow(sig,3)*(1/pow(zCoord,4)) - 9 * pow(sig,9)*(1/pow(zCoord,10)));
      }
      double rRight = box_l[2] - zCoord; //dist from right-hand wall
                  if(rRight < m_cut){
        double sig = sigmas[mols[i].bds[j].Symbol()];
                                double ep = epsilons[mols[i].bds[j].Symbol()];
                                force += 4 * ep * (3*pow(sig,3)*(1/pow(rRight,4)) - 9 * pow(sig,9)*(1/pow(rRight,10)));
      }
    }
  }
  force /= (box_l[0]*box_l[1]*2); //force per unit area 
  force *= (16.02); //convert from messy units to units of 10^10 Pa
  return -force; 

}


