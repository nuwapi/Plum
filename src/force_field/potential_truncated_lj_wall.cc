#include "potential_truncated_lj_wall.h" 

#include "../utilities/constants.h"

using namespace std; 

PotentialTruncatedLJWall::PotentialTruncatedLJWall(string potential_name)
                                      : PotentialExternal(potential_name) {
  ReadParameters(); 

}

void PotentialTruncatedLJWall::ReadParameters() {
  cout << setw(35) << "[eP] External potential type: " << "Truncated LJ wall"
       << endl;

  string flag; 
  string symbol; 
  double sigma, epsilon; 
  cin >> flag >> m_cut 
      >> flag >> m_sigWall 
      >> flag >> m_epWall;
  cout << setw(35) << "[eP] LJ cutoff              : " << m_cut << endl;
  cout << setw(35) << "[eP] Sigma for wall         : " << m_sigWall << endl;
  cout << setw(35) << "[eP] Epsilon for wall       : " << m_epWall << endl;

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
      cout << setw(35) << "[eP] Sigma for bead type    : " << symbol << " - "
           << sigma << endl;
      cout << setw(35) << "[eP] Epsilon for bead type  : " << symbol << " - "
           << epsilon << endl;
    }
  }

}

double PotentialTruncatedLJWall::BeadEnergy(Bead& bead, double box_l[]) {
  double energy = 0;
  // Get the z coordinate from the trial coordinate array.
  double z = bead.GetCrd(1, 2);
  double sigma = sigmas[bead.Symbol()];
  double epsilon = epsilons[bead.Symbol()];
  double R0 = 3*k213*sigma;  // For FENE potential. A bead can only be one bead
                             // away from the surface.
  double K = 1;              // For FENE potential.
                             // http://lammps.sandia.gov/doc/bond_fene.html

  if (epsilon == 0) {
    energy = 0;
  }
  // If the bead falls outside of the box in the z-direction.
  else if (z <= 0 || z >= box_l[2]) {
    return kVeryLargeEnergy;
  }
  else {
    // Use purely repulsive wall for non-grafted beads.
    if (m_cut < 0 && bead.Symbol() != "R" && bead.Symbol() != "L") {
      double r3_ref = pow((1.0/k213), 3);
      double energy_ref = 2.59807621135 * epsilon * (r3_ref*r3_ref - r3_ref);
      // 1.20093695518 = 3^(1/6)
      if (z < k213*sigma) {
        double r3 = pow((sigma/z), 3);
        // 2.59807621135 = 3*3^0.5/2
        energy += 2.59807621135 * epsilon * (r3*r3 - r3) - energy_ref;
      }

      if (box_l[2] - z < k213*sigma) {
        double r3 = pow((sigma/(box_l[2] - z)), 3);
        energy += 2.59807621135 * epsilon * (r3*r3 - r3) - energy_ref;
      }
    }
    // For left (z=0) grafted beads.
    else if (bead.Symbol() == "L") {
      // The left attractive well.
      double r3 = pow((sigma/z), 3);
      energy += 2.59807621135 * epsilon * (r3*r3 - r3);
      energy += -0.5*K*R0*R0*log(1-pow(z/R0,2));  // FENE part.

      // The right purely repulsive wall.
      if (box_l[2] - z < k213*sigma) {
        double r3 = pow((sigma/(box_l[2] - z)), 3);
        energy += 2.59807621135 * epsilon * (r3*r3 - r3);
      }
      else {
        double r3 = pow((1.0/k213), 3);
        energy += 2.59807621135 * epsilon * (r3*r3 - r3);
      }
    }
    // For right (z=box_l[2]) grafted beads.
    else if (bead.Symbol() == "R") {
      // The right attractive well.
      double r3 = pow((sigma/(box_l[2] - z)), 3);
      energy +=  2.59807621135 * epsilon * (r3*r3 - r3);
      energy += -0.5*K*R0*R0*log(1-pow((box_l[2]-z)/R0,2));  // FENE part.

      // The left purely repulsive wall.
      if (z < k213*sigma) {
        double r3 = pow((sigma/z), 3);
        // 2.59807621135 = 3*3^0.5/2
        energy += 2.59807621135 * epsilon * (r3*r3 - r3);
      }
      else {
        double r3 = pow((1.0/k213), 3);
        energy += 2.59807621135 * epsilon * (r3*r3 - r3);
      }
    }
    // Full LJ potential for non-grafted beads.
    else {
      double r3_ref = pow((1.0/m_cut), 3);
      double energy_ref = 2.59807621135 * epsilon * (r3_ref*r3_ref - r3_ref);
      if (z < m_cut) {
        double r3 = pow((sigma/z), 3);
        energy +=  2.59807621135 * epsilon * (r3*r3 - r3) - energy_ref;
      }

      if (box_l[2] - z < m_cut) {
        double r3 = pow((sigma/(box_l[2] - z)), 3);
        energy +=  2.59807621135 * epsilon * (r3*r3 - r3) - energy_ref;
      }
    }
  }

  return energy;

}

double PotentialTruncatedLJWall::BeadForceOnWall(Bead& bead, double box_l[]) {
  double force = 0;
  // Get the z coordinate from the trial coordinate array.
  double z = bead.GetCrd(1, 2);
  double sigma = sigmas[bead.Symbol()];
  double epsilon = epsilons[bead.Symbol()];
  double R0 = 3*k213*sigma;
  double K = 1;

  if (epsilon == 0) {
    force = 0;
  }
  // If the bead falls outside of the box in the z-direction.
  else if (z <= 0 || z >= box_l[2]) {
    return kVeryLargeEnergy;
  }
  else {
    // Use purely repulsive wall for non-grafted beads.
    if (m_cut < 0 && bead.Symbol() != "R" && bead.Symbol() != "L") {
      if (z < k213*sigma) {
        double r3 = pow((sigma/z), 3);
        force += 2.59807621135 * epsilon * (r3*r3*6/z - r3*3/z);
      }
      if (box_l[2] - z < k213*sigma) {
        double r3 = pow((sigma/(box_l[2] - z)), 3);
        force += 2.59807621135 * epsilon * (r3*r3*6/(box_l[2]-z) -
                 r3*3/(box_l[2]-z));
      }
    }
    // For left (z=0) grafted beads.
    else if (bead.Symbol() == "L") {
      // The left attractive well.
      double r3 = pow((sigma/z), 3);
      force += 2.59807621135 * epsilon * (r3*r3*6/z - r3*3/z);
      force += - K*R0*z/(1-pow(z/R0, 2));  // FENE part.

      // The right purely repulsive wall.
      if (box_l[2] - z < k213*sigma) {
        double r3 = pow((sigma/(box_l[2] - z)), 3);
        force += 2.59807621135 * epsilon * (r3*r3*6/(box_l[2]-z) -
                 r3*3/(box_l[2]-z));
      }
    }
    // For right (z=box_l[2]) grafted beads.
    else if (bead.Symbol() == "R") {
      // The right attractive well.
      double r3 = pow((sigma/(box_l[2] - z)), 3);
      force += 2.59807621135 * epsilon * (r3*r3*6/(box_l[2]-z) -
               r3*3/(box_l[2]-z));
      force += - K*R0*(box_l[2]-z)/(1-pow((box_l[2]-z)/R0, 2));  // FENE part.

      // The left purely repulsive wall.
      if (z < k213*sigma) {
        double r3 = pow((sigma/z), 3);
        force += 2.59807621135 * epsilon * (r3*r3*6/z - r3*3/z);
      }
    }
    // Full LJ potential for non-grafted beads.
    else {
      if (z < m_cut) {
        double r3 = pow((sigma/z), 3);
        force += 2.59807621135 * epsilon * (r3*r3*6/z - r3*3/z);
      }
      else {
        double r3 = pow((sigma/m_cut), 3);
        force += 2.59807621135 * epsilon * (r3*r3*6/z - r3*3/z);
      }
      if (box_l[2] - z < m_cut) {
        double r3 = pow((sigma/(box_l[2] - z)), 3);
        force += 2.59807621135 * epsilon * (r3*r3*6/(box_l[2]-z) - 
                 r3*3/(box_l[2]-z));
      }
      else {
        double r3 = pow((sigma/m_cut), 3);
        force += 2.59807621135 * epsilon * (r3*r3*6/z - r3*3/z);
      }
    }
  }

  return force;

}


