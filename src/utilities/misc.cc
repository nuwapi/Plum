#include "misc.h"

#include <stdlib.h>
#include <cmath>
#include <random> 

#include "../molecules/bead.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

double getDist(Bead* b1, Bead* b2, double box_l[], int npbc) {
  double dist = 0;
  for (int i = 0; i < 3; i++) {
    double di = b1->GetCrd(1, i) - b2->GetCrd(1,i);
    if (i < npbc) {
      // Periodic boundary conditions
      di -= box_l[i] * round(di / box_l[i]);
    }
    dist += (di * di);
  }

  return sqrt(dist);

}

double getDist(Bead& b1, Bead& b2, double box_l[], int npbc) {
  double dist = 0;
  for (int i = 0; i < 3; i++) {
    double di = b1.GetCrd(1, i) - b2.GetCrd(1,i);
    if (i < npbc) {
      // Periodic boundary conditions
      di -= box_l[i] * round(di / box_l[i]);
    }
    dist += (di * di);
  }

  return sqrt(dist);

}

void GetDistVector(Bead& b1, Bead& b2, double box_l[], int npbc,
                   double (&dist)[3]) {
  // Vector pointing from Bead1 to Bead2.
  for (int i = 0; i < 3; i++) {
    double di = b2.GetCrd(1, i) - b1.GetCrd(1, i);

    // Periodic boundary conditions.
    if (i < npbc) {
      di -= box_l[i] * round(di / box_l[i]);
    }
    // Finding the shortest distance between periodic images.
    if (abs(di) > box_l[i]/2.0) {
      if (di < 0) {
        di += box_l[i];
      }
      else {
        di -= box_l[i];
      }
    }
    dist[i] = di;
  }

}

void GetDistVectorC(Bead& b1, Bead& b2, double box_l[], int npbc,
                    double (&dist)[3]) {
  // Vector pointing from Bead1 to Bead2.
  for (int i = 0; i < 3; i++) {
    double di = b2.GetCrd(0, i) - b1.GetCrd(0, i);

    // Periodic boundary conditions.
    if (i < npbc) {
      di -= box_l[i] * round(di / box_l[i]);
    }
    // Finding the shortest distance between periodic images.
    if (abs(di) > box_l[i]/2.0) {
      if (di < 0) {
        di += box_l[i];
      }
      else {
        di -= box_l[i];
      }
    }
    dist[i] = di;
  }

}

void randSphere(double vec[], mt19937& rand_gen) {
  double rand_square = 2;
  double r1, r2;
  while (rand_square > 1) {
    // Random numbers between -1 and 1.
    r1 = 1 - 2 * ((double)rand_gen() / rand_gen.max()); 
    r2 = 1 - 2 * ((double)rand_gen() / rand_gen.max());
    rand_square = r1 * r1 + r2 * r2;
  }
  double ranh = 2 * sqrt(1 - rand_square);
  vec[0] = r1 * ranh; 
  vec[1] = r2 * ranh; 
  vec[2] = (1 - 2*rand_square); 

}

double gasdev(double mean, double stdev, mt19937& ranGen) {
  static int iset = 0;
  static double gset;
  double fac, rsq, v1, v2;
  if (iset == 0) {
    do {
      v1=2.0*((double)ranGen() / ranGen.max()) - 1.0;
      v2=2.0*((double)ranGen() / ranGen.max()) - 1.0;
      rsq = v1*v1 + v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0 );
      fac = sqrt(-2.0 * log(rsq)/rsq);

      gset = v1*fac;
      iset = 1;
      return v2*fac * stdev + mean;
  }
  else {
    iset = 0;
    return gset * stdev + mean;
  }

  return 0; 

}

string YesOrNo(bool input) {
  if (input) {
    return "yes";
  }
  else {
    return "no";
  }

}

int factorial(int n) {
  if (n != 1) {
     return n*factorial(n-1);
  }
  else {
    return 1;
  }

}

/** This is a simple least square fit procedure. */
double Interpolate(double * bin_pos, double * f, int len, double sigma) {
  /*                     INDEX 2
     X = | bin1_pos 1 |  B = | b1 |  Y = | g(bin1_pos) | d
         | bin2_pos 1 |      | b2 |      | g(bin2_pos) | i     INDEX 1
         |    ...     |                  |     ...     | m e
   */
  MatrixXf X(len, 2);
  VectorXf Y(len);
  for (int i = 0; i < len; i++) {
    X(i, 0) = bin_pos[i];
    X(i, 1) = 1.0;
    Y(i) = f[i];
  }
  VectorXf B = (X.transpose() * X).ldlt().solve(X.transpose() * Y);

  /*
  cout << sigma << '\t' << sigma * B(0) + B(1) << endl;
  for (int i = 0; i < len; i++)
  cout << bin_pos[i] << '\t' << f[i] << endl;
  cout << endl;
  */

  return sigma * B(0) + B(1);

}

double Interpolate2(double * bin_pos, double * f, int len, double sigma) {
  MatrixXf X(len, 2);
  VectorXf Y(len);
  for (int i = 0; i < len; i++) {
    X(i, 0) = bin_pos[i];
    X(i, 1) = 1.0;
    Y(i) = f[i];
  }
  VectorXf B = (X.transpose() * X).ldlt().solve(X.transpose() * Y);
  return bin_pos[len-1] * B(0) + B(1);

}


