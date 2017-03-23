#include <iostream>

#include "simulation/simulation.h"

using namespace std; 

int main(int argc, char * argv[]) {
  cout << " " << endl;
  cout << "     _ _       " << endl;
  cout << "   /. ' \\  -----=>> Plum v1.0 <<=-----" << endl;
  cout << "  | ~^ ^~| ---=>> Nuo Wang 2017 <<=---" << endl;
  cout << "   \\__\\_/  \n" << endl;
  cout << "  All energies are in kBT." << endl;
  cout << "  All lengths are in unit length (ul)." << endl;
  cout << "\n  Starting program!" << endl;

  // Initialize the simulation.
  Simulation simulation;
  // Run the simulation.
  simulation.Run();

  cout << "  END." << endl;

  return 0;

}


