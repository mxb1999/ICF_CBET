#include "include/implSim.h"

using namespace std;
//main function called in program
int main(int argc, char const *argv[]) {
  /*
    Basic steps of the program:
    1. Initialize the Arrays
    2. Launch rays and track for overlapping beams
    3. Perform CBET calculations and update arrays
    4. Update/write HDF5 file with desired output arrays
    5. Plot arrays (performed by matplotting.py)
  */
  auto start1 = chrono::high_resolution_clock::now();
  initialize();
  auto stop1 = chrono::high_resolution_clock::now();
  auto start2 = chrono::high_resolution_clock::now();
  launchRays();
  auto stop2 = chrono::high_resolution_clock::now();
  auto start3 = chrono::high_resolution_clock::now();
  cbet();
  auto stop3 = chrono::high_resolution_clock::now();
  auto start4 = chrono::high_resolution_clock::now();
  updateH5();
  auto stop4 = chrono::high_resolution_clock::now();
  cout << "Initialize CPU Time: " << chrono::duration_cast<chrono::milliseconds>(stop1-start1).count() << " ms" << endl;
  cout << "Ray Launch CPU Time: " << chrono::duration_cast<chrono::milliseconds>(stop2-start2).count() << " ms" << endl;
  cout << "CBET CPU Time: " << chrono::duration_cast<chrono::milliseconds>(stop3-start3).count() << " ms" << endl;
  cout << "HDF5 Write CPU Time: " << chrono::duration_cast<chrono::milliseconds>(stop4-start4).count() << " ms" << endl;
  cout << "_____________________________________________" << endl;
  cout << "Total CPU Time: " << chrono::duration_cast<chrono::milliseconds>(stop4-start1).count() << " ms" << endl;
  return 0;
}
