#include "implSim.hpp"
#include <Python.h>
using namespace std;

void benchmark()
{
  int trials = 10;
  /*calcCBET = 0;
  for(int t = 0; t < trials;t++)
  {
    for(int tr = 1; tr <= 12; tr++)
    {
      threads = tr;
      for(int i = 2; i < 100; i+=2)
      {
        rays_per_zone = i;
        uray_mult = intensity*(courant_mult)*pow(double(rays_per_zone),-1.0); //multiplier which determines intensity deposited in a given zone
        nrays= int(rays_per_zone*(beam_max_z-beam_min_z)/dz)+0;//number of rays per beam
        launchRays();
        cbet();
        freeCBETArrs();
        freeTraceArrs();
      }
    }
  }*/
  optimize = 0;
  cudaCalc = 1;
  
  for(int i = 2; i < 100; i+=2)
  {
    rays_per_zone = i;
    for(int t = 0; t < trials;t++)
    {
      uray_mult = intensity*(courant_mult)*pow(double(rays_per_zone),-1.0); //multiplier which determines intensity deposited in a given zone
      nrays= int(rays_per_zone*(beam_max_z-beam_min_z)/dz)+0;//number of rays per beam
      launchRays();
      cbet();
      freeTraceArrs();
      freeCBETArrs();
    }
  }
}
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
 
  threads = 12;
  if(argc > 2)
  {
     pyPlot = argv[1][0] - 48;
     std::string str(argv[2]);
     threads = std::stoi(str);
  }else if(argc > 1)
  {
    pyPlot = argv[1][0] - 48;
  }else
  {
    pyPlot = 0;
  }
  printf("Threads: %d\n", threads);
  auto start1 = chrono::high_resolution_clock::now();
  initialize();
  optimize = 0;
  int b = 0;
  if(b)
  {
    benchmark();
    return 0;
  }
  auto stop1 = chrono::high_resolution_clock::now();
  auto start2 = chrono::high_resolution_clock::now();
  launchRays();
  auto stop2 = chrono::high_resolution_clock::now();
  auto start3 = chrono::high_resolution_clock::now();
  optimize = 0;
  if(calcCBET)
  {
    cbet();
  }

  auto stop3 = chrono::high_resolution_clock::now();
  auto start4 = chrono::high_resolution_clock::now();
  updateH5();

  auto stop4 = chrono::high_resolution_clock::now();
  if(printTimings)
  {
    cout << "Initialize CPU Time: " << chrono::duration_cast<chrono::milliseconds>(stop1-start1).count() << " ms" << endl;
    cout << "Ray Launch CPU Time: " << chrono::duration_cast<chrono::milliseconds>(stop2-start2).count() << " ms" << endl;
    cout << "CBET CPU Time: " << chrono::duration_cast<chrono::milliseconds>(stop3-start3).count() << " ms" << endl;
    cout << "HDF5 Write CPU Time: " << chrono::duration_cast<chrono::milliseconds>(stop4-start4).count() << " ms" << endl;
    cout << "_____________________________________________" << endl;
    cout << "Total CPU Time: " << chrono::duration_cast<chrono::milliseconds>(stop4-start1).count() << " ms" << endl;
  }
  if(pyPlot)
  {
    char filename[] = "matplotting.py";
    FILE* fp;
    Py_Initialize();
    cout << "Hello" << endl;
    fp = _Py_fopen(filename, "r");
    PyRun_SimpleFile(fp, filename);
    Py_Finalize();
  }
  return 0;
}
