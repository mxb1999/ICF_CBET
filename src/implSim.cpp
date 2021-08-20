#include "implSim.hpp"
#include <Python.h>
#include <fstream>
using namespace std;

void benchmark()
{
  output = new std::ofstream();
  output->open("output.txt", std::ofstream::out);
  int trials = 1;
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
  int numpoints = 14;
  int lim = pow(2,numpoints);
  for(int i = 2; i <= lim; i*=2)
  {
    nrays = i;
    for(int t = 0; t < trials;t++)
    {
      launchRays();
      cbet();
      freeTraceArrs();
      freeCBETArrs();
    }
  }
  cudaCalc = 0;
  threads = 1;
  for(int i = 2; i <= lim; i*=2)
  {
    nrays = i;
    for(int t = 0; t < trials;t++)
    {
      launchRays();
      cbet();
      freeTraceArrs();
      freeCBETArrs();
    }
  }
  cudaCalc = 0;
  threads = 28;
  for(int i = 2; i <= lim; i*=2)
  {
    nrays = i;
    for(int t = 0; t < trials;t++)
    {
      launchRays();
      cbet();
      freeTraceArrs();
      freeCBETArrs();
    }
  }
  
  output->close();
}
void importFromH5(hid_t file, char* set_name, void** target, size_t unitSize, hid_t h5DataType)
{
  hid_t space, dset;
  herr_t status;
  dset = H5Dopen(file, set_name, H5P_DEFAULT);
  space = H5Dget_space(dset);
  const int ndims = H5Sget_simple_extent_ndims(space);
  hsize_t dims[ndims];
  H5Sget_simple_extent_dims(space, dims, NULL);
  hsize_t size = unitSize;
  hid_t type = H5Dget_type(dset);
  printf("Size of type %ld\n", H5Tget_size(type));
  //size = H5Tget_size(type);
  for(int i = 0; i < ndims; i++)
  {
    size *= dims[i];
  }
  printf("%d\n", ndims);
  if(ndims == 3)
  {
    printf("%d %d %d\n", dims[0], dims[1], dims[2]);
  }else if(ndims == 2)
  {
    printf("%d %d\n", dims[0], dims[1]);
  }

  *target = malloc(size);
  printf("Hello %s\n", set_name);
  status = H5Dread (dset, h5DataType, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                *target);
  H5Dclose(dset);
  H5Sclose(space);
}
void importFromH5()
{
  //import ray info (areas, trajectories, and vectors) from MATLAB results
  free(boxes);
  free(areas);
  char* path = "matlabcbet.h5";
  hid_t file;
  herr_t status;
  file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
  double* tempArea;
  importFromH5(file, "trajectories", (void**)&boxes, sizeof(int), H5T_NATIVE_INT);
  importFromH5(file, "areaRatios", (void**)&tempArea, sizeof(double), H5T_NATIVE_DOUBLE);
  importFromH5(file, "rayVectors", (void**)&ray_k, sizeof(double), H5T_NATIVE_DOUBLE);
  importFromH5(file, "ds", (void**)&dkmag, sizeof(double), H5T_NATIVE_DOUBLE);
  for(int i = 0; i < CROSS*3; i++)
  {
    //printf("%f\n", ray_k[i]);
  }
  dkx = new double[CROSS]{0.0};
  dkz = new double[CROSS]{0.0};
  areas = new double[CROSS]{0.0};
  for(int i = 0; i < nbeams; i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      for(int k = 0; k < ncrossings; k++)
      {
        double* value = vec4DP(ray_k, i, j, k, 0, nrays, 12, 3);
        double area = vec3D(tempArea, i, j, k, nrays, 12);
        double kx = value[0], kz = value[1];

        int rayindex = j;
        if(i == 1)
        {
          //rayindex == 32-j-1;
        }
        vec3DW(dkx, i, rayindex, k, nrays, ncrossings, kx);
        vec3DW(dkz, i, rayindex, k, nrays, ncrossings, kz);
        vec3DW(areas, i, rayindex, k, nrays, ncrossings, area);
      }
    }
  }
  getchar();
  free(tempArea);
  H5Fclose(file);
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
  
  initialize();
  auto start1 = chrono::high_resolution_clock::now();
  if(argc > 4)
  {
    pyPlot = argv[1][0] - 48;
    std::string str(argv[2]);
    threads = std::stoi(str);
    cudaCalc = std::stoi(string(argv[3]));
    rays_per_zone = std::stoi(string(argv[4]));
    uray_mult = intensity*(courant_mult)*pow(double(rays_per_zone),-1.0); //multiplier which determines intensity deposited in a given zone
    nrays= int(rays_per_zone*(beam_max_z-beam_min_z)/dz)+0;//number of rays per beam
    /*
    nx = std::stoi(string(argv[4]));
    nz = nx;
    dz = (zmax-zmin)/(nz-1);//dimensions of a single cell
    dx = (xmax-xmin)/(nx-1);
    nrays= int(rays_per_zone*(beam_max_z-beam_min_z)/dz)+0;//number of rays per beam
    dt=courant_mult*fmin(dx,dz)/c;//time stepping
    nt=int(pow(courant_mult,-1.0)*fmax(nx,nz)*2.0)+1;//number of time steps to track for a given ray
    numstored = nx*6;//number of rays stored per grid zone
    ncrossings = nx * 3;//max number of ray crossings that can be stored per ray
    freq = c/lambda;		// frequency of light, in Hz
    */
  }
  else if(argc > 2)
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
  if(!cudaCalc)
  {
    //printf("Threads: %d\n", threads);
  }
  optimize = 0;
  int b = 0;
  if(b)
  {
    benchmark();
    return 0;
  }
  auto stop1 = chrono::high_resolution_clock::now();
  auto start2 = chrono::high_resolution_clock::now();
  importFromH5();
  fillTraceArrays();
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
  printf("check %d\n", 3);
  fflush(stdout);
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
  printf("check %d\n", 4);
  fflush(stdout);
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
  printf("check %d\n", 5);
  fflush(stdout);
  return 0;
}
