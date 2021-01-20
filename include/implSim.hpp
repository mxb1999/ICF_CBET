#include <queue>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <cmath>
#include <cuda_runtime.h>

#include "customMath.hpp"
//#include "LinkedList.h"

#ifndef IMPLSIM_H_
#define IMPLSIM_H_
  template <typename T>
  inline T* vec4D(T* arr, int a, int b, int c, int d, int d2, int d3, int d4)
  {
    return &arr[(((a)*d2+b)*d3+c)*d4+d];//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }
    //Functions
  extern void initialize();
  extern void launchRays();
  extern void launch_ray_XZ(double x_init, double z_init, double kx_init, double kz_init,double urayinit, int raynum, int beam);
  extern void cbet();
  extern void updateH5();
  typedef struct Intersection Intersection;
  typedef struct Ray Ray;
  struct Crossing
  {
      int x;
      int z;
  };
  struct Ray
  {
    int index;
    //LinkedList* crossings;
    double intensity;
  };

  //Values needed throughout simulation
  extern double maxDev;//stores maximum change in a given iteration, determines convergence
  //extern int beam;//stores which beam is currently being tracked
  extern int iter;//stores iteration number, not currently used
  extern int count;//Used to track a specific grid square along with counter, tracking change over time
  extern int** counter;
  //Launch Ray Values
  extern double cs;
  extern double injected;
  extern int gridcount;
  extern double** intersections;
  extern int ray1num;
  extern double maxInc;
  //Pointers for necessary arrays
  //extern double** intersections; //nx nz
  //can store both marked and boxes in an adjacency list of Intesection objects?
  /*capabilities needed:
    track future path of ray across other zones -> boxes
    track which rays have entered a given zone -> marked
  */
  /*
  Use crossings adjacency list to get crossings, can use/update sampleIntensity
  */
  extern short** rayAdjList;//ith level has (nbeams-1-i)*nrays*nrays ints Try short for now, doubtful that over 32,000 rays will be used atm for each beam
  extern Ray** beamIndex;//nbeams nrays
  extern Ray*** spatialIndex; //nbeams nx nz
  extern int* marked; //nx nz nbeams Stored as a 1D array, saved SIGNIFICANT amount of time in initialization
  extern double** dedendx; //nx nz
  extern double** dedendz; //nx nz
  extern int** raypath; //nx nz store single ray path
  extern double* x; //nx
  extern double* z; //nz
  extern double** eden; //nx nz
  //Marked Stores which rays from each beam have passed through each zone
  //Boxes Stores whi
  extern double*** edep; //nx+2 nz+2 nbeams
  extern int*** present; //nx nz nbeams
  extern double** machnum; //nx nz
  extern int* boxes; //nbeams nrays ncrossings 2
  extern bool**** boxTrack;
  extern double**** W_storage; //nbeams nx nz nrays
  extern double** u_flow; //nx nz
  extern double*** dkx; //nbeams nrays 2
  extern double*** dkz; //nbeams nrays 2
  extern double*** dkmag; //nbeams nrays 2

  //Launch_Ray_XZ specific arrays (all have a length of nt)
  //CBET specific arrays
  extern double*** W;//nx nz
  extern double*** W_init;//nx nz
  extern double*** W_new;//nx nz
  extern double*** i_b;//nx nz
  extern double*** i_b_prev;//nbeams nx nz
  extern double*** i_b_new;//nx nz
  extern double** wpe; //nx nz
  extern double*** crossesz; //nbeams nrays ncrossings
  extern double*** crossesx; //nbeams nrays ncrossings
  extern int*** ints; //nbeams nrays ncrossings
  extern int iteration;
  //arrays used only for plotting
  extern double** gain2arr;
  extern double** gain1arr;
  extern double** mag;
  extern double** i_bplot;//nx nz
  extern double** orderplot1; //nx nz
  extern double** orderplot2; //nx nz
  extern double** i_b1Error; //nx nz
  extern double* convergeplot;//nx nz
  extern double** i_b2Error;//nx nz
  extern double** i_b_newplot;//nx nz
  extern double** edenplot; //the array is eden/ncrit,  nx nz
  extern double** edepplot; //nx nz
  extern int** raytrace;//nx nz
  extern double** ib_orig;//nx nz
  extern int** anyInt;//nx nz
  extern double** perturbation;//nx nz
  extern double* mult;//nbeams nrays ncrossings
  //SIMULATION CONSTANTS
  extern int printUpdates;
  extern int printTimings;
  extern int printCBETDiagnostics;
  extern int printRayTrackDiagnostics;
  extern int printHDF5Diagnostics;
  extern int printSpecificTimings;
  extern int pyPlot;
  extern int iterate;
  extern int calcCBET;

  extern int switchvar;
  extern double lambda;
  extern double estat;
  extern double mach;
  extern double Z;
  extern double mi;
  extern double mi_kg;
  extern double Te;
  extern double Te_eV;
  extern double Ti;
  extern double Ti_eV;
  extern double iaw;
  extern double ncrit;
  extern double freq;
  extern double omega;

  //spatial information
  extern int nx;
  extern int nz;
  extern double xmin;
  extern double xmax;
  extern double zmin;
  extern double zmax;
  extern double dx;
  extern double dz;

  extern double maxIncr;
  extern double converge;
  extern int maxIter;
  extern int threads;


  extern int nbeams;
  extern int nrays;
  extern int nt;
  extern int numstored;
  extern int rays_per_zone;
  extern int ncrossings;

  extern double courant_mult;
  extern double intensity;
  extern double offset;
  extern double uray_mult;
  extern double beam_min_z;
  extern double beam_max_z;
  extern double dt;

  //Fundamental Constants
  const double sigma = 1.7e-4;
  const double e0 =8.85418782e-12;
  const double me =9.10938356e-31;
  const double pi =3.14159265359;
  const double kb = 1.3806485279e-16;   //Boltzmann constant in erg/K
  const double kb2 = 1.3806485279e-23;   //Boltzmann constant in J/K
  const double ec = 1.60217662e-19;
  const double c = 29979245800.0;              // Speed of light in cm/s

  /*const double Ti_eV = 1.0e3;
  const double iaw = 0.2;                      // ion-acoustic wave energy-damping rate (nu_ia/omega_s)!!

  const double freq = c/lambda;		// frequency of light, in Hz
  const double omega = 2*pi*freq;	// frequency of light, in rad/s
  const double ncrit = 1e-6*(pow(omega,2.0)*me*e0/pow(ec,2.0));

  //Grid Constants for X by Z grid
  const int nx=201;
  const float xmin = -5.0e-4;
  const float xmax=5.0e-4;
  const float dx = (xmax-xmin)/(nx-1);
  const int nz=201;
  const float zmin = -5.0e-4;
  const float zmax=5.0e-4;
  const float dz = (zmax-zmin)/(nz-1);

  //Constants pertaining to iteration and parallelization
  const double maxIncrement = 0.2;
  const int maxIterations = 100*iterate;
  //Number of parallel threads
  const int threads = 12;
  //Fractional convergence cutoff (when the field stops changing by (converge*100) %)
  const double converge = 1e-3;

  //Beam/Ray Tracking parameters
  const int nbeams = 2; //number of interacting beams
  const int doubleer_zone = 5; //Rays launched per grid zone
  const double intensity = 1.0e17; // W/cm^2
  const float courant_mult = 0.2; // 0.37 // 0.25 // 0.36 // 0.22; Multiplier used to determine time stepping
  const double offset = 0;//Determines offset of beam along axis
  const double beam_max_z = 3.0e-4; const double beam_min_z = -3.0e-4;//determines the width of the beam
  const int nrays= int(rays_per_zone*(beam_max_z-beam_min_z)/dz)+0;//number of rays per beam
  const double dt=courant_mult*fmin(dx,dz)/c;//time stepping
  const int nt=int(pow(courant_mult,-1.0)*fmax(nx,nz)*2.0)+1;//number of time steps to track for a given ray
  const int numstored = nx*6;//number of rays stored per grid zone
  const int ncrossings = nx * 3;//max number of ray crossings that can be stored per ray

#endif
//Used to determine print statements, 1 = on, 0 = off
const int printUpdates = 1;
const int printTimings = 1;
const int printCBETDiagnostics = 1;
const int printRayTrackDiagnostics = 0;
const int printHDF5Diagnostics = 0;
const int printSpecificTimings = 1;
const int pyPlot = 1; //Determines whether the python script will automatically
const int iterate = 1;//determine whether the simulation will iterate
const int calcCBET = 1;

const double lambda = 3.51e-5;	// wavelength of light, in cm. This is frequency-tripled "3w" or "blue" (UV) light
const double estat=4.80320427e-10; 	       // electron charge in statC
const double mach = -1.0*sqrt(2);                 // Mach number for max resonance
const double Z = 3.1;                        // ionization state
const double Te = 2.0e3*11604.5052;          // Temperature of electron in K
const double Te_eV = 2.0e3;
const double Ti = 1.0e3*11604.5052;          // Temperature of ion in K
const double mi_kg = 10230.0*me;	   // Mass of ion in kg
const double mi = 10230*(1.0e3*me);          // Mass of ion in g
const double uray_mult = intensity*(courant_mult)*pow(double(rays_per_zone),-1.0); //multiplier which determines intensity deposited in a given zone
const int nrays= int(rays_per_zone*(beam_max_z-beam_min_z)/dz)+0;//number of rays per beam
const double dt=courant_mult*fmin(dx,dz)/c;//time stepping
const int nt=int(pow(courant_mult,-1.0)*fmax(nx,nz)*2.0)+1;//number of time steps to track for a given ray
const int numstored = nx*6;//number of rays stored per grid zone
const int ncrossings = nx * 3;//max number of ray crossings that can be stored per ray
const float dz = (zmax-zmin)/(nz-1);
const float dx = (xmax-xmin)/(nx-1);
const double freq = c/lambda;		// frequency of light, in Hz
const double omega = 2*pi*freq;	// frequency of light, in rad/s
const double ncrit = 1e-6*(pow(omega,2.0)*me*e0/pow(ec,2.0));
*/
#endif
