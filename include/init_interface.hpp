#ifndef INIT
#define INIT
  #include "implSim.hpp"
  #include <vector>
  #include <queue>
  #include <map>
  int printUpdates,printTimings,printCBETDiagnostics,printRayTrackDiagnostics,printHDF5Diagnostics,printSpecificTimings, pyPlot,iterate,calcCBET;
  double lambda,estat,mach,Z,mi,mi_kg,Te,Te_eV,Ti,Ti_eV,iaw,ncrit,freq,omega;
  int switchvar;
  //spatial information
  int nx, nz;
  double xmin, xmax, zmin, zmax, dx, dz;

  double maxIncr, converge;
  int threads,maxIter;

  int** raypath;
  int nbeams, nrays, nt, numstored, rays_per_zone, ncrossings;

  double intensity, offset, uray_mult, beam_min_z, beam_max_z, dt, courant_mult;
  //debugging
  double** gain2arr;
  double** gain1arr;
  double** mag;
  short** rayAdjList;

//Used to instantiate the variables declared in implSim.h without taking up half of the Initialize file
double maxDev;

//Launch Ray Values


double cs;
int ray1num;
double maxInc;
int iter;
int count = 0;
double injected = 0;
int gridcount = 0;
int** counter;
//Pointers for necessary arrays

double** intersections; //nx nz
int* marked; //nx nz numstored nbeams
double** dedendx; //nx nz
double** dedendz; //nx nz
double* x; //nx
double* z; //nz
double** eden; //nx nz

double*** edep; //nx+2 nz+2 nbeams
int*** present; //nx nz nbeams
double** machnum; //nx nz
int* boxes; //nbeams nrays nx*3 2
double** u_flow; //nx nz
double*** dkx; //nbeams nrays 2
double*** dkz; //nbeams nrays 2
double*** dkmag; //nbeams nrays 2

//Launch_Ray_XZ specific arrays (all have a length of nt)

//CBET specific arrays
double*** W;//nx nz
double*** W_init;//nx nz
double*** W_new;//nx nz
double*** i_b;//nx nz
double*** i_b_prev; //nbeams nx nz
double*** i_b_new;//nx nz
double** wpe; //nx nz
double*** crossesz; //nbeams nrays ncrossings
double*** crossesx; //nbeams nrays ncrossings
int*** ints; //nbeams nrays ncrossings
int iteration = 0;
//arrays used only for plotting
double** i_bplot;//nx nz
double* convergeplot; //maxWrites
double** orderplot1; //nx nz
double** orderplot2; //nx nz
double** i_b1Error;
double** i_b2Error;
double* mult;//nbeams nbeams-1 nrays ncrossings
double** i_b_newplot;//nx nz
double** edenplot; //the array is eden/ncrit,  nx nz
double** edepplot; //nx nz
int** raytrace;//nx nz
double** ib_orig;//nx nz
int** anyInt;//nx nz
double** perturbation;//nx nz
  //derivative values
  //double mi_kg,mi,uray_mult,dt,dx,dz,freq,omega,ncrit;
  //int nrays,nt,numstored,ncrossings;

#endif