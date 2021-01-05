#include <queue>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <cmath>
#include "customMath.hpp"
#include "simConst.hpp"
//#include "LinkedList.h"

#ifndef IMPLSIM_H_
#define IMPLSIM_H_
    //Functions
  void initialize();
  void launchRays();
  void launch_ray_XZ(double x_init, double z_init, double kx_init, double kz_init,double urayinit, int raynum);
  void cbet();
  void updateH5();
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
  extern int beam;//stores which beam is currently being tracked
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
  extern Ray** beamIndex;//nbeams nrays
  extern Ray*** spatialIndex; //nbeams nx nz
  extern queue<int>* marked; //nx nz nbeams Stored as a 1D array, saved SIGNIFICANT amount of time in initialization
  extern double** dedendx; //nx nz
  extern double** dedendz; //nx nz
  extern double* x; //nx
  extern double* z; //nz
  extern double** eden; //nx nz
  //Marked Stores which rays from each beam have passed through each zone
  //Boxes Stores whi
  extern double*** edep; //nx+2 nz+2 nbeams
  extern int*** present; //nx nz nbeams
  extern double** machnum; //nx nz
  extern int**** boxes; //nbeams nrays ncrossings 2
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
  extern double** i_bplot;//nx nz
  extern double** orderplot1; //nx nz
  extern double** orderplot2; //nx nz
  extern double** i_b1Error; //nx nz
  extern double* convergeplot;//nx nz
  extern double** i_b2Error;//nx nz
  extern double** i_b_newplot;//nx nz
  extern double** edenplot; //the array is eden/ncrit,  nx nz
  extern double** edepplot; //nx nz
  extern double** raytrace;//nx nz
  extern double** ib_orig;//nx nz
  extern int** anyInt;//nx nz
  extern double** perturbation;//nx nz
#endif
