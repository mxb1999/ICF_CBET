#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <omp.h>
#include <chrono>
#include <cmath>
#include "customMath.h"
#include "simConst.h"
#include <hdf5/serial/H5Cpp.h>
#ifndef IMPLSIM_H_
#define IMPLSIM_H_
    //Functions
  void initialize();
  void launchRays();
  void launch_ray_XZ(double x_init, double z_init, double kx_init, double kz_init,double urayinit, int raynum);
  void cbet();
  void updateH5();
  //Values needed throughout simulation
  extern double maxDev;
  extern int beam;
  extern int iter;
  extern int count;
  extern int** counter;
  //Launch Ray Values


  extern double cs;
  extern double injected;
  extern int gridcount;
  extern int ray1num;
  extern double maxInc;
  //Pointers for necessary arrays
  extern double** intersections; //nx nz
  //marked stores the trajectory of a given ray
  extern int** marked; //nx nz numstored nbeams
  extern double** dedendx; //nx nz
  extern double** dedendz; //nx nz
  extern double* x; //nx
  extern double* z; //nz
  extern double** eden; //nx nz

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
  extern double** i_b1Error;
  extern double* convergeplot;
  extern double** i_b2Error;
  extern double** i_b_newplot;//nx nz
  extern double** edenplot; //the array is eden/ncrit,  nx nz
  extern double** edepplot; //nx nz
#endif
