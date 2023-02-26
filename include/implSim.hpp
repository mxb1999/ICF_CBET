






#ifndef IMPLSIM_H_
#define IMPLSIM_H_
  #include "parallelConfig.hpp"
  #include "customMath.hpp"
  #include "dataFields.hpp"
  #include "Trace_interface.hpp"
  #include "io_interface.hpp"
  #include "CBET_Interface.hpp"
  #define VEC2D(arr, i, j, s2) (arr)[(i)*(s2) + (j)]
  #define VEC3D(arr, i, j, k, s2, s3) (arr)[((i)*(s2) + (j))*(s3) + (k)]
  #define VEC4D(arr, i, j, k, l, s2, s3, s4) (arr)[(((i)*(s2) + (j))*(s3) + (k))*(s4) + (l)]
  #define GRID nx*nz*ny
  #define RAYS nbeams*nrays
  #define CROSS nbeams*nrays*ncrossings

  //Functions
  extern void initialize();
  extern void launchRays();
void cbet(double* wMult, double* i_b, double* eden, double* areas, double* machnum, int* boxes, double* kvec, double medianDS, int* marked, double tolerance);
  extern void updateH5();


#endif
