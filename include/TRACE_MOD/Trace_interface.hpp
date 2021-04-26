#ifndef TRACE
#define TRACE
  #include <iostream>
  #include "dataFields.hpp"
  #include "customMath.hpp"
  #include "parallelConfig.hpp"
  #include "cuda_help.hpp"
  #define Q1 4
  #define Q2 3
  #define Q3 1
  #define Q4 2
  struct rayinit//composite struct, compactly store necessary initial information
  {
      double xinit;
      double zinit;
      double kxinit;
      double kzinit;
      int beam;
      double urayinit;
      double wpeinit;
  };
  class RayBundle
  {
    static int numRays;
    private:
      int* rays;
      int intesin;
  };

  extern void launchRays();
  extern void track(int raynum, int xinit, int zinit, int kxinit, int kzinit, double urayinit, int beam);
  extern void LaunchCUDARays(rayinit* rays);
  extern void fillTraceArrays();
  extern void launch_ray_XZ(rayinit initSettings, int raynum);
  extern void freeTraceArrs();
  extern void optimizePostTraceArrs();
  extern void optimizePreTraceArrs();
#endif
