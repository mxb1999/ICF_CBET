#ifndef TRACE
#define TRACE
  #include "implSim.hpp"
  #include "customMath.hpp"
  extern void launchRays();
  extern void track(int raynum, int xinit, int zinit, int kxinit, int kzinit, double urayinit, int beam);
  extern omp_lock_t writelock;
#endif
