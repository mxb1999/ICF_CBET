






#ifndef IMPLSIM_H_
#define IMPLSIM_H_ 
  #include "parallelConfig.hpp"
  #include "customMath.hpp"
  #include "dataFields.hpp"
  #include "Trace_interface.hpp"
  #include "io_interface.hpp"
  #include "CBET_Interface.hpp"

  #define GRID nx*nz
  #define RAYS nbeams*nrays
  #define CROSS nbeams*nrays*ncrossings
  
  //Functions
  extern void initialize();
  extern void launchRays();
  extern void cbet();
  extern void updateH5();

  
#endif
