

#ifndef IODEF
#define IODEF

  #include <hdf5/serial/hdf5.h>
  #include <hdf5/serial/H5Cpp.h>
  #include <vector>
  #include <iostream>
  #include <queue>
  #include <map>
  #include "dataFields.hpp"
  #include "parallelConfig.hpp"
  extern void updateH5();//Export simulation data to HDF5 file
  extern void initialize();
using namespace std;
#endif
