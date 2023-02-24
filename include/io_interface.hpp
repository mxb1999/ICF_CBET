

#ifndef IODEF
#define IODEF
  #include <hdf5.h>
  #include "H5Cpp.h"
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
