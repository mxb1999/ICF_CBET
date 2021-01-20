#include "implSim.h"
#include "io_interface.hpp"
//define device variables
//define other functions necessary for GPU operation
#ifndef CUDAHELP
#define CUDAHELP
    extern CGPUMemCpy
    extern CGPUMalloc()
#endif