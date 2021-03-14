#ifndef CBET
#define CBET
  #include "implSim.hpp"
  #include "customMath.hpp"
  #include "cuda_help.hpp"
  #include <cmath>
  #include <iomanip>
  #include <cstring>
  extern void initArrays();
  extern void launchCBETKernel();

  struct CBETVars
  {
    double estat_cu;
    double me_cu;
    double c_cu;
    double kb_cu;
    double Te_cu;
    double Ti_cu;
    int nbeams_cu;
    int nx_cu;
    int nz_cu;
    double dx_cu;
    double dz_cu;
    int nrays_cu;
    int ncrossings_cu;
    int numstored_cu;
    double ne_cu;
    double ncrit_cu;
    double pi_cu;
    double iaw_cu;
    double cs_cu;
    double Z_cu;
    double omega_cu;
  };
  struct CBETArrs
  {
    double* eden_cu;
    double* dkx_cu;
    double* dkz_cu;
    double* dkmag_cu;
    double* uflow_cu;
    double* i_b_cu;
    double* i_b_new_cu;
    double* W_new_cu;
    double* W_cu;
    double* x_cu;
    double* z_cu;
    int* numrays_cu;
    int* ints_cu;
    int* boxes_cu;
    int* present_cu;

  };
  extern CBETArrs* new_cbetArrs();
  extern CBETVars* new_cbetVars();
#endif
