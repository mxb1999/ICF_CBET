#ifndef CBET
#define CBET
  #include "dataFields.hpp"
  #include "customMath.hpp"
  #include "cuda_help.hpp"
  #include <cmath>
  #include <iomanip>
  #include <cstring>
  extern void initArrays();
  extern void launchCBETKernel();
  extern void freeCBETArrs();
  extern void cbetOptimize();
//Structs for passing constants to kernel launches
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
    double* dk_cu;
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
  inline double beamProfile(double w, double n, double r)
  {
    return exp(-2*pow(abs(r/w),n));
  }
  //
  extern CBETArrs* new_cbetArrs();
  extern CBETVars* new_cbetVars();
  extern double* direct_map(double* cbetAmp);
#endif
