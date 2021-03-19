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

  //better life through macros, used frequently in device CBET functions
 /* #define CBET_IMPORT(CBETVars constants, CBETArrs arrays) ({\
  //imported constants
          int ncrossings_cu = constants.ncrossings_cu;\
          int nx_cu = constants.nx_cu;\
          int nz_cu = constants.nz_cu;\
          int numstored_cu = constants.numstored_cu;\
          double dx_cu = constants.dx_cu;\
          double dz_cu = constants.dz_cu;\
          double ncrit_cu = constants.ncrit_cu;\
          double c_cu = constants.c_cu;\
          double pi_cu = constants.pi_cu;\
          double iaw_cu = constants.iaw_cu;\
          double cs_cu = constants.cs_cu;\
          double estat_cu = constants.estat_cu;\
          double Ti_cu = constants.Ti_cu;\
          double Te_cu = constants.Te_cu;\
          double Z_cu = constants.Z_cu;\
          double omega_cu = constants.omega_cu;\
          double kb_cu = constants.kb_cu;\
          double me_cu = constants.me_cu;\

          //imported arrays
          double* i_b_cu = arrays.i_b_cu;\
          double* i_b_new_cu = arrays.i_b_new_cu;\

          double* W_cu = arrays.W_cu;\

          double* x_cu = arrays.x_cu;\
          
          double* z_cu = arrays.z_cu;\
          double* W_new_cu = arrays.W_new_cu;\
          double* dkx_cu = arrays.dkx_cu;\
          double* dkz_cu = arrays.dkz_cu;\
          double* dkmag_cu = arrays.dkmag_cu;\
          double* uflow_cu = arrays.uflow_cu;\

          int* ints_cu = arrays.ints_cu;\
          double* eden_cu = arrays.eden_cu;\
          int* boxes_cu = arrays.boxes_cu;\
          int* numrays_cu = arrays.numrays_cu;\
          int* present_cu = arrays.present_cu;\
      })*/
  extern CBETArrs* new_cbetArrs();
  extern CBETVars* new_cbetVars();
#endif
