//define device variables
//define other functions necessary for GPU operation
#ifndef CUDAHELP
#define CUDAHELP
    #include "implSim.hpp"
    
    extern void CGPUMemCpy(void* dest, void* source, size_t size, int direction);
    extern void CGPUAlloc(void* p, size_t size);
    extern void CGPUDataTrans(void** hostP, void** deviceP, size_t* sizes, int n, int dir);//function to transfer large quantities of variables from host to GPU, makes bulk transfers easier


    extern void* rays;
    extern __device__ void* rays_cu;//rayinit struct, to be declared later
    extern __device__ double* x_cu;
    extern __device__ double* z_cu;
    extern __device__ double* eden_cu;
    extern __device__ double* wpe_cu;
    extern __device__ double* dedendx_cu;
    extern __device__ double* dedendz_cu;
    extern __device__ bool* marked_cu;//nbeams nrays nx nz
    extern __device__ int* boxes_cu;
    extern __device__ double c_cu, lambda_cu, __device__ estat_cu,mach_cu,Z_cu,mi_cu,mi_kg_cu,Te_cu,Te_eV_cu,Ti_cu,Ti_eV_cu,iaw_cu,ncrit_cu,freq_cu,omega_cu;
    //spatial information
    extern __device__ int nx_cu, nz_cu;
    extern __device__ double xmin_cu, xmax_cu, zmin_cu, zmax_cu, dx_cu, dz_cu;

    extern __device__ double maxIncr_cu, converge_cu;

    extern __device__ int nbeams_cu, nrays_cu, nt_cu, numstored_cu, rays_per_zone_cu, ncrossings_cu;

    extern __device__ double intensity_cu, offset_cu, uray_mult_cu, beam_min_z_cu, beam_max_z_cu, dt_cu, courant_mult_cu;
    extern __device__ double* edep_cu, *crossx_cu, *crossz_cu;
    

  template <typename T>
  __forceinline__ __device__ 
  T vec4D_cu(T* arr, int a, int b, int c, int d, int d2, int d3, int d4)
  {
    return arr[(((a)*d2+b)*d3+c)*d4+d];//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }

  template <typename T>
  __forceinline__ __device__ 
  T vec3D_cu(T* arr, int a, int b, int c, int d2, int d3)
  {
    return arr[((a)*d2+b)*d3+c];
  }
  template <typename T>
  __forceinline__ __device__ 
  T vec2D_cu(T* arr, int a, int b, int d2)
  {
    return arr[(a)*d2+b];//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }
  template <typename T>
   __forceinline__ __device__
  T* vec4DP_cu(T* arr, int a, int b, int c, int d, int d2, int d3, int d4)
  {
    return arr + (((a)*d2+b)*d3+c)*d4+d;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }

  template <typename T>
  __forceinline__ __device__
  T* vec3DP_cu(T* arr, int a, int b, int c, int d2, int d3)
  {
    return arr + ((a)*d2+b)*d3+c;
  }
  template <typename T>
  __forceinline__ __device__
  T* vec2DP_cu(T* arr, int a, int b, int d2)
  {
    return arr + (a)*d2+b;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }

  template <typename T>
  __forceinline__ __device__
  void vec4DW_cu(T* arr, int a, int b, int c, int d, int d2, int d3, int d4, T val)
  {
    arr[(((a)*d2+b)*d3+c)*d4+d] = val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }

  template <typename T>
  __forceinline__ __device__
  void vec3DW_cu(T* arr, int a, int b, int c, int d2, int d3, T val)
  {
    arr[((a)*d2+b)*d3+c] = val;
  }
  template <typename T>
  __forceinline__ __device__
  void vec2DW_cu(T* arr, int a, int b, int d2, T val)
  {
    arr[(a)*d2+b] = val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }

  template <typename T>
  __forceinline__ __device__
  void vec4DI_cu(T* arr, int a, int b, int c, int d, int d2, int d3, int d4, T val)
  {
    arr[(((a)*d2+b)*d3+c)*d4+d] += val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }
  template <typename T>
  __forceinline__ __device__
  void vec3DI_cu(T* arr, int a, int b, int c, int d2, int d3, T val)
  {
    arr[((a)*d2+b)*d3+c] += val;
  }
  template <typename T>
  __forceinline__ __device__
  void vec2DI_cu(T* arr, int a, int b, int d2, T val)
  {
    arr[(a)*d2+b] += val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }

  
#endif