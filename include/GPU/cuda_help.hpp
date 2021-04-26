
  
//define device variables
//define other functions necessary for GPU operation
#ifndef CUDAHELP
#define CUDAHELP
  #include <vector>
  #include <string>
  #include <cuda_runtime.h>
  //extern __device__ double c_cu, lambda_cu, estat_cu,mach_cu,Z_cu,mi_cu,mi_kg_cu,Te_cu,Te_eV_cu,Ti_cu,Ti_eV_cu,iaw_cu,ncrit_cu,freq_cu,omega_cu;
  //spatial information
  //extern __device__ int nx_cu, nz_cu;
  //extern __device__ double xmin_cu, xmax_cu, zmin_cu, zmax_cu, dx_cu, dz_cu;

  //extern __device__ double maxIncr_cu, converge_cu;

  //extern __device__ int nbeams_cu, nrays_cu, nt_cu, numstored_cu, rays_per_zone_cu, ncrossings_cu;

  //extern __device__ double intensity_cu, offset_cu, uray_mult_cu, beam_min_z_cu, beam_max_z_cu, dt_cu, courant_mult_cu;
  extern std::vector<void*> hostP;
  extern std::vector<void*> devP;
  extern std::vector<std::string> varNames;
  extern std::vector<size_t> sizes;
  


  //define a struct to hold TrackRay Constants->passing once by reference should suffice
  struct TrackConst
  {
    //spatial information
    double dx_cu;
    double dz_cu;
    int nx_cu;
    int nz_cu;
    double xmax_cu;
    double xmin_cu;
    double zmax_cu;
    double zmin_cu;
    //ray params
    int nrays_cu;
    int ncrossings_cu;
    int nbeams_cu;
    int numstored_cu;
    int nt_cu;
    double dt_cu;
    double omega_cu;
    double c_cu;
    double ncrit_cu;
  };
  struct LinkCross
  {
    void* next;
    double vals[9];
    unsigned int cross;
    unsigned int location;//save space, index by modular arithmetic
  };
  extern __device__ void add_LinkCross(LinkCross* entry, LinkCross* newVal);
  extern  void add_LinkCrossHost(LinkCross* entry, LinkCross* newVal);

  extern __device__ LinkCross* new_LinkCross(double* values, unsigned int  ix, unsigned int  iz, unsigned int cross, int nz);
  extern LinkCross* new_LinkCrossHost(double* values,  unsigned int  ix, unsigned int  iz,unsigned int cross);
  //define a struct to hold TrackRay Constants->passing once by reference should suffice
  struct TrackArrs
  {
    //spatial information
    double* dedendx_cu;
    double* dedendz_cu;
    //char* rayZones_cu;
    double* x_cu;
    double* z_cu;
    double* crossesx_cu;
    double* crossesz_cu;
    LinkCross** edep_cu;
    double* wpe_cu;
    int* boxes_cu;
    int* ints_cu;
  };
  extern TrackArrs* deviceTrackArrs(int device);
  extern TrackConst* deviceTrackConst(int device);
  
  //inlined flattened array accessor/mutator methods
  template <typename T>
  __forceinline__ __device__ 
  T vec4D_cu(T* arr, int a, int b, int c, int d, int d2, int d3, int d4)
  {
    return arr[(((a)*d2+b)*d3+c)*d4+d];//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }
  template <typename T>
  __forceinline__ __device__
  void vec5D_cu(T* arr, int a, int b, int c, int d,int e, int d2, int d3, int d4,int d5, T val)
  {
    return arr[((((a)*(d2)+b)*(d3)+c)*(d4)+d)*(d5)+e];//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
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
  void vec5DP_cu(T* arr, int a, int b, int c, int d,int e, int d2, int d3, int d4,int d5, T val)
  {
    return arr+((((a)*(d2)+b)*(d3)+c)*(d4)+d)*(d5)+e;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
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
  void vec5DW_cu(T* arr, int a, int b, int c, int d,int e, int d2, int d3, int d4,int d5, T val)
  {
    arr[((((a)*(d2)+b)*(d3)+c)*(d4)+d)*(d5)+e] = val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
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
    arr[(a)*(d2)+(b)] = val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }
  template <typename T>
  __forceinline__ __device__
  void vec5DI_cu(T* arr, int a, int b, int c, int d,int e, int d2, int d3, int d4,int d5, T val)
  {
    arr[((((a)*(d2)+b)*(d3)+c)*(d4)+d)*(d5)+e] += val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }
  template <typename T>
  __forceinline__ __device__
  void vec4DI_cu(T* arr, int a, int b, int c, int d, int d2, int d3, int d4, T val)
  {
    arr[(((a)*(d2)+(b))*(d3)+(c))*(d4)+(d)] += val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }
  template <typename T>
  __forceinline__ __device__
  void vec3DI_cu(T* arr, int a, int b, int c, int d2, int d3, T val)
  {
    arr[((a)*(d2)+(b))*(d3)+(c)] += val;
  }
  template <typename T>
  __forceinline__ __device__
  void vec2DI_cu(T* arr, int a, int b, int d2, T val)
  {
    arr[(a)*(d2)+(b)] += val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }
  template <typename T>
  __forceinline__ __device__
  void vec5DM_cu(T* arr, int a, int b, int c, int d,int e, int d2, int d3, int d4,int d5, T val)
  {
    arr[((((a)*(d2)+b)*(d3)+c)*(d4)+d)*(d5)+e] *= val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }
  template <typename T>
  __forceinline__ __device__
  void vec4DM_cu(T* arr, int a, int b, int c, int d, int d2, int d3, int d4, T val)
  {
    arr[(((a)*(d2)+(b))*(d3)+(c))*(d4)+(d)] *= val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }
  template <typename T>
  __forceinline__ __device__
  void vec3DM_cu(T* arr, int a, int b, int c, int d2, int d3, T val)
  {
    arr[((a)*(d2)+(b))*(d3)+(c)] *= val;
  }
  template <typename T>
  __forceinline__ __device__
  void vec2DM_cu(T* arr, int a, int b, int d2, T val)
  {
    arr[(a)*(d2)+(b)] *= val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
  }

  
#endif