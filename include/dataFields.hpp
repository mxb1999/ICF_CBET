  #ifndef DATAFIELDS
  #define DATAFIELDS
    #include "DEV_MOD/parallelConfig.hpp"
    #include <queue>
    #include <iostream>
    #include <vector>
    #include <string>
    #include <fstream>
    #include <omp.h>
    #include <chrono>
    #include <cmath>
    #define GRID nx*nz
    #define RAYS nbeams*nrays
    #define CROSS nbeams*nrays*ncrossings
    extern GConfig* deviceConfiguration;
    template <typename T>
    inline T vec4D(T* arr, int a, int b, int c, int d, int d2, int d3, int d4)
    {
        return arr[(((a)*d2+b)*d3+c)*d4+d];//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }

    template <typename T>
    inline T vec3D(T* arr, int a, int b, int c, int d2, int d3)
    {
        return arr[((a)*(d2)+b)*(d3)+c];
    }
    template <typename T>
    inline T vec2D(T* arr, int a, int b, int d2)
    {
        return arr[(a)*(d2)+b];//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }
    template <typename T>
    inline T* vec4DP(T* arr, int a, int b, int c, int d, int d2, int d3, int d4)
    {
        return arr + (((a)*(d2)+b)*(d3)+c)*(d4)+d;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }

    template <typename T>
    inline T* vec3DP(T* arr, int a, int b, int c, int d2, int d3)
    {
        return arr + ((a)*(d2)+b)*(d3)+c;
    }
    template <typename T>
    inline T* vec2DP(T* arr, int a, int b, int d2)
    {
        return arr + (a)*d2+b;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }
    template <typename T>
    inline void vec4DW(T* arr, int a, int b, int c, int d, int d2, int d3, int d4, T val)
    {
        arr[(((a)*d2+b)*d3+c)*d4+d] = val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }

    template <typename T>
    inline void vec3DW(T* arr, int a, int b, int c, int d2, int d3, T val)
    {
        arr[((a)*d2+b)*d3+c] = val;
    }
    template <typename T>
    inline void vec2DW(T* arr, int a, int b, int d2, T val)
    {
        arr[(a)*d2+b] = val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }

    template <typename T>
    inline void vec4DI(T* arr, int a, int b, int c, int d, int d2, int d3, int d4, T val)
    {
        arr[(((a)*d2+b)*d3+c)*d4+d] += val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }
    template <typename T>
    inline void vec3DI(T* arr, int a, int b, int c, int d2, int d3, T val)
    {
        arr[((a)*d2+b)*d3+c] += val;
    }
    template <typename T>
    inline void vec2DI(T* arr, int a, int b, int d2, T val)
    {
        arr[(a)*d2+b] += val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }

    template <typename T>
    inline void vec4DM(T* arr, int a, int b, int c, int d, int d2, int d3, int d4, T val)
    {
        arr[(((a)*d2+b)*d3+c)*d4+d] *= val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }
    template <typename T>
    inline void vec3DM(T* arr, int a, int b, int c, int d2, int d3, T val)
    {
        arr[((a)*d2+b)*d3+c] *= val;
    }
    template <typename T>
    inline void vec2DM(T* arr, int a, int b, int d2, T val)
    {
        arr[(a)*d2+b] *= val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }
    template <typename T>
    inline void vec4DWA(T* arr, int a, int b, int c, int d, int d2, int d3, int d4, T val)
    {
        #pragma omp atomic write
        arr[(((a)*d2+b)*d3+c)*d4+d] = val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }

    template <typename T>
    inline void vec3DWA(T* arr, int a, int b, int c, int d2, int d3, T val)
    {
        #pragma omp atomic write
        arr[((a)*d2+b)*d3+c] = val;
    }
    template <typename T>
    inline void vec2DWA(T* arr, int a, int b, int d2, T val)
    {
        #pragma omp atomic write
        arr[(a)*d2+b] = val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }

    template <typename T>
    inline void vec4DIA(T* arr, int a, int b, int c, int d, int d2, int d3, int d4, T val)
    {
        #pragma omp atomic update
        arr[(((a)*d2+b)*d3+c)*d4+d] += val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }
    template <typename T>
    inline void vec3DIA(T* arr, int a, int b, int c, int d2, int d3, T val)
    {
        #pragma omp atomic update
        arr[((a)*d2+b)*d3+c] += val;
    }
    template <typename T>
    inline void vec2DIA(T* arr, int a, int b, int d2, T val)
    {
        #pragma omp atomic update
        arr[(a)*d2+b] += val;//(arr.data() + (((a)*d2+b)*d3+c)*d4+d);
    }

    typedef struct Intersection Intersection;
    typedef struct Ray Ray;
    struct Crossing
    {
        int x;
        int z;
    };
    struct Ray
    {
        int index;
        //LinkedList* crossings;
        double intensity;
    };
    
    //Values needed throughout simulation
    extern double maxDev;//stores maximum change in a given iteration, determines convergence
    //extern int beam;//stores which beam is currently being tracked
    extern int iter;//stores iteration number, not currently used
    extern int count;//Used to track a specific grid square along with counter, tracking change over time
    extern int* counter;
    //Launch Ray Values
    extern double cs;
    extern double injected;
    extern int gridcount;
    extern int* intersections;
    extern int ray1num;
    extern double maxInc;
    extern int optimize;
    //Pointers for necessary arrays
    extern short* rayAdjList;//ith level has (nbeams-1-i)*nrays*nrays ints Try short for now, doubtful that over 32,000 rays will be used atm for each beam
    extern Ray* beamIndex;//nbeams nrays
    extern Ray* spatialIndex; //nbeams nx nz
    extern int* marked; //nx nz nbeams Stored as a 1D array, saved SIGNIFICANT amount of time in initialization
    extern double* dedendx; //nx nz
    extern double* dedendz; //nx nz
    extern double* wMult;

    extern int* raypath; //nx nz store single ray path
    extern double* x; //nx
    extern double* z; //nz
    extern double* eden; //nx nz
    extern int* numrays;
    //Marked Stores which rays from each beam have passed through each zone
    //Boxes Stores whi
    extern double* edep; //nx+2 nz+2 nbeams
    extern int* present; //nx nz nbeams
    extern double* machnum; //nx nz
    extern int* boxes; //nbeams nrays ncrossings 2
    extern double* u_flow; //nx nz
    extern double* dkx; //nbeams nrays 2
    extern double* dkz; //nbeams nrays 2
    extern double* dkmag; //nbeams nrays 2

    //Launch_Ray_XZ specific arrays (all have a length of nt)
    //CBET specific arrays
    extern double* W;//nx nz
    extern double* W_init;//nx nz
    extern double* edep_flat;
    extern double* W_new;//nx nz
    extern double* i_b;//nx nz
    extern double* i_b_prev;//nbeams nx nz
    extern double* i_b_new;//nx nz
    extern double* wpe; //nx nz
    extern double* crossesz; //nbeams nrays ncrossings
    extern double* crossesx; //nbeams nrays ncrossings
    extern int* ints; //nbeams nrays ncrossings
    extern int iteration;
    //arrays used only for plotting
    extern double* gain2arr;
    extern double* gain1arr;
    extern double* mag;
    extern double* i_bplot;//nx nz
    extern double* orderplot1; //nx nz
    extern double* orderplot2; //nx nz
    extern double* i_b1Error; //nx nz
    extern double* convergeplot;//nx nz
    extern double* i_b2Error;//nx nz
    extern double* i_b_newplot;//nx nz
    extern double* edenplot; //the array is eden/ncrit,  nx nz
    extern double* edepplot; //nx nz
    extern int* raytrace;//nx nz
    extern double* ib_orig;//nx nz
    extern int* anyInt;//nx nz
    extern double* perturbation;//nx nz
    extern double* mult;//nbeams nrays ncrossings
    //SIMULATION CONSTANTS
    extern int printUpdates;
    extern int printTimings;
    extern int printCBETDiagnostics;
    extern int printRayTrackDiagnostics;
    extern int printHDF5Diagnostics;
    extern int printSpecificTimings;
    extern int pyPlot;
    extern int iterate;
    extern int calcCBET;
    extern int cudaCalc;

    extern int switchvar;
    extern double lambda;
    extern double estat;
    extern double mach;
    extern double Z;
    extern double mi;
    extern double mi_kg;
    extern double Te;
    extern double Te_eV;
    extern double Ti;
    extern double Ti_eV;
    extern double iaw;
    extern double ncrit;
    extern double freq;
    extern double omega;

    //spatial information
    extern int nx;
    extern int nz;
    extern double xmin;
    extern double xmax;
    extern double zmin;
    extern double zmax;
    extern double dx;
    extern double dz;

    extern double maxIncr;
    extern double converge;
    extern int maxIter;
    extern int threads;


    extern int nbeams;
    extern int nrays;
    extern int nt;
    extern int numstored;
    extern int rays_per_zone;
    extern int ncrossings;

    extern double courant_mult;
    extern double intensity;
    extern double offset;
    extern double uray_mult;
    extern double beam_min_z;
    extern double beam_max_z;
    extern double dt;

    //Fundamental Constants
    extern  double sigma; 
    extern  double e0;
    extern  double me; 
    extern  double pi;
    extern  double kb;  //Boltzmann constant in erg/K
    extern  double kb2;   //Boltzmann constant in J/K
    extern  double ec;
    extern  double c;         // Speed of light in cm/s
#endif