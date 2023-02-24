#ifndef VARINITCU
#define VARINITCU
    #include "cuda_help.hpp"
    #define DSIZE sizeof(double)
    #define ISIZE sizeof(int)


    __device__ const int printUpdates_cu=1;
    __device__ const int printTimings_cu=0;
    __device__ const int printCBETDiagnostics_cu=1;
    __device__ const int printRayTrackDiagnostics_cu=0;
    __device__ const int printHDF5Diagnostics_cu=0;
    __device__ const int printSpecificTimings_cu=0;
    __device__ const int iterate_cu=0;
    __device__ const int calcCBET_cu=0;
    __device__ const double lambda_cu=3.51e-5;
    __device__ const double estat_cu=4.80320427e-10;
    __device__ const double mach_cu=-1.414213562;
    __device__ const double Z_cu=3.1;
    __device__ const double Te_cu=23209010.4;
    __device__ const double Te_eV_cu=2.0e3;
    __device__ const double Ti_cu=11604505.2;
    __device__ const double Ti_eV_cu=2.0e3;
    __device__ const double iaw_cu=0.2;
    __device__ const int nx_cu=50;
    __device__ const int nz_cu=50;
    __device__ const double xmin_cu=-5.0e-4;
    __device__ const double xmax_cu=5.0e-4;
    __device__ const double zmin_cu=-5.0e-4;
    __device__ const double zmax_cu=5.0e-4;
    __device__ const double maxIncr_cu=0.2;
    __device__ const int maxIterations_cu=100;
    __device__ const double converge_cu=1e-4;
    __device__ const int threads_cu=12;
    __device__ const int nbeams_cu=2;
    __device__ const int rays_per_zone_cu=5;
    __device__ const double courant_mult_cu=0.05;
    __device__ const double intensity_cu=2.0e15;
    __device__ const double offset_cu=0;
    __device__ const double beam_max_z_cu=3.0e-4;
    __device__ const double beam_min_z_cu=-3.0e-4;

    __device__ double c_cu, lambda_cu,  estat_cu,mach_cu,Z_cu,mi_cu,mi_kg_cu,Te_cu,Te_eV_cu,Ti_cu,Ti_eV_cu,iaw_cu,ncrit_cu,freq_cu,omega_cu;
    //spatial information
    __device__ int nx_cu, nz_cu;
    __device__ double xmin_cu, xmax_cu, zmin_cu, zmax_cu, dx_cu, dz_cu;

    __device__ double maxIncr_cu, converge_cu;

    __device__ int nbeams_cu, nrays_cu, nt_cu, numstored_cu, rays_per_zone_cu, ncrossings_cu;

    __device__ double intensity_cu, offset_cu, uray_mult_cu, beam_min_z_cu, beam_max_z_cu, dt_cu, courant_mult_cu;
    /*
    __device__ double* edep_cu, *crossx_cu, *crossz_cu;

    //store device pointers to copy to in an array
    std::vector<void*> devP = {&c_cu, &lambda_cu, &estat_cu, &mach_cu, &Z_cu, &mi_cu, &mi_kg_cu, &Te_cu, &Te_eV_cu, &Ti_cu, &Ti_eV_cu, &iaw_cu, &ncrit_cu, &freq_cu, &omega_cu
        , &xmin_cu, &xmax_cu, &zmin_cu, &zmax_cu, &dx_cu, &dz_cu, &intensity_cu, &offset_cu, &uray_mult_cu, &beam_min_z_cu, &beam_max_z_cu, &dt_cu, &courant_mult_cu
        , &nx_cu, &nz_cu, &nbeams_cu, &nrays_cu, &nt_cu, &numstored_cu, &rays_per_zone_cu, &ncrossings_cu};
    //store pointers for data to be copied
    std::vector<void*> hostP =  {&c, &lambda, &estat, &mach, &Z, &mi, &mi_kg, &Te, &Te_eV, &Ti, &Ti_eV, &iaw, &ncrit, &freq, &omega//double const values
        , &xmin, &xmax, &zmin, &zmax, &dx, &dz, &intensity, &offset, &uray_mult, &beam_min_z, &beam_max_z, &dt, &courant_mult//more double const values
        , &nx, &nz, &nbeams, &nrays, &nt, &numstored, &rays_per_zone, &ncrossings};//
    std::vector<std::string> varNames =  {"c","lambda","estat","mach","Z","mi","mi_kg","Te","Te_eV","Ti","Ti_eV","iaw","ncrit","freq","omega"//double const values
    ,"xmin","xmax","zmin","zmax","dx","dz","intensity","offset","uray_mult","beam_min_z","beam_max_z","dt","courant_mult"//more double const values
    ,"nx","nz","nbeams","nrays","nt","numstored","rays_per_zone","ncrossings"};//
    std::vector<size_t> sizes = { DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE, DSIZE, //c lambda, estat ,mach, Z, mi, mi_kg, Te, Te_eV, Ti, Ti_eV, iaw, ncrit, freq, omega
        DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE, //xmin xmax zmin zmax dx dz intensity offset uray_mult beam_min_z beam_max_z, dt, courant_mult
        ISIZE,ISIZE,ISIZE,ISIZE,ISIZE,ISIZE,ISIZE,ISIZE};//nx nz nbeams nrays nt numstored, rays_per_zone, ncrossings
    std::vector<bool> pointers = {true,
        true, true, true,true,true, true,//x z edep crossesx crossesz boxes
        true, true, true, true,true, //eden wpe dedendx dedendz ints
        false,false,false,false,false,false,false,false,false,false,false,false,false,false, false, //c lambda, estat ,mach, Z, mi, mi_kg, Te, Te_eV, Ti, Ti_eV, iaw, ncrit, freq, omega
        false,false,false,false,false,false,false,false,false,false,false,false,false, //xmin xmax zmin zmax dx dz intensity offset uray_mult beam_min_z beam_max_z, dt, courant_mult
        false,false,false,false,false,false,false,false};


*/

#endif