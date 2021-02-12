#ifndef VARINITCU
#define VARINITCU
    #include "cuda_help.hpp"
    __device__ void* rays_cu;//rayinit struct, to be declared later
    __device__ double* x_cu;
    __device__ double* z_cu;
    __device__ double* eden_cu;
    __device__ double* wpe_cu;
    __device__ double* dedendx_cu;
    __device__ double* dedendz_cu;
    __device__ bool* marked_cu;//nbeams nrays nx nz
    __device__ int* boxes_cu;
    __device__ double c_cu, lambda_cu, __device__ estat_cu,mach_cu,Z_cu,mi_cu,mi_kg_cu,Te_cu,Te_eV_cu,Ti_cu,Ti_eV_cu,iaw_cu,ncrit_cu,freq_cu,omega_cu;
    //spatial information
    __device__ int nx_cu, nz_cu;
    __device__ double xmin_cu, xmax_cu, zmin_cu, zmax_cu, dx_cu, dz_cu;

    __device__ double maxIncr_cu, converge_cu;

    __device__ int nbeams_cu, nrays_cu, nt_cu, numstored_cu, rays_per_zone_cu, ncrossings_cu;

    __device__ double intensity_cu, offset_cu, uray_mult_cu, beam_min_z_cu, beam_max_z_cu, dt_cu, courant_mult_cu;
    __device__ double* edep_cu, *crossx_cu, *crossz_cu;
    
    //store device pointers to copy to in an array
    void* devP[] = {&rays_cu, &x_cu, &z_cu, &edep_cu, &crossx_cu, &crossz_cu, &boxes_cu,&eden_cu, &wpe_cu, &dedendx_cu, &dedendz_cu
        , &c_cu, &lambda_cu, &estat_cu, &mach_cu, &Z_cu, &mi_cu, &mi_kg_cu, &Te_cu, &Te_eV_cu, &Ti_cu, &Ti_eV_cu, &iaw_cu, &ncrit_cu, &freq_cu, &omega_cu
        , &xmin_cu, &xmax_cu, &zmin_cu, &zmax_cu, &dx_cu, &dz_cu, &intensity, &offset, &uray_mult, &beam_min_z, &beam_max_z, &dt, &courant_mult
        , &nx_cu, &nz_cu, &nbeams_cu, &nrays_cu, &nt_cu, &numstored_cu, &rays_per_zone_cu, &ncrossings_cu};
    //store pointers for data to be copied
    void* hostP[] =  {NULL, //space for rayinit array, to be filled in later
        &x, &z, &edep, &crossesx, &crossesz, &boxes, &eden, &wpe, &dedendx, &dedendz //double arrays
        , &c, &lambda, &estat, &mach, &Z, &mi, &mi_kg, &Te, &Te_eV, &Ti, &Ti_eV, &iaw, &ncrit, &freq, &omega//double const values
        , &xmin, &xmax, &zmin, &zmax, &dx, &dz, &intensity, &offset, &uray_mult, &beam_min_z, &beam_max_z, &dt, &courant_mult//more double const values
        , &nx, &nz, &nbeams, &nrays, &nt, &numstored, &rays_per_zone, &ncrossings};//
    size_t sizes[] = {0,
        nx*sizeof(double), nz*sizeof(double), (nx+2)*(nz+2)*nbeams*sizeof(double),(nrays)*(ncrossings)*nbeams*sizeof(double),(nrays)*(ncrossings)*nbeams*sizeof(double), nbeams*ncrossings*nrays*2*sizeof(int),
        sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),
        sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),sizeof(double),
        sizeof(int),sizeof(int),sizeof(int),sizeof(int),sizeof(int),sizeof(int),sizeof(int),sizeof(int)};
#endif