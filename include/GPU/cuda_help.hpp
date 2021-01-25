#include "implSim.hpp"
#include "io_interface.hpp"
//define device variables
//define other functions necessary for GPU operation
#ifndef CUDAHELP
#define CUDAHELP

    __device__ double* x_cu;
    __device__ double* z_cu;
    __device__ bool* marked_cu;//nbeams nrays nx nz
    __device__ double lambda,estat,mach,Z,mi,mi_kg,Te,Te_eV,Ti,Ti_eV,iaw,ncrit,freq,omega;
  //spatial information
    __device__ int nx, nz;
    __device__ double xmin, xmax, zmin, zmax, dx, dz;

    __device__ double maxIncr, converge;
    __device__ int threads,maxIter;

    __device__ int nbeams, nrays, nt, numstored, rays_per_zone, ncrossings;

    __device__ double intensity, offset, uray_mult, beam_min_z, beam_max_z, dt, courant_mult;
#endif