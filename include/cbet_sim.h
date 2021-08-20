/*Variation of Follett Algorithm utilizing ray-space problem decomposition

Struct: Grid
    -EDen
    -MachNum
    -nx, ny, nz
    -dx, dy, dz
    -Courant #

Arrays of interest:
    -Multiplier (nrays*ncrossings, initialized to all 1's)
    -Initial intensity (nrays)
    -Wave vector (nrays*ncrossings)
    -Ray areas (nrays*ncrossings, may be able to store initial as double precision and subsequent values as a single precision multiplier)
    -Boxes (nrays*ncrossings, can typedef)

Need to track convergence, can use a device function to set a value for each recurrance->stored in register

1. initialize EDen and Machnum
2. Decompose rays across devices
3. Perform ray trace -> fill boxes
4. Perform Divide and Conquer CBET
    4a: Divide-> If greater than 2, divide seed range by 2 and the pump range by 2, call recurrance function 4 times.
    4b: Conquer-> Check number of rays being compared. If 2, then base case, find intersection point and calculate multiplier
    4c. Combine-> For all pump rays considered, multiply over the ray-crossing multiplier space
    4d. Continue until convergence
*/
#include <stdio.h>
#include <stdlib.h>
#ifndef CBETSIM
    #define CBETSIM
    typedef double spatial_t;
    typedef double plasma_t;
    typedef double laser_t;

    typedef struct
    {
        spatial_t *eden, *machnum, *wpe, *beam_centers, *beam_vectors, dims[6];
        spatial_t dx, dy, dz;
        int nx, ny, nz, ncrossings, nrays, nbeams, ndims, rays_per_zone;
        spatial_t courant;
        plasma_t Te, Ti, Z, iaw, ncrit, mi, A, cs;
        laser_t intensity, lambda, omega, *beam_radii;
    }Grid;
    #define ACCESS3D(arr, i, j, k, s2, s3) arr[((i)*(s2)+(j))*(s3)+(k)]
    #define ACCESS2D(arr, i, j, s2) arr[(i)*(s2)+(j)]
    #define GRID grid
    #define NX GRID->nx
    #define NZ GRID->nz
    #define NY GRID->ny
    #define NRAYS GRID->nrays
    #define RAYS_PER_ZONE GRID->rays_per_zone
    #define NBEAMS GRID->nbeams
    #define CROSSINGS GRID->ncrossings
    #define DX GRID->dx
    #define DY GRID->dy
    #define DZ GRID->dz
    #define XMAX ACCESS2D(GRID->dims, 0, 0, 2)
    #define XMIN ACCESS2D(GRID->dims, 0, 1, 2)
    #define YMAX ACCESS2D(GRID->dims, 1, 0, 2)
    #define YMIN ACCESS2D(GRID->dims, 1, 1, 2)
    #define ZMAX ACCESS2D(GRID->dims, 2, 0, 2)
    #define ZMIN ACCESS2D(GRID->dims, 2, 1, 2)
    #define OMEGA GRID->omega
    #define COURANT GRID->courant
    #define T_E GRID->Te
    #define T_I GRID->Ti
    #define Z_VAL GRID->Z
    #define IAW GRID->iaw
    #define A_VAL GRID->A
    #define NCRIT GRID->ncrit
    #define MI GRID->mi
    #define CS GRID->cs
    //Fundamental Constants
    const double sigma = 1.7e-4;
    const double e0 =8.85418782e-12;
    const double me =9.10938291e-31;
    const double pi =3.14159265359;
    const double kb= 1.3806485279e-16;   //Boltzmann constant in erg/K
    const double kb2= 1.3806485279e-23;   //Boltzmann constant in J/K
    const double ec= 1.60217662e-19;
    const double c= 29979245800.0;              // Speed of light in cm/s
    #define NDIMS GRID->ndims
    #define PEAK_INTENSITY GRID->intensity
    #define BEAM_RAD(beam) GRID->beam_radii[beam]
    #define BEAM_CENTER(beam, dim) ACCESS2D(GRID->beam_min, beam, dim, NDIMS)
    #define EDEN(i, j, k) ACCESS3D(GRID->eden, i, j, k, NY, NZ)
    #define MACH(i, j, k) ACCESS3D(GRID->machnum, i, j, k, NY, NZ)
    #define WPE(i, j, k) ACCESS3D(GRID->wpe, i, j, k, NY, NZ)
    #define EDEN_2D(i, j) ACCESS2D(GRID->eden, i, j, NZ)
    #define MACH_2D(i, j) ACCESS2D(GRID->machnum, i, j, NZ)
    #define WPE_2D(i, j) ACCESS2D(GRID->wpe, i, j, NZ)
    //Boxes, only need 3 bits per crossing in 2D, 5 bits in 2D
    //store 4 crossings
    typedef union
    {
        char path;
    }path_t;

    #define UNCHANGED 0x0
    #define INCREASED 0x1
    #define DECREASED 0x2

    #define XMASK 0x03
    #define YMASK 0x0C
    #define ZMASK 0x30

    #define X(path) (path & XMASK)
    #define Y(path) (path & YMASK) >> 2
    #define Z(path) (path & ZMASK) >> 4
    #define IFNOMEM(pointer){\
        if(pointer == NULL)\
        {\
            printf("No free memory for allocation");\
            return NULL;\
        }\
    };
    Grid* setupGrid();
#endif


