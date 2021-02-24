#include "cuda_var_init.hpp"
__device__ int test_cu;

void gpuInit()
{
    size_t sizeStore[] = {sizes[0],nx*DSIZE, nz*DSIZE, (nx+2)*(nz+2)*nbeams*DSIZE,CROSS*DSIZE,CROSS*DSIZE, CROSS*2*ISIZE,//x z edep crossesx crossesz boxes
    GRID*DSIZE, GRID*DSIZE, GRID*DSIZE, GRID*DSIZE,CROSS*numstored*ISIZE, //eden wpe dedendx dedendz ints
    DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE, DSIZE, //c lambda, estat ,mach, Z, mi, mi_kg, Te, Te_eV, Ti, Ti_eV, iaw, ncrit, freq, omega
    DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE,DSIZE, //xmin xmax zmin zmax dx dz intensity offset uray_mult beam_min_z beam_max_z, dt, courant_mult
    ISIZE,ISIZE,ISIZE,ISIZE,ISIZE,ISIZE,ISIZE,ISIZE};//nx nz nbeams nrays nt numstored, rays_per_zone, ncrossings
    printf("Initializing GPU Pointers\n");
    deviceConfiguration = new GConfig(CUDA);
    int len = varNames.size();
    printf("%d %d %d %d %d\n", hostP.size(), 0, varNames.size(), sizes.size(), pointers.size());
    for(int i = 0; i < len; i++)
    {
        printf("size %d\n", sizeStore[i]);
        if(pointers[i])
        {
            deviceConfiguration->addData(hostP[i], devP[i], varNames[i], sizeStore[i], pointers[i]);
        }else
        {
            deviceConfiguration->addData(hostP[i], NULL, varNames[i], sizeStore[i], pointers[i]);
        }
    }
    int* test;
    printf("Attempt to Copy\n");
    cudaError_t error = cudaGetSymbolAddress((void**)&test, "test_cu");
    //cudaMemcpyToSymbol(name, &nx, sizeof(int), 0, cudaMemcpyHostToDevice);
    printf("%s\n", cudaGetErrorString(error));
}