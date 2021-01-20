#include "io_interface.hpp"
#include "cuda_help.h"
#include "opencl_help.h"

//GPU IO handlers
int GPUMemCpy(gconfig* gpu, void* hostmem, void* devmem, size_t size, int direction)//copy memory to gpu
{
    switch (gpu->type)
    {
    case CUDA
        CGPUMemCpy(gpu, hostmem, devmem, size, direction);
        break;
    
    case OPENCL
        OCLGPUMemCpy(gpu, hostmem, devmem, size, direction);
        break;
    }
};
int GPUAlloc(gconfig* gpu, void* p, size_t size)//allocate memory on gpu in pointer p
{
        switch (gpu->type)
    {
    case CUDA
        CGPUMalloc(gpu, p, size);
        break;
    
    case OPENCL
        OCLGPUMalloc(gpu, p, size);
        break;
    }
};
