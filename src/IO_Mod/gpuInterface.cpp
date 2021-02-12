#include "parallelConfig.hpp"
#include "cuda_help.hpp"
#include "opencl_help.hpp"

//GPU IO handlers
void GConfig::gpuMalloc(void** p, size_t size)//wrapper function for GPU memory allocation
{
    switch(this->type)
    {
        case CUDA :
            cudaMalloc(&p, size);//default cuda malloc function
            break;
        case OPENCL :
            break;
    }
};
void GConfig::gpuMemcpy(void* device, void* host, size_t size, int direction)//wrapper function for GPU memcpy functionality
{
    switch(this->type)
    {
        case CUDA :
        {
            cudaMemcpyKind dir = (direction) ? cudaMemcpyHostToDevice : cudaMemcpyDeviceToHost;//1 to copy host to device, 0 for the reverse
            cudaMemcpy(device, host, size, dir);//default cuda memcpy function
            break;
        }
        case OPENCL :
            break;
    }
}

void GConfig::deviceDataTransfer(void** hostP, void** devP, size_t* sizes, int n, int dir)//wrapper function for GPU 
{
    switch(this->type)
    {
        case CUDA :
            CGPUDataTrans(hostP, devP, sizes, n, dir);//calls custom "mass transfer" function. Assumes that no cuda pointers in devP have already been allocated for Host->Device
            break;
        case OPENCL :
            break;
    }
}
