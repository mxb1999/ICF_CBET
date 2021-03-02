#include "parallelConfig.hpp"
#include "cuda_help.hpp"
#include "opencl_help.hpp"

//GPU IO handlers
void GConfig::gpuMalloc(void** p, size_t size)//wrapper function for GPU memory allocation
{
    switch(this->type)
    {
        case CUDA :
            printf("Allocating CUDA Memory\n");
            cudaMalloc(p, size);//default cuda malloc function
            break;
        case OPENCL :
            printf("Allocating OpenCL Memory\n");

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
void GConfig::addData(void* h, void* d,void* r, const std::string name, size_t size, bool isP)
{
    std::map<entryPair>* map = this->addressHash;
    (*map)[name] = new DeviceDataEntry(h,d,r,name,size, isP);//this->addressHash->insert_or_assign(name, new DeviceDataEntry(h,d,name, size));
};
DeviceDataEntry* GConfig::getDataEntry(const std::string name)
{
    printf("name");
    return this->addressHash->at(name);
};
GConfig::GConfig(int t)
{
    this->type = t;
    this->addressHash = new std::map<entryPair>();
}
void GConfig::deviceDataTransfer(int dir)//wrapper function for GPU 
{
    switch(this->type)
    {
        case CUDA :
            CGPUDataTrans(this, dir);//calls custom "mass transfer" function. Assumes that no cuda pointers in devP have already been allocated for Host->Device
            break;
        case OPENCL :
            break;
    }
}
