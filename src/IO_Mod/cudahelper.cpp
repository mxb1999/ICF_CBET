#include "GPU/cuda_help.h"


//GPU IO handlers
int GConfig::CGPUMemCpy(void* dest, void* source, size_t size, int direction)
{
    if(this.type != CUDA)
    {
        printf("Non-CUDA GPU called CUDA function\n");
        return 0;
    }
    int dir = (direction) ? cudaMemcpyHostToDevice : cudaMemcpyDeviceToHost;
    cudaMemcpy(dest, source, size, dir);
    return 1;
    //Allocate host memory

};//allocate memory on gpui
int GConfig::CGPUAlloc(void* p, size_t size)
{
    if(this.type != CUDA)
    {
        printf("Non-CUDA GPU called CUDA function\n");
        return 0;
    }
    cudaMalloc(&p, size);
    return 1;
};//allocate memory on gpu in pointer p
void GConfig::CGPUDataTrans(void** hostP, void** deviceP, size_t* sizes, int n, int dir)//function to transfer large quantities of variables from host to GPU, makes bulk transfers easier
{
    if(this.type != CUDA)
    {
        printf("Non-CUDA GPU called CUDA function\n");
        return;
    }
    if(!hostP || !deviceP || !sizes)
    {
        printf("Null Input Array: hostP: %p, deviceP: %p, sizes: %p", hostP, deviceP, sizes);
        return;
    }
    //for each element in hostP and deviceP, copy data from host to device
    if(dir = cudaMemcpyHostToDevice)
    {
        for(int i = 0; i < n; i++)
        {
            void* currH = hostP[i];//host pointer
            void* currD = devP[i];//device pointer
            if(!currH || !currP)//returns if arrays are different sizes or if variable n is incorrect, prevents segfaults
            {
                printf("Memory Transfer Error: incorrect size information\n");
                exit(0);
            }
            cudaMalloc(&currD, sizes[i]);
            cudaMemcpy(currD, currH, sizes[i], dir);
        }
    }else//for each element in hostP and deviceP, copy data from device to host
    {
        for(int i = 0; i < n; i++)
        {
            void* currH = hostP[i];
            void* currP = devP[i];
            if(!currH || !currP)//returns if arrays are different sizes or if variable n is incorrect, prevents segfaults
            {
                printf("Memory Transfer Error: incorrect size information\n");
                exit(0);
            }
            cudaMemcpy(currH, currD, sizes[i], dir);
        }
    }
    
}
