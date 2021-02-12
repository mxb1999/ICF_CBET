#include "GPU/cuda_help.hpp"





//GPU IO handlers
void CGPUMemCpy(void* dest, void* source, size_t size, int direction)
{
    cudaMemcpyKind dir = (direction) ? cudaMemcpyHostToDevice : cudaMemcpyDeviceToHost;
    cudaMemcpy(dest, source, size, dir);
    //Allocate host memory
};//allocate memory on gpui
void CGPUAlloc(void* p, size_t size)
{
    cudaMalloc(&p, size);
};//allocate memory on gpu in pointer p
void CGPUDataTrans(void** hostP, void** deviceP, size_t* sizes, int n, int dir)//function to transfer large quantities of variables from host to GPU, makes bulk transfers easier
{
    if(!hostP || !deviceP || !sizes)
    {
        printf("Null Input Array: hostP: %p, deviceP: %p, sizes: %p", hostP, deviceP, sizes);
        return;
    }
    //for each element in hostP and deviceP, copy data from host to device
    if(dir = 0)
    {
        for(int i = 0; i < n; i++)
        {
            void* currH = hostP[i];//host pointer
            void* currD = deviceP[i];//device pointer
            if(!currH || !currD)//returns if arrays are different sizes or if variable n is incorrect, prevents segfaults
            {
                printf("Memory Transfer Error: incorrect size information\n");
                exit(0);
            }
            cudaMalloc(&currD, sizes[i]);
            cudaMemcpy(currD, currH, sizes[i], cudaMemcpyHostToDevice);
        }
    }else//for each element in hostP and deviceP, copy data from device to host
    {
        for(int i = 0; i < n; i++)
        {
            void* currH = hostP[i];
            void* currD = deviceP[i];
            if(!currH || !currD)//returns if arrays are different sizes or if variable n is incorrect, prevents segfaults
            {
                printf("Memory Transfer Error: incorrect size information\n");
                exit(0);
            }
            cudaMemcpy(currH, currD, sizes[i], cudaMemcpyDeviceToHost);
        }
    }
    
}
