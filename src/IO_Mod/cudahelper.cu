#include "GPU/cuda_help"


//GPU IO handlers
int CGPUMemCpy(gconfig* gpu, void* dest, void* source, size_t size, int direction)
{
    if(gpu->type != CUDA)
    {
        printf("Non-CUDA GPU called CUDA function\n");
        return 0;
    }
    int dir = (direction) ? cudaMemcpyHostToDevice : cudaMemcpyDeviceToHost;
    cudaMemcpy(dest, source, size, dir);
    return 1;
    //Allocate host memory

};//allocate memory on gpui
int CGPUAlloc(gconfig* gpu, void* p, size_t size)
{
    if(gpu->type != CUDA)
    {
        printf("Non-CUDA GPU called CUDA function\n");
        return 0;
    }
    cudaMalloc(&p, size);
    return 1;
};//allocate memory on gpu in pointer p

