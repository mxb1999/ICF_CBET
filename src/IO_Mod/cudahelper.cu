#include "cuda_help.hpp"





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
void CGPUDataTrans(GConfig* device, int dir)//function to transfer large quantities of variables from host to GPU, makes bulk transfers easier
{
    printf("Check 3\n");
    fflush(stdout);
    std::map<entryPair> *map = device->getMap();
    //for each element in hostP and deviceP, copy data from host to device
    if(dir == 0)
    {
        
        std::map<entryPair>::const_iterator it;
        for(it = map->begin(); it != map->end(); ++it)
        {
            
            DeviceDataEntry* entry = it->second;
            void* currH = entry->getHostAddress();//host pointer
            void* currD = entry->getDeviceAddress();//device pointer
            size_t size = entry->getSize();
            bool p = entry->isPointer();
            //std::cout << entry->getName() + " " << p << std::endl;
            fflush(stdout);
            if(!currH)//returns if arrays are different sizes or if variable n is incorrect, prevents segfaults
            {
                printf("Memory Transfer Error: incorrect size information\n");
                //exit(0);
            }
            if(p)
            {
                std::string deviceName = entry->getName()+"_cu";
                
                //cudaMalloc(&currD, size);
                cudaError_t error = cudaMemcpyToSymbol(deviceName.c_str(), currH, size,0, cudaMemcpyHostToDevice);
                printf("Malloc Error: %s\n", cudaGetErrorString(error));

            }else
            {
                std::string deviceName = entry->getName()+"_cu";
                cudaError_t error = cudaMemcpyToSymbol(deviceName.c_str(),currH,size, 0, cudaMemcpyHostToDevice);
                printf("Malloc Error: %s\n", cudaGetErrorString(error));

                fflush(stdout);
            }

        }
    }else//for each element in hostP and deviceP, copy data from device to host
    {
        std::map<entryPair>::const_iterator it;
        for(it = map->begin(); it != map->end(); ++it)
        {

            DeviceDataEntry* entry = it->second;
            //printf("%s\n", entry->getName().c_str());
            //fflush(stdout);
            void* currH = entry->getHostAddress();//host pointer
            void* currD = entry->getDeviceAddress();//device pointer
            size_t size = entry->getSize();
            bool p = entry->isPointer();
            if(!currH || !currD)//returns if arrays are different sizes or if variable n is incorrect, prevents segfaults
            {

                printf("Memory Transfer Error: incorrect size information: %s %d, %d :: %p %p\n",entry->getName().c_str(),!currH,!currD, currH, currD);
                exit(0);
            }
            if(p)//only copy back for pointers
            {
                int q = cudaMemcpy(currH, currD, size, cudaMemcpyDeviceToHost);
                if(q)
                {
                    printf("Memory Transfer Error: CUDA Memcpy did not like that: %s %ld,:: %p %p\n",entry->getName().c_str(),size, currH, currD);
                }
            }
            
        }
    }
    
}
