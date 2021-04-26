#include "GPU/cuda_help.hpp"
#include "implSim.hpp"


__device__ LinkCross* new_LinkCross(double* values, unsigned int ix, unsigned int iz, unsigned int cross, int nz)//creates a new link in a "crossing" chain, keeps track of 
{
    LinkCross* result;
    result = (LinkCross*)malloc(sizeof(LinkCross));
    if(!result)
    {
        return NULL;
    }
    if(values != NULL)
    {
        for(int i = 0; i < 9;i++)
        {
            result->vals[i] = values[i]; 
        }
    }
    result->location = ix*nz+iz;
    result->cross = cross;
    result->next = NULL;
    return result;
};
LinkCross* new_LinkCrossHost(double* values, unsigned int ix, unsigned int iz, unsigned int cross)//creates a new link in a "crossing" chain, keeps track of 
{
    LinkCross* result;
    if(cudaCalc)
    {
        cudaError_t err = cudaMallocManaged(&result, sizeof(LinkCross));
        if(err)
        {
            printf("%s\n", cudaGetErrorString(err));
        }
    }else
    {
        result = new LinkCross;
    }
    if(!result)
    {
        return NULL;
    }
    if(values != NULL)
    {
        for(int i = 0; i < 9;i++)
        {
            result->vals[i] = values[i]; 
        }
    }
    result->location = ix*nz+iz;
    result->cross = cross;
    result->next = NULL;
    return result;
};
__device__ int add_LinkCrossAtomic(LinkCross* entry, LinkCross* newVal)//only adds if current pointer is null, uses CUDA atomicCAS
{
    if(entry == NULL)
    {
        return 0;
    }
    unsigned long long int result = atomicCAS((unsigned long long int*)entry, 0x0, (unsigned long long int)newVal);
    if(result != 0x0)
    {
        return 0;
    }
    return 1;
};
__device__ void add_LinkCross(LinkCross* entry, LinkCross* newVal)//only adds if current pointer is null, uses CUDA atomicCAS
{
    if(entry == NULL)
    {
        return;
    }
    entry->next = newVal;
};
void add_LinkCrossHost(LinkCross* entry, LinkCross* newVal)//only adds if current pointer is null, uses CUDA atomicCAS
{
    if(entry == NULL)
    {
        return;
    }
    entry->next = newVal;
};
TrackArrs* deviceTrackArrs(int device)
{
  TrackArrs* constArrs;
  cudaMallocManaged(&constArrs, sizeof(TrackArrs));
  if(constArrs == NULL)
  {
    return NULL;
  }
  constArrs->dedendx_cu = dedendx;
  constArrs->dedendz_cu = dedendz;
  constArrs->x_cu = x;
  constArrs->z_cu = z;
  constArrs->crossesx_cu = crossesx;
  constArrs->crossesz_cu = crossesz;
  constArrs->edep_cu = edep;
  constArrs->wpe_cu = wpe;
  constArrs->boxes_cu = boxes;
  constArrs->ints_cu = ints;

  /*cudaError_t error = cudaMemPrefetchAsync(dedendx, sizeof(double)*nx*nz, device);
  cudaMemPrefetchAsync(dedendz, sizeof(double)*GRID, device);
  cudaMemPrefetchAsync(edep, sizeof(double)*(nx+2)*(nz+2)*RAYS, device);
  cudaMemPrefetchAsync(wpe, sizeof(double)*GRID, device);
  printf("%s\n", cudaGetErrorString(error));
  cudaMemAdvise(x, sizeof(int)*nx, cudaMemAdviseSetReadMostly, device);
  cudaMemAdvise(z, sizeof(int)*nz, cudaMemAdviseSetReadMostly, device);
  cudaMemAdvise(dedendx, sizeof(double)*GRID, cudaMemAdviseSetReadMostly, device);
  cudaMemAdvise(dedendz, sizeof(double)*GRID, cudaMemAdviseSetReadMostly, device);
  cudaMemAdvise(edep, sizeof(double)*(nx+2)*(nz+2)*RAYS, cudaMemAdviseSetPreferredLocation, device);
  cudaMemAdvise(wpe, sizeof(double)*GRID, cudaMemAdviseSetReadMostly, device);

  cudaMemAdvise(constArrs, sizeof(TrackArrs), cudaMemAdviseSetReadMostly, device);*/
  return constArrs;
}
TrackConst* deviceTrackConst(int device)
{
  TrackConst* constVals;
  cudaMallocManaged(&constVals, sizeof(TrackConst));
  if(constVals == NULL)
  {
    return NULL;
  }
  constVals->dx_cu = dx;
  constVals->dz_cu = dz;
  constVals->nx_cu = nx;
  constVals->nz_cu = nz;
  constVals->xmax_cu = xmax;
  constVals->xmin_cu = xmin;
  constVals->zmax_cu = zmax;
  constVals->zmin_cu = zmin;
  constVals->nrays_cu = nrays;
  constVals->nbeams_cu = nbeams;

  constVals->ncrossings_cu = ncrossings;
  constVals->numstored_cu = numstored;
  constVals->nt_cu = nt;
  constVals->dt_cu = dt;
  constVals->omega_cu = omega;
  constVals->ncrit_cu = ncrit;

  constVals->c_cu = c;
  //cudaMemAdvise(constVals, sizeof(TrackConst), cudaMemAdviseSetReadMostly, device);
  return constVals;
}

/*
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
            double* ref = entry->getRefAddress();
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
                //cudaMalloc(&currD, size);
                cudaError_t error = cudaMemcpyToSymbol(*ref, currH, size,0, cudaMemcpyHostToDevice);
                printf("Malloc Error: %s\n", cudaGetErrorString(error));

            }else
            {
                cudaError_t error = cudaMemcpyToSymbol(*ref,currH,size, 0, cudaMemcpyHostToDevice);
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
*/