#include "cuda_help.hpp"
#include "implSim.hpp"
#include "dataFields.hpp"


_
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
            double* addr = (double*)entry->getAddress();
            size_t size = entry->getSize();
            //std::cout << entry->getName() + " " << p << std::endl;
            fflush(stdout);
            if(!addr)//returns if arrays are different sizes or if variable n is incorrect, prevents segfaults
            {
                printf("Memory Transfer Error: incorrect size information\n");
                //exit(0);
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
            void* addr = entry->getAddress();//host pointer
            size_t size = entry->getSize();
            if(!addr)//returns if arrays are different sizes or if variable n is incorrect, prevents segfaults
            {

                printf("Memory Transfer Error: incorrect size information: %s %d, %d :: %p %p\n",entry->getName().c_str(),!addr);
                exit(0);
            }
            
        }
    }
    
}
