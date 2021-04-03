#include <stdio.h>
#include <stdlib.h>
#include "Trace_interface.hpp"
#include "io_interface.hpp"
#include "cuda_help.hpp"
#include <stdarg.h>

void optimizePreTraceArrs()
{
  if(cudaCalc)
  {
    cudaMemAdvise(dedendx, sizeof(double)*GRID, cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(dedendz, sizeof(double)*GRID, cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(x, sizeof(double)*nx, cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(z, sizeof(double)*nz, cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(eden, sizeof(double)*GRID, cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(marked, sizeof(int)*GRID*numstored*nbeams, cudaMemAdviseSetPreferredLocation, 0);
    cudaMemAdvise(boxes, sizeof(int)*CROSS*2, cudaMemAdviseSetPreferredLocation, 0);
    cudaMemAdvise(present, sizeof(int)*GRID*nbeams, cudaMemAdviseSetPreferredLocation, 0);
    cudaMemAdvise(crossesx, sizeof(double)*CROSS, cudaMemAdviseSetPreferredLocation, 0);
    cudaMemAdvise(crossesz, sizeof(double)*CROSS, cudaMemAdviseSetPreferredLocation, 0);
  }
}
void optimizePostTraceArrs()
{
  if(cudaCalc)
  {
    cudaMemAdvise(marked, sizeof(int)*GRID*numstored*nbeams, cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(boxes, sizeof(int)*CROSS*2, cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(present, sizeof(int)*GRID*nbeams, cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(crossesx, sizeof(double)*CROSS, cudaMemAdviseSetReadMostly, 0);
    cudaMemAdvise(crossesz, sizeof(double)*CROSS, cudaMemAdviseSetReadMostly, 0);
  }
}

void freeTraceArrs()
{
  if(cudaCalc)
  {
    cudaFree(dedendx);
    cudaFree(dedendz);
    cudaFree(x);
    cudaFree(z);
    cudaFree(eden);
    cudaFree(marked);
    cudaFree(present);
    cudaFree(boxes);
    cudaFree(wpe);
    cudaFree(crossesx);
    cudaFree(crossesz);
  }else
  {
    delete [] dedendx;
    delete [] dedendz;
    delete [] x;
    delete [] z;
    delete [] eden;
    delete [] marked;
    delete [] present;
    delete [] boxes;
    delete [] wpe;
    delete [] crossesx;
    delete [] crossesz;
  }
}