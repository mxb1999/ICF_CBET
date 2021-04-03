#include "CBET_Interface.hpp"
#include "parallelConfig.hpp"

void cbetOptimize()
{
    if(cudaCalc)
    {
        cudaMemAdvise(machnum, sizeof(double)*GRID, cudaMemAdviseSetReadMostly, 0);
        cudaMemAdvise(u_flow, sizeof(double)*GRID, cudaMemAdviseSetReadMostly, 0);

        cudaMemAdvise(dkx, sizeof(double)*CROSS, cudaMemAdviseSetReadMostly, 0);
        cudaMemAdvise(dkz, sizeof(double)*CROSS, cudaMemAdviseSetReadMostly, 0);
        cudaMemAdvise(dkmag, sizeof(double)*CROSS, cudaMemAdviseSetReadMostly, 0);

        cudaMemAdvise(i_b_new, sizeof(double)*CROSS, cudaMemAdviseSetPreferredLocation, 0);
        cudaMemAdvise(i_b, sizeof(double)*CROSS, cudaMemAdviseSetPreferredLocation, 0);


        cudaMemPrefetchAsync(machnum, sizeof(double)*GRID, 0);
        cudaMemPrefetchAsync(u_flow, sizeof(double)*GRID, 0);

        cudaMemPrefetchAsync(dkx, sizeof(double)*CROSS, 0);
        cudaMemPrefetchAsync(dkz, sizeof(double)*CROSS, 0);
        cudaMemPrefetchAsync(dkmag, sizeof(double)*CROSS, 0);

        cudaMemPrefetchAsync(i_b_new, sizeof(double)*CROSS, 0);
        cudaMemPrefetchAsync(i_b, sizeof(double)*CROSS, 0);
    }
}


void freeCBETArrs()
{
  if(cudaCalc)
  {
    cudaFree(i_b_new);
    cudaFree(i_b);
    cudaFree(machnum);
    cudaFree(u_flow);
    cudaFree(dkx);
    cudaFree(dkz);
    cudaFree(dkmag);
  }else
  {
    delete [] i_b_new;
    delete [] i_b;
    delete [] machnum;
    delete [] u_flow;
    delete [] dkx;
    delete [] dkz;
    delete [] dkmag;
  }
}