#include "CBET_Interface.hpp"
#include "parallelConfig.hpp"
#include "cuda_help.hpp"

void cbetOptimize()
{
    if(cudaCalc)
    {
        cudaMemAdvise(machnum, sizeof(double)*GRID, cudaMemAdviseSetReadMostly, 0);
        cudaMemAdvise(u_flow, sizeof(double)*GRID, cudaMemAdviseSetReadMostly, 0);

        cudaMemAdvise(dk, sizeof(double)*CROSS*3, cudaMemAdviseSetReadMostly, 0);
        cudaMemAdvise(dkmag, sizeof(double)*CROSS, cudaMemAdviseSetReadMostly, 0);

        cudaMemAdvise(i_b_new, sizeof(double)*CROSS, cudaMemAdviseSetPreferredLocation, 0);
        cudaMemAdvise(i_b, sizeof(double)*CROSS, cudaMemAdviseSetPreferredLocation, 0);


        cudaMemPrefetchAsync(machnum, sizeof(double)*GRID, 0);
        cudaMemPrefetchAsync(u_flow, sizeof(double)*GRID, 0);

        cudaMemPrefetchAsync(dk, sizeof(double)*CROSS*3, 0);
        cudaMemPrefetchAsync(dkmag, sizeof(double)*CROSS, 0);

        cudaMemPrefetchAsync(i_b_new, sizeof(double)*CROSS, 0);
        cudaMemPrefetchAsync(i_b, sizeof(double)*CROSS, 0);
    }
}


void freeCBETArrs()
{

  if(cudaCalc)
  {
    cudaError_t err = cudaFree(i_b_new);
    err = cudaFree(i_b);
    err = cudaFree(machnum);
    err = cudaFree(u_flow);
    err = cudaFree(dk);
    err = cudaFree(dkmag);
  }else
  {
    delete [] i_b_new;
    delete [] i_b;
    delete [] machnum;
    delete [] u_flow;
    delete [] dk;
    delete [] dkmag;
  }

}