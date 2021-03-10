#include "CBET_Interface.hpp"
#include "parallelConfig.hpp"

void initArrays()
{
    //Allocate global memory for relevant arrays
    cudaMallocManaged(&i_b, sizeof(double)*GRID);
    cudaMallocManaged(&i_b_prev, sizeof(double)*GRID);
    cudaMallocManaged(&i_b_new, sizeof(double)*GRID);
    cudaMallocManaged(&W_new, sizeof(double)*GRID);
    cudaMallocManaged(&W, sizeof(double)*GRID);

      //Initialize CBET array values
  for(int m = 0; m < nbeams; m++)
  {
    #pragma omp parallel for num_threads(threads)
    for(int i = 0; i < nx;i++)
    {
      for(int j = 0; j < nz;j++)
      {
        double energyDep = vec3D(edep_flat,m,i,j,nx,nz);//break up pipeline to prevent RAW dependencies
        double initEnergy = sqrt(1.0-vec2D(eden,i,j,nz)/ncrit)/double(rays_per_zone);
        double flownum = vec2D(machnum,i,j,nz)*cs;
        vec3DW(i_b,m,i,j,nx,nz, energyDep);///present[m][i][j];
        vec3DW(i_b_new,m,i,j,nx,nz, energyDep);///present[m][i][j];
        vec3DW(W,m,i,j,nx,nz, initEnergy);///present[m][i][j];
        vec3DW(W_new,m,i,j,nx,nz, initEnergy);///present[m][i][j];
        vec2DW(u_flow,i,j,nz, flownum); //
      }

    }
  }
  
    for(int i = 0; i < nbeams;i++)
    {
        #pragma omp parallel for num_threads(threads)
        for(int j = 0; j < nrays; j++)
        {
            for(int m = 0; m < ncrossings-1; m++)
            {
                double delX = vec3D(crossesx,i,j,m+1,nrays,ncrossings)-vec3D(crossesx,i,j,m,nrays,ncrossings);
                double delZ = vec3D(crossesz,i,j,m+1,nrays,ncrossings)-vec3D(crossesz,i,j,m,nrays,ncrossings);
                double mag = sqrt(pow(delX,2.0)+pow(delZ,2.0));
                vec3DW(dkx, i,j,m, nrays, ncrossings, delX);
                vec3DW(dkz, i,j,m, nrays, ncrossings, delZ);
                vec3DW(dkmag, i,j,m, nrays, ncrossings, mag);
            }
        }
    }
  
}