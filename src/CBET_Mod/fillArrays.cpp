#include "CBET_Interface.hpp"
#include "parallelConfig.hpp"


CBETVars* new_cbetVars()
{
  CBETVars* curr;
  cudaMallocManaged(&curr, sizeof(CBETVars));
  if(curr == NULL)
  {
    return NULL;
  }
  curr->estat_cu = estat;
  curr->me_cu = me;
  curr->c_cu = c;
  curr->kb_cu = kb;
  curr->Te_cu = Te;
  curr->Ti_cu = Ti;
  curr->nbeams_cu = nbeams;
  curr->nx_cu = nx;
  curr->numstored_cu = numstored;
  curr->nz_cu = nz;
  curr->dx_cu = dx;
  curr->dz_cu = dz;
  curr->nrays_cu = nrays;
  curr->ncrossings_cu = ncrossings;
  curr->ncrit_cu = ncrit;
  curr->pi_cu = pi;
  curr->iaw_cu = iaw;
  curr->cs_cu = cs;
  curr->estat_cu = estat;
  curr->Z_cu = Z;
  curr->omega_cu = omega;
  return curr;
}

CBETArrs* new_cbetArrs()
{
  CBETArrs* curr;
  cudaMallocManaged(&curr, sizeof(CBETArrs));
  if(curr == NULL)
  {
    return NULL;
  }
  curr->eden_cu = eden;
  curr->dkx_cu = dkx;  
  curr->dkz_cu = dkz;
  curr->dkmag_cu = dkmag;
  curr->uflow_cu = u_flow;
  curr->i_b_cu = i_b;
  curr->i_b_new_cu = i_b_new;
  curr->W_new_cu = W_new;
  curr->W_cu = W;
  curr->numrays_cu = numrays;
  curr->ints_cu = ints;
  curr->boxes_cu = boxes;
  curr->x_cu = x;
  curr->z_cu = z;
  curr->present_cu = present;
  return curr;
}
void initArrays()
{
    //Allocate global memory for relevant arrays
    cudaMallocManaged(&i_b_prev, sizeof(double)*CROSS*RAYS);
    cudaMallocManaged(&i_b_new, sizeof(double)*CROSS);
    cudaMallocManaged(&W_new, sizeof(double)*CROSS);
    cudaMallocManaged(&W, sizeof(double)*CROSS);
    cudaMallocManaged(&i_b, sizeof(double)*CROSS);
      //Initialize CBET array values
  double area = dx*dz;
  for(int m = 0; m < nbeams; m++)
  {
    //#pragma omp parallel for num_threads(threads)
    for(int i = 0; i < nrays;i++)
    {
      for(int j = 0; j < ncrossings;j++)
      {
        int boxx = vec4D(boxes, m,i,j,0,nrays,ncrossings,2);
        int boxz = vec4D(boxes, m,i,j,1,nrays,ncrossings,2);
        if(!boxx || !boxz)
        {
          continue;
        }
        boxx--;
        boxz--;
        double energyDep = vec3D(edep_flat,m,boxx,boxz,nx+2,nz+2)/vec3D(present,m,boxx,boxz, nx,nz);//break up pipeline to prevent RAW dependencies
        double initEnergy = sqrt(1.0-vec2D(eden,boxx,boxz,nz)/ncrit);
        vec3DW(i_b_new,m,i,j,nrays, ncrossings, energyDep);///present[m][i][j];
        vec3DW(W,m,i,j,nrays,ncrossings, initEnergy);///present[m][i][j];
        vec3DW(i_b,m,i,j,nrays, ncrossings, energyDep);///present[m][i][j];
        vec3DW(W_new,m,i,j,nrays,ncrossings, initEnergy);///present[m][i][j];
        double flownum = vec2D(machnum,boxx,boxz,nz)*cs;
        vec2DW(u_flow,boxx,boxz,nz, flownum); //
      }

    }
    #pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nrays; j++)
    {
        for(int q = 0; q < ncrossings-1; q++)
        {
            double delX = vec3D(crossesx,m,j,q+1,nrays,ncrossings)-vec3D(crossesx,m,j,q,nrays,ncrossings);
            double delZ = vec3D(crossesz,m,j,q+1,nrays,ncrossings)-vec3D(crossesz,m,j,q,nrays,ncrossings);
            double mag = sqrt(pow(delX,2.0)+pow(delZ,2.0));
            vec3DW(dkx, m,j,q, nrays, ncrossings, delX);
            vec3DW(dkz, m,j,q, nrays, ncrossings, delZ);
            vec3DW(dkmag, m,j,q, nrays, ncrossings, mag);
        }
    }
  }
  
}