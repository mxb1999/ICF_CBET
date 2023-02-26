#include "CBET_Interface.hpp"
#include "parallelConfig.hpp"
#include "CBETVars.hpp"


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
  curr->dk_cu = dk;
  curr->dkmag_cu = dkmag;
  curr->uflow_cu = u_flow;
  curr->i_b_cu = i_b;
  curr->i_b_new_cu = i_b_new;
  curr->W_new_cu = W_new;
  curr->W_cu = W;
  curr->ints_cu = ints;
  curr->boxes_cu = boxes;
  curr->x_cu = x;
  curr->z_cu = z;
  curr->present_cu = present;
  return curr;
}
void initArrays()
{

    printf("SIZE %d %d %d\n", nbeams, nrays, ncrossings);
  //Allocate global memory for relevant arrays
  if(cudaCalc)
  {
    cudaMallocManaged(&i_b_new, sizeof(double)*CROSS);
    cudaError_t err = cudaMallocManaged(&i_b, sizeof(double)*CROSS);
    cudaMallocManaged(&machnum, sizeof(double)*GRID);
    cudaMallocManaged(&u_flow, sizeof(double)*GRID);
    cudaMallocManaged(&dk, sizeof(double)*CROSS*3);
    cudaMallocManaged(&dkmag, sizeof(double)*CROSS);
  }else
  {
    i_b_new = new double[CROSS]{0.0};
    i_b = new double[CROSS]{0.0};
    machnum = new double[GRID]{0.0};
    u_flow = new double[GRID]{0.0};
    //dkx = new double[CROSS]{0.0};
    //dkz = new double[CROSS]{0.0};
    //dkmag = new double[CROSS]{0.0};
  }

  double* phase_x = new double[nrays];//phase of ray
  double* pow_x = new double[nrays];//power delivery of ray
  double*  initIntensities = new double[nrays];
  span(phase_x,beam_min_z, beam_max_z, nrays);
  #pragma omp parallel for num_threads(threads)
  for(int i = 0; i < nrays;i++)
  {
    pow_x[i] = exp(-1*pow(pow(phase_x[i]/sigma,2.0),4.0/2.0));
  }
  double dbRay = (beam_max_z-beam_min_z)/nrays;
  double offsets[nrays];
  span(offsets, beam_min_z, beam_max_z, nrays);
  //#pragma omp parallel for num_threads(threads)
  int actualIntensity = intensity/1e14;
  for(int i = 0; i < nrays;i++)
  {
    //printf("%e\n",interp(phase_x, pow_x, beam_min_z+dbRay*i, nrays));
    initIntensities[i] = actualIntensity*exp(-2*pow(abs(offsets[i]/2e-4), 4.0));//interp(phase_x, pow_x, beam_min_z+dbRay*i, nrays)*intensity;
  }
  //Initialize CBET array values

  for(int m = 0; m < nbeams; m++)
  {

    //#pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nrays; j++)
    {
        double thisinit = initIntensities[j];
        printf("%e\n", thisinit);
        vec3DW(i_b, m,j,0,nrays,ncrossings,thisinit);
        vec3DW(i_b_new, m,j,0,nrays,ncrossings,thisinit);
        for(int q = 0; q < ncrossings-1; q++)
        {

            vec3DW(i_b, m,j,q+1,nrays,ncrossings,thisinit);
            //double delX = vec3D(crossesx,m,j,q+1,nrays,ncrossings)-vec3D(crossesx,m,j,q,nrays,ncrossings);
            //double delZ = vec3D(crossesz,m,j,q+1,nrays,ncrossings)-vec3D(crossesz,m,j,q,nrays,ncrossings);
            //double mag = sqrt(pow(delX, 2.0)+pow(delZ, 2.0));
           // vec3DW(dkx, m,j,q, nrays, ncrossings, delX);
           // vec3DW(dkz, m,j,q, nrays, ncrossings, delZ);
            //vec3DW(dkmag, m,j,q, nrays, ncrossings, mag);

        }
    }
  }
  double area = dx*dz;
  double* spanVal = new double[nz];
  span(spanVal,1.4,2.2,nz);
  //#pragma omp parallel for num_threads(threads)
  for(int i = 0; i < nz; i++)
  {
    double val = -1*spanVal[i];
    for(int j = 0; j < nx; j++)
    {
      double temp = val;
      vec2DW(machnum,i,j,nz, val);
      vec2DW(u_flow,i,j,nz, val*cs); //
    }
  }
  delete [] pow_x;
  delete [] phase_x;
  delete [] initIntensities;
}

