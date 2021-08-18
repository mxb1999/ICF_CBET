#include "CBET_Interface.hpp"
#include "parallelConfig.hpp"




//define CBET Gain Function
__global__ void 
cbetGain(CBETVars* constants, CBETArrs* arrays,int* marked, double* wMult, double* wMultOld, double mi, double mi_kg, double* maxDelta,int beam = 0)
{

  int nbeams_cu = constants->nbeams_cu;
  int nrays_cu = constants->nrays_cu;
  int index = blockIdx.x*blockDim.x+threadIdx.x;
  int ncrossings_cu = constants->ncrossings_cu;
  if(index >= nbeams_cu*nrays_cu*ncrossings_cu)
  {
      return;
  }
  
  beam = index / (nrays_cu*ncrossings_cu);
  /*if(index >= (nrays_cu))
  {
      return;
  }*/
  int raynum = index/ncrossings_cu % nrays_cu;
  //imported constants
  int nx_cu = constants->nx_cu;
  int nz_cu = constants->nz_cu;
  int numstored_cu = constants->numstored_cu;
  double dx_cu = constants->dx_cu;
  double dz_cu = constants->dz_cu;
  double ncrit_cu = constants->ncrit_cu;
  double c_cu = constants->c_cu;
  double pi_cu = constants->pi_cu;
  double iaw_cu = constants->iaw_cu;
  double cs_cu = constants->cs_cu;
  double estat_cu = constants->estat_cu;
  double Ti_cu = constants->Ti_cu;
  double Te_cu = constants->Te_cu;
  double Z_cu = constants->Z_cu;
  double omega_cu = constants->omega_cu;
  double kb_cu = constants->kb_cu;
  double me_cu = constants->me_cu;

  //imported arrays
  double* i_b_cu = arrays->i_b_cu;
  double* i_b_new_cu = arrays->i_b_new_cu;
  double* W_cu = arrays->W_cu;

  double* x_cu = arrays->x_cu;

  double* z_cu = arrays->z_cu;
  double* W_new_cu = arrays->W_new_cu;
  double* dkx_cu = arrays->dkx_cu;
  double* dkz_cu = arrays->dkz_cu;
  double* dkmag_cu = arrays->dkmag_cu;
  double* uflow_cu = arrays->uflow_cu;
  
  double* eden_cu = arrays->eden_cu;
  int* boxes_cu = arrays->boxes_cu;
  int* numrays_cu = arrays->numrays_cu;
  int* present_cu = arrays->present_cu; 
  double constant1 = (pow(estat_cu,2.0))/(4*(1.0e3*me_cu)*c_cu*omega_cu*kb_cu*Te_cu*(1+3*Ti_cu/(Z_cu*Te_cu)));
  double coullim1 = Ti_cu*me_cu/mi;
  double coullim2 = 10*pow(Z_cu,2);
  int logswitch = 0;
  if(Te_cu < coullim1)
  {
    logswitch = 1;
  }else if(coullim1 < Te_cu && Te_cu < coullim2)
  {
    logswitch = 2;
  }else if(coullim1 < coullim2 && coullim2 < Te_cu)
  {
    logswitch = 3;
  }
  if(logswitch == 0)
  {
    printf("Error in inputs\n");
    return;
  }   
  int m = index % ncrossings_cu;
  //iterate over each ray beam (excepting the last one)
  //each beam will be treated as a pump beam for those preceeding, as a seed beam for those following
  double i0 = vec3D_cu(i_b_cu, beam,raynum,0, nrays_cu, ncrossings_cu);
  int max = INT32_MIN;
  
  int ix = vec4D_cu(boxes_cu, beam,raynum,m, 0, nrays_cu, ncrossings_cu, 2);
  int iz = vec4D_cu(boxes_cu, beam,raynum,m, 1, nrays_cu, ncrossings_cu, 2);
  if(!ix || !iz)
  {
    return;
  }
  ix--;
  iz--;
  int icnt = vec3D_cu(present_cu, beam, ix,iz,nx_cu,nz_cu);
  double omega1 = omega_cu;
  double mag1 = vec3D_cu(dkmag_cu, beam, raynum, m, nrays_cu, ncrossings_cu);
  double ne = vec2D_cu(eden_cu, ix,iz,nz_cu);
  double epsilon = 1.0-ne/ncrit_cu;
  double coullog;//coulomb logarithm for electron-electron collisions
  if(logswitch == 1)
  {
    double mp_kg = 1.6726219e-27;//mass of proton in kg 
    coullog = 30-log(sqrt(ne)/pow(Ti_cu, 3/2)*pow(Z_cu, 3/2)*mi_kg/mp_kg);
  }else if(logswitch == 2)
  {
    coullog = 23-log(sqrt(ne)/pow(Ti_cu, 3/2)*Z_cu);
  }else
  {
    coullog = 24-log(sqrt(ne)/Ti_cu);
  }
  double vei = 4*sqrt(2*pi_cu)*pow(estat_cu,4)*pow(Z_cu,2)*coullog/(3*sqrt(me_cu)*pow(Te_cu, 3/2));
  double L_aij = c_cu*sqrt(epsilon)*ncrit_cu/(vei*ne);
  double prevVal = vec3D_cu(wMultOld,beam,raynum,m, nrays_cu, ncrossings_cu);
  double prev = exp(-1*mag1/L_aij);
  vec3DW_cu(wMult,beam,raynum,m, nrays_cu, ncrossings_cu, exp(-1*mag1/L_aij));
  double limmult = prev;
  for(int q = 0; q < nbeams_cu;q++)
  {
    if(q == beam)
    {
      continue;
    }
    if(!vec4D_cu(marked,q,ix,iz,0, nx_cu,nz_cu,numstored_cu))
    {
      continue;
    }
    int qcnt = vec3D_cu(present_cu, q, ix,iz,nx_cu,nz_cu);
    //int qcross[qcnt]{0};
    int lim = (qcnt < icnt) ? qcnt : icnt;
    for(int l = 0; l < lim; l++)
    {
    int rayCross = 0;
    int r = vec4D_cu(marked,q,ix,iz,l, nx_cu,nz_cu,numstored_cu)-1;
    double multAcc = vec3D_cu(wMultOld, q,r,0, nrays_cu, ncrossings_cu);
    for(int p = 0; p < ncrossings_cu; p++)
    {
      int ox = vec4D_cu(boxes_cu, q,r,p, 0, nrays_cu, ncrossings_cu, 2);
      int oz = vec4D_cu(boxes_cu, q,r,p, 1, nrays_cu, ncrossings_cu, 2);
      if(!ox || !oz)
      {
        break;
      }
      ox--;
      oz--;
      multAcc *= vec3D_cu(wMultOld, q,r,p, nrays_cu, ncrossings_cu);
      if(ox == ix && oz == iz)
      {
        rayCross = p;
        break;
      }
    }

    double mag2 = vec3D_cu(dkmag_cu, q, r, rayCross, nrays_cu, ncrossings_cu);
    if(mag2 < 1e-10)
    {
      continue;
    }
    double cumProd = vec3D_cu(wMult, q, r, rayCross, nrays_cu, ncrossings_cu);
    double kmag = (omega_cu/c_cu)*sqrt(epsilon);
    double kx1 = kmag * vec3D_cu(dkx_cu, beam, raynum, m, nrays_cu, ncrossings_cu) / mag1;
    double kx2 = kmag * vec3D_cu(dkx_cu, q, r, rayCross, nrays_cu, ncrossings_cu) / mag2;
    double kz1 = kmag * vec3D_cu(dkz_cu, beam, raynum, m, nrays_cu, ncrossings_cu) / mag1;
    double kz2 = kmag * vec3D_cu(dkz_cu, q, r, rayCross, nrays_cu, ncrossings_cu) / mag2;
    double kiaw = sqrt(pow(kx2-kx1,2.0)+pow(kz2-kz1,2.0));
    double ws = kiaw*cs_cu;
    double omega2 = omega_cu;
    double eta = ((omega2-omega1)-(kx2-kx1)*vec2D_cu(uflow_cu,ix,iz,nz_cu))/(ws+1.0e-10);
    /*if(beam == 1)
    {
      printf("eta %e\n", eta);
    }*/
    double efield2 = sqrt(8.*pi_cu*1.0e7*vec3D_cu(i_b_cu, q, r, rayCross, nrays_cu, ncrossings_cu)/c_cu);   
    double P = (pow(iaw_cu,2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw_cu),2)*pow(eta,2));  
    double gain1 = constant1*pow(efield2,2)*(ne/ncrit_cu)*(1/iaw_cu)*P/icnt;               //L^-1 from Russ's paper
    double oldEnergy2 = multAcc;
    double newEnergy1Mult = exp(oldEnergy2*mag1*gain1/sqrt(epsilon));
    limmult*=newEnergy1Mult;
    //printf("multAcc %e\n", multAcc);
    }
    double curr = limmult;
    if(beam == 1  && limmult > 1.5)
    {
      //printf("Limmult %d %d %e\n",raynum, m, limmult);
    }       
    vec3DW_cu(wMult, beam,raynum,m, nrays_cu, ncrossings_cu, limmult);

    double currDev = abs(prevVal-curr)/prevVal;
    maxDelta[index] = (currDev > maxDelta[index]) ? currDev : maxDelta[index];       
    
  }
  
}

__global__ void
updateIterVals(double* wMultOld, double* wMult, double* i_b, double* i_b_new, int nbeams, int nrays, int ncrossings)
{
  int index = blockDim.x*blockIdx.x + threadIdx.x;
  if(index >= nbeams*nrays)
  {
      return;
  }
  int beam = index / (nrays);
  /*if(index >= (nrays_cu))
  {
      return;
  }*/
  int raynum = (index) % nrays;
  //double i0 = vec3D_cu(i_b, beam, raynum, 0, nrays, ncrossings);
  for(int cross = 1; cross < ncrossings;cross++)
  {
    double newMult = vec3D_cu(wMult, beam, raynum, cross, nrays, ncrossings);
    vec3DW_cu(wMultOld, beam, raynum, cross, nrays,ncrossings, newMult);
    //vec3DW_cu(i_b, beam, raynum, cross, nrays,ncrossings, newMult*i0);
    //i0 = vec3D_cu(i_b, beam, raynum, cross, nrays,ncrossings);
  }

}
void updateIterValsSerial(double* wMultOld)
{
  /*if(index >= (nrays_cu))
  {
      return;
  }*/
  for(int beam = 0; beam < nbeams; beam++)
  {
    for(int raynum = 0; raynum < nrays; raynum++)
    {
      for(int i = 0; i < ncrossings;i++)
      {
        double newMult = vec3D(wMult, beam, raynum, i, nrays, ncrossings);
        double newIntensity = vec3D(i_b_new, beam, raynum, i, nrays, ncrossings);
        double oldIntensity = vec3D(i_b, beam, raynum, i, nrays, ncrossings);
        double oldMultVal = vec3D(wMultOld, beam, raynum, i, nrays, ncrossings);
        vec3DW(wMultOld, beam, raynum, i, nrays,ncrossings, newMult);
        vec3DW(i_b, beam, raynum, i, nrays,ncrossings, newIntensity);
    
      }  
    }
  }

}

__global__ void 
cbetUpdate(int nbeams, int nrays, int ncrossings, double* wMult, double* intensity, int* boxes)
{
    
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int beam = index/nrays;
    if(beam >= (nbeams))
    {
        return;
    }
    int raynum = index % nrays;
    //imported constants
    //iterate over each ray beam (excepting the last one)
    //each beam will be treated as a pump beam for those preceeding, as a seed beam for those following
    double i0 = vec3D_cu(intensity, beam,raynum,0, nrays, ncrossings);
    for(int m = 1; m < ncrossings; m++)
    {
      double mult = vec3D_cu(wMult, beam,raynum,m, nrays, ncrossings);
      vec3DW_cu(intensity, beam,raynum,m, nrays, ncrossings, i0*mult);
      i0 = vec3D_cu(intensity, beam,raynum,m, nrays, ncrossings);
    }
}
void freeIntermediateTraceArrs()
{
  
  cudaError_t err = cudaFree(dedendx);
  if(err)
  {
    printf("%s\n", cudaGetErrorString(err));
  }
  err = cudaFree(dedendz);
  if(err)
  {
    printf("%s\n", cudaGetErrorString(err));
  }
  
  err = cudaFree(crossesx);
  if(err)
  {
    printf("%s\n", cudaGetErrorString(err));
  }
  err = cudaFree(crossesz);
  if(err)
  {
    printf("%s\n", cudaGetErrorString(err));
  }
}
void launchCBETKernel()
{
  
    if(optimize)
    {
      cbetOptimize();
    }
    freeIntermediateTraceArrs();
    CBETVars* vars = new_cbetVars();
    CBETArrs* arrays = new_cbetArrs();
    double* cbetKill; 
    double* wMultOld;   
    cudaMallocManaged(&wMult, sizeof(double)*CROSS);//store the cumulative product of the normalized ray energies

    cudaMallocManaged(&wMultOld, sizeof(double)*CROSS);//store the cumulative product of the normalized ray energies
    for(int i = 0; i < CROSS; i++)
    {
        wMultOld[i] = 1.0;

        wMult[i] = 1.0;
    }
    int T = 256;
    int B = nrays*nbeams/T+1;
    int B2 = nrays*nbeams*ncrossings/T+1;
    int cnt = 0;

    double* maxDelta;
    auto startKernel = std::chrono::high_resolution_clock::now();
    cudaMallocManaged(&maxDelta, sizeof(double)*CROSS);//store the cumulative product of the normalized ray energies
    for(int i = 0; i < maxIter; i++)
    {
      //printf("Test1\n");
      cbetGain<<<B2, T>>>(vars, arrays, marked,wMult, wMultOld, mi, mi_kg,maxDelta);
      cudaDeviceSynchronize();
      //printf("Test2\n");
      updateIterVals<<<B, T>>>(wMultOld, wMult, i_b, i_b_new, nbeams, nrays, ncrossings);
      cudaDeviceSynchronize();
      //updateIterValsSerial(wMultOld);
      double max = 0.0;
      #pragma omp parallel for num_threads(threads)
      for(int j = 0; j < CROSS; j++)
      {
        max = fmax(maxDelta[j], max);
        maxDelta[j] = 0.0;
      }
      printf("Max %e\n", max);
      if(max <= converge)
      {
        break;
      }

    }
    cbetUpdate<<<B, T>>>(nbeams, nrays, ncrossings, wMult, i_b_new,boxes);
    cudaDeviceSynchronize();
    auto stopKernel = std::chrono::high_resolution_clock::now();    
    cudaFree(wMult);
    cudaFree(wMultOld);
    cudaFree(maxDelta);
    //*output << "CUCBET "<< nrays << " " << chrono::duration_cast<chrono::milliseconds>(stopKernel-startKernel).count() << std::endl;
  
}

