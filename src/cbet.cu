#include "CBET_Interface.hpp"
#include "parallelConfig.hpp"




//define CBET Gain Function
#define ABSLIM 1000.0
#ifndef NORM
  #define NORM 1e14
#endif
//main CBET routine
__global__
void getCBETGain(double* wMultOld, double* conv, double* logger, double medianDS, int iteration, int* raystore)
{
  double e = 4.8032e-10;
  me = 9.109384e-28;
  kb = 1.38065e-16;
  c=2.99792458e10;
  double cbetConst = (8*pi*1e7*NORM/c)*(e*e/(4*me*c*kb*1.1605e7))*1e-4;
  int i, j, k, q, p;//loop variables
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int beam = index / nrays;
  int ray = index % nrays;
  for(k = 0; k < ncrossings-1; k++)//for all crossings of ray j
  {
    int ix, iz;
    int id = vec3D(boxes, i, j, k, nrays, ncrossings);//get the spatial location
    if(!id)//no more crossings, break
    {
      break;
    }
    ix = (id - 1) % nx;//convert to x and z coordinates
    iz = (id - ix - 1)/nx;
    int islast = 0;//store whether this is the last crossing, needed for area step size calculations
    if(k != ncrossings - 1)
    {
      int idnext = vec3D(boxes, i, j, k+1, nrays, ncrossings);
      if(!idnext)//this is the last crossing of ray j
      {
        islast = 1;
      }
    }else
    {
      islast = 1;
    }
    double neOverNc = vec3D(neovernc, i, j, k, nrays, ncrossings);//
    double cbetSum = 0.0;
    double ds = vec3D(dkmag, i, j, k, nrays, ncrossings);

    for(q = 0; q < nbeams; q++)//for every other beam
    {
      if(q == i)
      {
        continue;
      }
      int ray_o = raystore[(q*nx + ix)*nz + iz] - 1;//ray determined to be the pump for this cell of beam q
      if(ray_o == -1)//no valid interactions with beam q
      {
        continue;
      }

      //RAYCROSS SUBROUTINE
      int raycross = 0;
      int q_location = vec3D(boxes, q, ray_o, raycross, nrays, ncrossings);//spatial index of the ray
      int xcurr = (q_location - 1) % nx; //x index (from MATLAB array scheme, will change based on C array structure)
      int zcurr = (q_location - xcurr - 1)/nx;//z index (from MATLAB array scheme, will change based on C array structure)
      while (xcurr != ix || zcurr != iz)//while in a different cell
      {
        raycross += 1;//increase rayCross
        if(raycross >= ncrossings)//should not occur, once all of the bugs are ironed out->can be removed to reduce branching
        {
          printf("Error: Ray %d did not pass through (%d, %d)\n", ray_o, ix, iz);
          break;
        }
        q_location = vec3D(boxes, q, ray_o, raycross, nrays, ncrossings);//new spatial index
        xcurr = (q_location - 1) % nx;//update indices
        zcurr = (q_location - xcurr - 1)/nx;
      };



      int islastq = !(vec3D(boxes, q, ray_o, raycross+1, nrays, ncrossings));//check if this is the last crossing for ray_o
      //INTERACTION MULTIPLIER: find ray_o's interaction multiplier
      double area1, area2;
      //get the average area of the ray across the zone
      area1 = vec3D(areas, q, ray_o, raycross+!islastq, nrays, ncrossings);
      area2 = vec3D(areas, q, ray_o, raycross, nrays, ncrossings);
      double areaAvg = (area1+area2)/2;
      double neOverNc = eden[ix*nz + iz]/ncrit;//NOTE: This neOverNc value can be taken straight from the grid
      double neOverNcCorrected = fmin(neOverNc, 1.0);//clamp the maximum neOverNc value to
      double neTerm = sqrt(1-neOverNcCorrected);//term used multiple times in calculations
      double epsilonEff = neTerm*neTerm;
      double interactionMult = 1/(areaAvg*neTerm)*1/sqrt(epsilonEff);


      //eta Routine
      double mag1 = ds, mag2;//get the magnitudes of the rays, can calculate here from dkx and dkz instead of storing explicitly
      mag2 = vec3D(dkmag, q, ray_o, raycross, nrays, ncrossings);
      double kx_seed = vec3D(dkx, i, j, k, nrays, ncrossings);//get the x and y components of the ray vectors
      double kz_seed = vec3D(dkz, i, j, k, nrays, ncrossings);

      double kx_pump = vec3D(dkx, q, ray_o, raycross, nrays, ncrossings);
      double kz_pump = vec3D(dkz, q, ray_o, raycross, nrays, ncrossings);

      double machx = machnum[ix*nz+iz];//the mach number of the plasma velocity
      double machz = 0.0;//TODO: Implement multidimensional plasma velocity

      double omega1 = omega, omega2 = omega; //assuming uniform wavelength/frequency
      //find ion acoustic wave vector, difference of ray trajectories scaled to omega/(sqrt(1-neOverNc)*c)

      double iawVector[] = {(omega1*kx_seed - omega2*kx_pump)*sqrt(1-neOverNc)/c,(omega1*kz_seed - omega2*kz_pump)*sqrt(1-neOverNc)/c};
      double k_iaw = sqrt(iawVector[0]*iawVector[0] + iawVector[1]*iawVector[1]);//magnitude of the iaw vector
      double etaNumerator = omega1-omega2 - (iawVector[0]*machx + iawVector[1]*machz)*cs;//numerator of eta
      double etaDenominator = k_iaw*cs;//denominator of eta
      double eta = etaNumerator/etaDenominator;//find eta
      //double eta = getEta(i, j, k, q, ray_o, raycross, ix, iz);//calculate eta from the interacting rays


      //FIND COUPLING MULTIPLIER
      double param1 = cbetConst/(omega*(Te_eV/1e3 + 3 * Ti_eV/1e3/Z)); //Split coupling multiplier into discrete chuncks->easier debugging
      double param2 = neOverNc/iaw*iaw*iaw*eta;
      double param3 = (pow(eta*eta - 1, 2) + iaw*iaw*eta*eta);
      double param4 = interactionMult;
      double couplingMult = param1*param2/param3*param4*ds;//get the coupling multiplier

      double otherIntensity1 = vec3D(i_b, q, ray_o, raycross, nrays, ncrossings);
      double otherIntensity2 = vec3D(i_b, q, ray_o, raycross + !islastq, nrays, ncrossings);
      double avgIntensity = (otherIntensity2 + otherIntensity1)/2;//average the intensity of the other ray across the cell
      cbetSum += couplingMult*avgIntensity;//add to the accumulator
    }
    double mult = exp(-1*cbetSum);//exponentiate the sum
    //LIMIT THE MULTIPLIER
    double stepRatio = ds/medianDS;//ratio of path length to median
    double multPerStep = (1-mult)/stepRatio;//ratio of multiplier to step length
    double multiplierLim = maxIncr*iteration;//maximum allowable multiplier for iteration
    multPerStep = (multPerStep > multiplierLim) ? multiplierLim : multPerStep;//abs(mult) greater than multiplier limit, clamp with appropriate sign
    multPerStep = (multPerStep < -1*multiplierLim) ? -1*multiplierLim : multPerStep;
    mult = 1 - multPerStep*stepRatio;//restore mult value
    mult = (mult > ABSLIM) ? ABSLIM : mult;//make sure mult is never above absolute maximum
    mult = (mult < 1/ABSLIM) ? 1/ABSLIM : mult;
    vec3DW(wMult, i, j, k, nrays, ncrossings, mult);//store the limited value
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
      //cbetGain<<<B2, T>>>(vars, arrays, marked,wMult, wMultOld, mi, mi_kg,maxDelta);
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

