#include "CBET_Interface.hpp"
#define ABSLIM 1000
//limit the energy change allowed in a given cell in a given iteration
#define NORM 1e14
//main CBET routine
void getCBETGain(double* wMult, double* i_b, double* eden,  double* areas, double* machnum, int* boxes, double* kvec, double medianDS, int* marked, int iteration, double* convergence)
{

    double cbetconst = (8*pi*1e7*NORM/c)*(estat*estat/(4*me*c*kb*1.1605e7))*1e-4;
    //#pragma omp parallel for num_threads(omp_get_max_threads())
    for(int i = 0; i < nbeams; i++)//for all beams in system
    {
        for(int j = 0; j < nrays; j++)//for all rays in beam i
        {
          if(fabs(vec3D(i_b, i, j, 0, nrays, ncrossings)) < 1e-12)
          {
            continue;
          }
          for(int k = 0; k < ncrossings; k++)
          {
            int ix, iy, iz;
            ix = vec4D(boxes, i, j, k, 0, nrays, ncrossings, 3) - 1;//get the spatial location
                    iy = vec4D(boxes, i, j, k, 1, nrays, ncrossings, 3) - 1;//get the spatial location
                    iz = vec4D(boxes, i, j, k, 2, nrays, ncrossings, 3) - 1;//get the spatial location
                    if(ix == -1)//no more crossings, break
                    {
                        break;
                    }
                    double cbetSum = 0.0;
                    double kxi, kyi, kzi, ds;
                    kxi = VEC4D(kvec, i, j, k, 0, nrays, ncrossings, 3);        //get the ray vector across the cell
                    kyi = VEC4D(kvec, i, j, k, 1, nrays, ncrossings, 3);
                    kzi = VEC4D(kvec, i, j, k, 2, nrays, ncrossings, 3);
                    ds = sqrt(kxi*kxi + kyi*kyi + kzi*kzi);
                    kxi/=ds;    //normalize to step size
                    kyi/=ds;
                    kzi/=ds;
                    for(int q = 0; q < nbeams; q++)//for every other i
                    {
                        if(q == i)
                        {
                          continue;
                        }
                        int ray_o = VEC4D(marked, q, ix, iy, iz, nx, ny, nz) - 1;//ray determined to be the pump for this cell of i q
                        if(ray_o == -1)//no valid interactions with i q
                        {
                          continue;
                        }
                        //RAYCROSS SUBROUTINE
                        int raycross = 0, xcurr, ycurr, zcurr;
                        xcurr = VEC4D(boxes, q, ray_o, raycross, 0, nrays, ncrossings, 3) - 1;//new spatial index
                        ycurr = VEC4D(boxes, q, ray_o, raycross, 1, nrays, ncrossings, 3) - 1;
                        zcurr = VEC4D(boxes, q, ray_o, raycross, 2, nrays, ncrossings, 3) - 1;
                        while (xcurr != ix || zcurr != iz || ycurr != iy)//while in a different cell
                        {
                            raycross += 1;//increase rayCross
                            if(raycross >= ncrossings || !(VEC4D(boxes, q, ray_o, raycross, 0, nrays, ncrossings, 3)))//should not occur, once all of the bugs are ironed out->can be removed to reduce branching
                            {
                                printf("Error: Ray %d of i %d did not pass through (%d, %d, %d)\n", ray_o,q, ix, iy, iz);
                                break;
                            }
                            xcurr = VEC4D(boxes, q, ray_o, raycross, 0, nrays, ncrossings, 3) - 1;//new spatial index
                            ycurr = VEC4D(boxes, q, ray_o, raycross, 1, nrays, ncrossings, 3) - 1;
                            zcurr = VEC4D(boxes, q, ray_o, raycross, 2, nrays, ncrossings, 3) - 1;
                        };


                        int islastq = !(VEC4D(boxes, q, ray_o, raycross+1, 0, nrays, ncrossings, 3));//check if this is the last crossing for ray_o
                        //INTERACTION MULTIPLIER: find ray_o's interaction multiplier
                        double area1, area2;
                        //get the average area of the ray across the zone
                        area1 = VEC3D(areas, q, ray_o, raycross+!islastq, nrays, ncrossings);
                        area2 = VEC3D(areas, q, ray_o, raycross, nrays, ncrossings);
                        double areaAvg = (area1+area2)/2;
                        double neOverNc = eden[ix*nz + iz]/ncrit;//NOTE: This neOverNc value can be taken straight from the grid
                        double neOverNcCorrected = fmin(neOverNc, 1.0);//clamp the maximum neOverNc value to
                        double neTerm = sqrt(1-neOverNcCorrected);//term used multiple times in calculations
                        double epsilonEff = neTerm*neTerm;
                        double interactionMult = 1/(areaAvg*neTerm)*1/sqrt(epsilonEff);
                        //eta Routine
                        double mag2, kx_pump, ky_pump, kz_pump;//get the magnitudes of the rays, can calculate here from dkx and dkz instead of storing explicitly
                        kx_pump = VEC4D(kvec, q, ray_o, raycross, 0, nrays, ncrossings, 3);        //get the ray vector across the cell
                        ky_pump = VEC4D(kvec, q, ray_o, raycross, 1, nrays, ncrossings, 3);
                        kz_pump = VEC4D(kvec, q, ray_o, raycross, 2, nrays, ncrossings, 3);
                        mag2 = sqrt(kx_pump*kx_pump + ky_pump*ky_pump + kz_pump*kz_pump);
                        kx_pump /= mag2;
                        ky_pump /= mag2;
                        kz_pump /= mag2;
                        double machx, machy, machz;

                        machx = machnum[ix*nz+iz];//VEC4D(machnum, ix, iy, iz, 0, ny, nz, 3);//the mach number of the plasma velocity
                        machy = 0.0;//VEC4D(machnum, ix, iy, iz, 1, ny, nz, 3);
                        machz = 0.0;//VEC4D(machnum, ix, iy, iz, 2, ny, nz, 3);
                        double omega1 = omega, omega2 = omega;
                        double omega1_iaw = omega*sqrt(1-neOverNc)/c;
                        double omega2_iaw = omega1_iaw; //assuming uniform wavelength/frequency
                        //find ion acoustic wave vector, difference of ray trajectories scaled to omega/(sqrt(1-neOverNc)*c)

                        double iawVector[] = {(omega1_iaw*kxi - omega2_iaw*kx_pump),\
                                              (omega1_iaw*kyi - omega2_iaw*ky_pump),\
                                              (omega1_iaw*kzi - omega2_iaw*kz_pump)};
                        double k_iaw = sqrt(iawVector[0]*iawVector[0] + iawVector[1]*iawVector[1] + iawVector[2]*iawVector[2]);//magnitude of the iaw vector
                        double etaNumerator = omega1 - omega2 - (iawVector[0]*machx + iawVector[1]*machy + iawVector[2]*machz)*cs;//numerator of eta
                        double etaDenominator = k_iaw*cs;//denominator of eta
                        double eta = etaNumerator/etaDenominator;//find eta
                        //double eta = getEta(i, j, k, q, ray_o, raycross, ix, iz);//calculate eta from the interacting rays


                        //FIND COUPLING MULTIPLIER
                        double param1 = cbetconst/(omega*(Te_eV/1e3 + 3 * Ti_eV/1e3/Z)); //Split coupling multiplier into discrete chuncks->easier debugging
                        double param2 = neOverNc/iaw*iaw*iaw*eta;
                        double param3 = (pow(eta*eta - 1, 2) + iaw*iaw*eta*eta);
                        double param4 = interactionMult;
                        double couplingMult = param1*param2/param3*param4*ds;//get the coupling multiplier

                        double otherIntensity1 = VEC3D(i_b, q, ray_o, raycross, nrays, ncrossings);
                        double otherIntensity2 = VEC3D(i_b, q, ray_o, raycross + !islastq, nrays, ncrossings);
                        double avgIntensity = (otherIntensity2 + otherIntensity1)/2;//average the intensity of the other ray across the cell
                        cbetSum += couplingMult*avgIntensity;//add to the accumulator
                    }
                    double mult = exp(-1*cbetSum);//exponentiate the sum
                    //LIMIT THE MULTIPLIER
                    double prev = VEC3D(wMult, i, j, k, nrays, ncrossings);
                    double denom = fabs(prev);
                    denom += (denom < 1e-10);
                    double numer = fabs(mult - prev);
                    double stepRatio = ds/medianDS;//ratio of path length to median
                    double multPerStep = (1-mult)/stepRatio;//ratio of multiplier to step length
                    double multiplierLim = 2e-4;//maximum allowable multiplier for iteration
                    multPerStep = (multPerStep > multiplierLim) ? multiplierLim : multPerStep;//abs(mult) greater than multiplier limit, clamp with appropriate sign
                    multPerStep = (multPerStep < -1*multiplierLim) ? -1*multiplierLim : multPerStep;
                    mult = 1 - multPerStep*stepRatio;//restore mult value
                    mult = (mult > ABSLIM) ? ABSLIM : mult;//make sure mult is never above absolute maximum
                    mult = (mult < 1/ABSLIM) ? 1/ABSLIM : mult;
                    VEC3D(wMult, i, j, k, nrays, ncrossings) = mult;//store the limited value

                  }
        }
    }
}
void updateIntensities(double* i_b, double* wMult, int iteration, double* convMax, double currMax)
{
  int tr = omp_get_max_threads();
  double convMaxTR[tr] = {0};
  int convMaxTRInd[tr] = {0};
  int convMaxTRCross[tr] = {0};
  #pragma omp parallel for num_threads(tr)
  for(int i = 0; i < nbeams; i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      double i0 = VEC3D(i_b, i, j, 0, nrays, ncrossings);
      if(fabs(i0) < 1e-12)
      {
        continue;
      }
      double multAcc = 1.0;



  for(int k = 1; k < ncrossings; k++)
  {
    double mult = VEC3D(wMult, i, j, k, nrays, ncrossings);
    multAcc *= mult;

    //VEC3D(spatialLog, i, j, k, nrays, ncrossings, new_intensity/1e14);
    double i_prev = VEC3D(i_b, i, j, k, nrays, ncrossings);//previous value in cell

    double i_curr = i0*multAcc;//"current" value, unclamped update
    //printf("%e %e\n", i_curr, i_prev);
    double fractionalChange = fabs(i_curr-i_prev)/i_prev;//the fractional change in energy from imposing the update as is
    //printf("(%e) ", i0);
    if(fractionalChange > 1e100) printf("%e\n", multAcc);
    //convMaxTR[omp_get_thread_num()] = fmax(fractionalChange, convMaxTR[omp_get_thread_num()]);//update the convergence check variable
    double curr = *convMax;
    if(convMaxTR[omp_get_thread_num()] < fractionalChange)
    {
      convMaxTR[omp_get_thread_num()] = fractionalChange;
    }
    /*unsigned long long currull = __double_as_longlong(curr);
    if(curr < fractionalChange)
    {
      unsigned long long varull;
      do{
        varull = currull;
        currull = atomicCAS((unsigned long long int*)convMax, currull, __double_as_longlong(fractionalChange));
      }while(currull != varull);
    }*/
    double correction = 0.0;
    int cond = (fractionalChange > currMax);
    int sign = (i_curr - i_prev > 0) ? 1 : -1;
    correction = 1 + (currMax*sign)*cond;
    i_curr = i_prev*correction*cond + i_curr*(!cond);
    VEC3D(i_b, i, j, k, nrays, ncrossings) = i_curr;
  }
    }
  }
  for(int i = 0; i < tr; i++)
  {
    *convMax = fmax(*convMax, convMaxTR[i]);
  }
}
#define CPUCBET

#define CBETCONVERGENCE 0.998
void cbet(double* wMult, double* i_b, double* eden, double* areas, double* machnum, int* boxes, double* kvec, double medianDS, int* marked, double tolerance)
{
  double* wMult_old;
  cudaMalloc(&wMult_old, sizeof(double)*nrays*ncrossings*nbeams);
  cudaMemcpy(wMult_old, wMult, sizeof(double)*nrays*ncrossings*nbeams, cudaMemcpyDeviceToDevice);
  double* convergence;
  cudaMallocManaged(&convergence, sizeof(double));
  int iteration = 0;
  int maxiter = 3000;
  double maxIncr = 0.02;
  double currmax = maxIncr;
  int threads_per_block = 256;
  int blocks = nbeams*nrays/threads_per_block + 1;
  printf("%d\n", blocks);
  cudaDeviceSynchronize();
  int sense_rev = 0;
    int temp = nz;
  nz = ny;
  ny = temp;
  printf("%e\n", intensity);
  do
  {
    iteration++;
    sense_rev++;
    *convergence = 0.0;
    #ifdef GPUCBET
    cbet_gain<<<blocks, threads_per_block>>>(wMult, i_b, eden, areas, machnum, boxes, kvec, medianDS, marked, iteration, convergence);
    update_i_cu<<<blocks, threads_per_block>>>(i_b, wMult, iteration, convergence, currmax);
    cudaDeviceSynchronize();
    #endif
    #ifdef CPUCBET
    getCBETGain(wMult, i_b, eden, areas, machnum, boxes, kvec, medianDS, marked, iteration, convergence);
    updateIntensities(i_b, wMult, iteration, convergence, currmax);
    #endif
    if(sense_rev % 1000 == 0) {
      sense_rev = 0;
    }
    double currmaxa = maxIncr*pow(CBETCONVERGENCE, iteration);
    double currmaxb = CBETCONVERGENCE*(*convergence);
    //printf("%e\n", *convergence);
    
    cudaMemcpy(wMult_old, wMult, sizeof(double)*nrays*ncrossings*nbeams, cudaMemcpyDeviceToDevice);
    printf("%d %e\n", iteration, *convergence);
    //currmax = fmin(currmaxa, currmaxb);
  }while(*convergence > tolerance);
  cudaDeviceSynchronize();
}
