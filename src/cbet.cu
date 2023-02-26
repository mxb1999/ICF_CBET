#include "CBET_Interface.hpp"
#include "parallelConfig.hpp"

#define ABSLIM 1000
//limit the energy change allowed in a given cell in a given iteration
#define NORM 1e14
CBETMultiplierFields* new_cbet_args() {
  CBETMultiplierFields* args = new CBETMultiplierFields;
  if(!args) {
    return NULL;
  }
  cudaMalloc(&args->boxes, sizeof(int)*CROSS);
  cudaMalloc(&args->marked, sizeof(int)*CROSS);
  cudaMalloc(&args->wMult, sizeof(int)*CROSS);
  cudaMalloc(&args->i_b, sizeof(int)*CROSS);
  cudaMalloc(&args->eden, sizeof(int)*CROSS);
  cudaMalloc(&args->areas, sizeof(int)*CROSS);
  cudaMalloc(&args->machnum, sizeof(int)*CROSS);
  cudaMalloc(&args->boxes, sizeof(int)*CROSS);
  cudaMalloc(&args->boxes, sizeof(int)*CROSS);
  cudaMalloc(&args->boxes, sizeof(int)*CROSS);
  cudaMalloc(&args->boxes, sizeof(int)*CROSS);
  cudaMalloc(&args->boxes, sizeof(int)*CROSS);
};


double cbet_conv(CBETMultiplierFields* target) {
 return 0.0;

};
void copy_to(CBETMultiplierFields* target) {


};
void copy_back(CBETMultiplierFields* target) {


};
//now we make kernel
__global__
void cbet_gain(double* wMult, double* i_b, double* eden, double* areas, double* machnum, int* boxes, double* kvec, double medianDS, int* marked, int iteration, double* convergence)
{
  int index = (blockIdx.x*blockDim.x + threadIdx.x);
  if(index >= nbeams*nrays_cu)
    return;
  int raynum = index % nrays_cu;
  int beam = index / nrays_cu;
  if(fabs(VEC3D(i_b, beam, raynum, 0, nrays_cu, ncrossings)) < 1e-12)
  {
    return;
  }
  double cbetconst = (8*pi_cu*1e7*NORM/c_cu)*(e*e/(4*me_cu*c_cu*kb_cu*1.1605e7))*1e-4;
  for(int k = 0; k < ncrossings; k++)
  {
    int ix, iy, iz;
    ix = VEC4D(boxes, beam, raynum, k, 0, nrays_cu, ncrossings, 3) - 1;//get the spatial location
            iy = VEC4D(boxes, beam, raynum, k, 1, nrays_cu, ncrossings, 3) - 1;//get the spatial location
            iz = VEC4D(boxes, beam, raynum, k, 2, nrays_cu, ncrossings, 3) - 1;//get the spatial location
            if(ix == -1)//no more crossings, break
            {
                break;
            }
            double cbetSum = 0.0;
            double kxi, kyi, kzi, ds;
            kxi = VEC4D(kvec, beam, raynum, k, 0, nrays_cu, ncrossings, 3);        //get the ray vector across the cell
            kyi = VEC4D(kvec, beam, raynum, k, 1, nrays_cu, ncrossings, 3);
            kzi = VEC4D(kvec, beam, raynum, k, 2, nrays_cu, ncrossings, 3);
            ds = sqrt(kxi*kxi + kyi*kyi + kzi*kzi);
            kxi/=ds;    //normalize to step size
            kyi/=ds;
            kzi/=ds;
            for(int q = 0; q < nbeams; q++)//for every other beam
            {
                if(q == beam)
                {
                  continue;
                }
                int ray_o = VEC4D(marked, q, ix, iy, iz, nx, ny, nz) - 1;//ray determined to be the pump for this cell of beam q
                if(ray_o == -1)//no valid interactions with beam q
                {
                continue;
                }

                //RAYCROSS SUBROUTINE
                int raycross = 0, xcurr, ycurr, zcurr;
                xcurr = VEC4D(boxes, q, ray_o, raycross, 0, nrays_cu, ncrossings, 3) - 1;//new spatial index
                ycurr = VEC4D(boxes, q, ray_o, raycross, 1, nrays_cu, ncrossings, 3) - 1;
                zcurr = VEC4D(boxes, q, ray_o, raycross, 2, nrays_cu, ncrossings, 3) - 1;
                while (xcurr != ix || zcurr != iz || ycurr != iy)//while in a different cell
                {
                    raycross += 1;//increase rayCross
                    if(raycross >= ncrossings || !(VEC4D(boxes, q, ray_o, raycross, 0, nrays, ncrossings, 3)))//should not occur, once all of the bugs are ironed out->can be removed to reduce branching
                    {
                        printf("Error: Ray %d of beam %d did not pass through (%d, %d, %d)\n", ray_o,q, ix, iy, iz);
                        break;
                    }
                    xcurr = VEC4D(boxes, q, ray_o, raycross, 0, nrays_cu, ncrossings, 3) - 1;//new spatial index
                    ycurr = VEC4D(boxes, q, ray_o, raycross, 1, nrays_cu, ncrossings, 3) - 1;
                    zcurr = VEC4D(boxes, q, ray_o, raycross, 2, nrays_cu, ncrossings, 3) - 1;
                };

                int islastq = !(VEC4D(boxes, q, ray_o, raycross+1, 0, nrays_cu, ncrossings, 3));//check if this is the last crossing for ray_o
                //INTERACTION MULTIPLIER: find ray_o's interaction multiplier
                double area1, area2;
                //get the average area of the ray across the zone
                area1 = VEC3D(areas, q, ray_o, raycross+!islastq, nrays_cu, ncrossings);
                area2 = VEC3D(areas, q, ray_o, raycross, nrays_cu, ncrossings);
                double areaAvg = (area1+area2)/2;
                double neOverNc = VEC3D(eden, ix, iy, iz, ny, nz)/ncrit;//NOTE: This neOverNc value can be taken straight from the grid
                //printf("%e %e\n", VEC3D(eden, ix, iy, iz, ny, nz), ncrit);
                double neOverNcCorrected = fmin(neOverNc, 1.0);//clamp the maximum neOverNc value to
                double neTerm = sqrt(1-neOverNcCorrected);//term used multiple times in calculations
                double epsilonEff = neTerm*neTerm;
                double interactionMult = 1/(areaAvg*neTerm)*1/sqrt(epsilonEff);
                //eta Routine
                double mag2, kx_pump, ky_pump, kz_pump;//get the magnitudes of the rays, can calculate here from dkx and dkz instead of storing explicitly
                kx_pump = VEC4D(kvec, q, ray_o, raycross, 0, nrays_cu, ncrossings, 3);        //get the ray vector across the cell
                ky_pump = VEC4D(kvec, q, ray_o, raycross, 1, nrays_cu, ncrossings, 3);
                kz_pump = VEC4D(kvec, q, ray_o, raycross, 2, nrays_cu, ncrossings, 3);
                mag2 = sqrt(kx_pump*kx_pump + ky_pump*ky_pump + kz_pump*kz_pump);
                kx_pump /= mag2;
                ky_pump /= mag2;
                kz_pump /= mag2;
                double machx, machy, machz;

                machx = VEC4D(machnum, ix, iy, iz, 0, ny, nz, 3);//the mach number of the plasma velocity
                machy = VEC4D(machnum, ix, iy, iz, 1, ny, nz, 3);
                machz = VEC4D(machnum, ix, iy, iz, 2, ny, nz, 3);
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
                  //double eta = getEta(beam, raynum, k, q, ray_o, raycross, ix, iz);//calculate eta from the interacting rays

                  //FIND COUPLING MULTIPLIER
                  double param1 = cbetconst/(omega*(Te_eV/1e3 + 3 * Ti_eV/1e3/Z)); //Split coupling multiplier into discrete chuncks->easier debugging
                  double param2 = neOverNc*iaw*eta;
                  double param3 = (pow(eta*eta - 1, 2) + iaw*iaw*eta*eta);
                  double param4 = interactionMult;
                  double couplingMult = param1*param2/param3*param4*ds;//get the coupling multiplier

                  double otherIntensity1 = VEC3D(i_b, q, ray_o, raycross, nrays_cu, ncrossings);
                  double otherIntensity2 = VEC3D(i_b, q, ray_o, raycross + !islastq, nrays_cu, ncrossings);
                  double avgIntensity = (otherIntensity2 + otherIntensity1)/(2);//average the intensity of the other ray across the cell
                  cbetSum += couplingMult*avgIntensity;//add to the accumulator
                  if(raynum == 659 && beam == 0 && k == 66)
                    printf("Coupling mult = %e: (%e %e %e %e)\nAvgIntensity = %e\n\n", couplingMult, param1, param2, param3, param4, avgIntensity);
                }
                double mult = exp(-1*cbetSum);//exponentiate the sum
                double tmp = mult;
                double prev = VEC3D(wMult, beam, raynum, k, nrays_cu, ncrossings);
                double denom = (fabs(prev) < 1e-12) ? 1.0 : prev;
                //need to check for convergence in the multipliers
                double fractionalChange = fabs(mult-prev)/denom;//the fractional change in energy from imposing the update as is
                double curr = *convergence;
                unsigned long long currull = __double_as_longlong(curr);
                if(curr < fractionalChange)
                {
                  unsigned long long varull;
                  do{
                    varull = currull;
                    currull = atomicCAS((unsigned long long int*)convergence, currull, __double_as_longlong(fractionalChange));
                  }while(currull != varull && fractionalChange > __longlong_as_double(currull));
                }
            double stepRatio = ds/medianDS;//ratio of path length to median
            double multPerStep = (1-mult)/stepRatio;//ratio of multiplier to step length
            double multiplierLim = 0.02;//maximum allowable multiplier for iteration
            multPerStep = (multPerStep > multiplierLim) ? multiplierLim : multPerStep;//abs(mult) greater than multiplier limit, clamp with appropriate sign
            multPerStep = (multPerStep < -1*multiplierLim) ? -1*multiplierLim : multPerStep;
            mult = 1 - multPerStep*stepRatio;//restore mult value
            mult = (mult > ABSLIM) ? ABSLIM : mult;//make sure mult is never above absolute maximum
            mult = (mult < 1/ABSLIM) ? 1/ABSLIM : mult;

            VEC3D(wMult, beam, raynum, k, nrays_cu, ncrossings) = mult;//store the limited value
      }


  }

  __global__
  void update_i_cu(double* i_b, double* wMult, int iteration, double* convMax, double currMax)
  {
  int index = (blockIdx.x*blockDim.x + threadIdx.x);
  if(index >= nbeams*nrays_cu)
    return;
  int raynum = index % nrays_cu;
  int beam = index / nrays_cu;
  double i0 = VEC3D(i_b, beam, raynum, 0, nrays_cu, ncrossings);
  double multAcc = 1.0;
  if(fabs(i0) < 1e-12)
  {
    return;
  }

  for(int k = 1; k < ncrossings; k++)
  {
    double mult = VEC3D(wMult, beam, raynum, k, nrays_cu, ncrossings);
    multAcc *= mult;
    //VEC3D(spatialLog, beam, raynum, k, nrays, ncrossings, new_intensity/1e14);
    double i_prev = VEC3D(i_b, beam, raynum, k, nrays_cu, ncrossings);//previous value in cell
    double i_curr = VEC3D(i_b, beam, raynum, k-1, nrays_cu, ncrossings)*mult;//"current" value, unclamped update
    //printf("%e %e\n", i_curr, i_prev);
    double fractionalChange = fabs(i_curr-i_prev)/i_prev;//the fractional change in energy from imposing the update as is
    //printf("(%e) ", i0);
    //convMaxTR[omp_get_thread_num()] = fmax(fractionalChange, convMaxTR[omp_get_thread_num()]);//update the convergence check variable
    double curr = *convMax;
    unsigned long long currull = __double_as_longlong(curr);
    if(curr < fractionalChange)
    {
      unsigned long long varull;
      do{
        varull = currull;
        currull = atomicCAS((unsigned long long int*)convMax, currull, __double_as_longlong(fractionalChange));
      }while(currull != varull && fractionalChange > __longlong_as_double(currull));
    }
    double correction = 0.0;
    int cond = (fractionalChange > currMax);
    int sign = (i_curr - i_prev > 0) ? 1 : -1;
    correction = 1 + (currMax*sign)*cond;
    i_curr = i_prev*correction*cond + i_curr*(!cond);
    VEC3D(i_b, beam, raynum, k, nrays_cu, ncrossings) = i_curr;
  }
}
//main CBET routine
void getCBETGain(double* wMult, double* i_b, double* eden,  double* areas, double* machnum, int* boxes, double* kvec, double medianDS, int* marked, int iteration)
{
    double cbetconst = (8*pi*1e7*NORM/c)*(e*e/(4*me*c*kb*1.1605e7))*1e-4;
    #pragma omp parallel for num_threads(omp_get_max_threads())
    for(int i = 0; i < nbeams; i++)//for all beams in system
    {
        for(int j = 0; j < nrays_cu; j++)//for all rays in beam i
        {
          if(fabs(VEC3D(i_b, i, j, 0, nrays_cu, ncrossings)) < 1e-12)
          {
            continue;
          }
          for(int k = 0; k < ncrossings; k++)
          {
            int ix, iy, iz;
            ix = VEC4D(boxes, i, j, k, 0, nrays_cu, ncrossings, 3) - 1;//get the spatial location
                    iy = VEC4D(boxes, i, j, k, 1, nrays, ncrossings, 3) - 1;//get the spatial location
                    iz = VEC4D(boxes, i, j, k, 2, nrays, ncrossings, 3) - 1;//get the spatial location
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

                        machx = VEC4D(machnum, ix, iy, iz, 0, ny, nz, 3);//the mach number of the plasma velocity
                        machy = VEC4D(machnum, ix, iy, iz, 1, ny, nz, 3);
                        machz = VEC4D(machnum, ix, iy, iz, 2, ny, nz, 3);
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
                        //printf("%e %e\n", couplingMult, avgIntensity);
                    }
                    double mult = exp(-1*cbetSum);//exponentiate the sum
                    //LIMIT THE MULTIPLIER
                    double stepRatio = ds/medianDS;//ratio of path length to median
                    double multPerStep = (1-mult)/stepRatio;//ratio of multiplier to step length
                    double multiplierLim = 0.2;//maximum allowable multiplier for iteration
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
/*//main CBET routine
void getCBETGain(double* wMultOld, double* conv, double* logger, double medianDS, int iteration, int* raystore)
{
  double cbetConst = getConstant();//main CBET constant
  int i, j, k, q, p;//loop variables

  for(i = 0; i < nbeams; i++)//for all beams in system
  {
    for(j = 0; j < nrays; j++)//for all rays in beam i
    {
      for(k = 0; k < ncrossings-1; k++)//for all crossings of ray j
      {
        int ix, iz;
        int id = VEC3D(boxes, i, j, k, nrays, ncrossings);//get the spatial location
        if(!id)//no more crossings, break
        {
          break;
        }
        ix = (id - 1) % nx;//convert to x and z coordinates
        iz = (id - ix - 1)/nx;
        int islast = 0;//store whether this is the last crossing, needed for area step size calculations
        if(k != ncrossings - 1)
        {
          int idnext = VEC3D(boxes, i, j, k+1, nrays, ncrossings);
          if(!idnext)//this is the last crossing of ray j
          {
            islast = 1;
          }
        }else
        {
          islast = 1;
        }
        double neOverNc = VEC3D(neovernc, i, j, k, nrays, ncrossings);//
        double cbetSum = 0.0;
        double ds = VEC3D(dkmag, i, j, k, nrays, ncrossings);

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
          int raycross = getCrossing(q, ray_o, ix, iz);//get the crossing index of the other ray
          int islastq = !(VEC3D(boxes, q, ray_o, raycross+1, nrays, ncrossings));//check if this is the last crossing for ray_o
          double interactionMultiplier = getInteractionMultiplier(q, ray_o, raycross, ix, iz, islast);//find ray_o's interaction multiplier
          double eta = getEta(i, j, k, q, ray_o, raycross, ix, iz);//calculate eta from the interacting rays
          double couplingMult = getCouplingMultiplier(i, j, k, ix, iz, cbetConst, eta, interactionMultiplier);//get the coupling multiplier
          couplingMult *= VEC3D(dkmag, i,  j, k, nrays, ncrossings);//scale the coupling multiplier by the magnitude of the seed ray's path length
          double otherIntensity1 = VEC3D(i_b, q, ray_o, raycross, nrays, ncrossings);
          double otherIntensity2 = VEC3D(i_b, q, ray_o, raycross + !islastq, nrays, ncrossings);
          double avgIntensity = (otherIntensity2 + otherIntensity1)/2;//average the intensity of the other ray across the cell
          cbetSum += couplingMult*avgIntensity;//add to the accumulator
          VEC3D(spatialLog, i,  j, k, nrays, ncrossings, couplingMult);//logging variable
        }
        double mult = exp(-1*cbetSum);//exponentiate the sum
        mult = limitMultiplier(mult, i, j, k, medianDS, iteration);//limit the multiplier
        VEC3D(wMult, i, j, k, nrays, ncrossings, mult);//store the limited value
      }
    }
  }
  getchar();
}
*/

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
#define GPUCBET
//#define CPUCBET
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
  int blocks = nbeams*nrays/threads_per_block + 1;
  printf("%d\n", blocks);
  printf("nindices %d\n", nindices);
  cudaDeviceSynchronize();
  int sense_rev = 0;
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
    getCBETGain(wMult, wMult_old, i_b, eden, areas, machnum, boxes, kvec, medianDS, marked, iteration, convergence);
    updateIntensities(i_b, wMult, iteration, convergence, currmax);
    #endif
    if(sense_rev % 1000 == 0) {
      sense_rev = 0;
      printf("%d: %e\n", iteration, *convergence);
    }
    double currmaxa = maxIncr*pow(CBETCONVERGENCE, iteration);
    double currmaxb = CBETCONVERGENCE*(*convergence);
    //printf("%e\n", *convergence);
    cudaMemcpy(wMult_old, wMult, sizeof(double)*nrays*ncrossings*nbeams, cudaMemcpyDeviceToDevice);

    //currmax = fmin(currmaxa, currmaxb);
  }while(*convergence > tolerance);
  printf("%e %d\n", iteration, *convergence);
  cudaDeviceSynchronize();
}