#include "CBET_Interface.hpp"
#define AVERAGEPUMP 0


//Different approach, break into unit-testable subroutines (inefficient, but should work)
#define NORM 1e14
double getInteractionMultiplier(int beam, int ray, int crossing, int ix, int iz, int islast)
{
  //calculate the interaction multiplier given the ray areas, the spatial position, the ray index, and the beam index
  double area1, area2;
  //get the average area of the ray across the zone
  area1 = vec3D(areas, beam, ray, crossing+!islast, nrays, ncrossings);
  area2 = vec3D(areas, beam, ray, crossing, nrays, ncrossings);
  double areaAvg = (area1+area2)/2;
  double neOverNc = eden[ix*nz + iz]/ncrit;//NOTE: This neOverNc value can be taken straight from the grid
  double neOverNcCorrected = fmin(neOverNc, 1.0);//clamp the maximum neOverNc value to
  double neTerm = sqrt(1-neOverNcCorrected);//term used multiple times in calculations
  double epsilonEff = neTerm*neTerm;
  double interactionMult = 1/(areaAvg*neTerm)*1/sqrt(epsilonEff);
  return interactionMult;
}
double getConstant()//calculate the CBET constant based on fundamental physical values, normalized to some relative intensity. Can easily be moved to top of main routine
{
  double e = 4.8032e-10;
  me = 9.109384e-28;
  kb = 1.38065e-16;
  c=2.99792458e10;
  double cbetconst = (8*pi*1e7*NORM/c)*(e*e/(4*me*c*kb*1.1605e7))*1e-4;
  return cbetconst;
}


//get the crossing index of ray passing through ix iz. NOTE: WILL SEGFAULT IF IT DOES NOT
int getCrossing(int beam, int ray, int ix, int iz)
{
  int rayCross = 0;
  int id = vec3D(boxes, beam, ray, rayCross, nrays, ncrossings);//spatial index of the ray
  int xcurr = (id - 1) % nx; //x index (from MATLAB array scheme, will change based on C array structure)
  int zcurr = (id - xcurr - 1)/nx;//z index (from MATLAB array scheme, will change based on C array structure)
  while (xcurr != ix || zcurr != iz)//while in a different cell
  {
    rayCross += 1;//increase rayCross
    if(rayCross >= ncrossings)//should not occur, once all of the bugs are ironed out->can be removed to reduce branching
    {
      printf("Error: Ray %d did not pass through (%d, %d)\n", ray, ix, iz);
      return -1;
    }
    id = vec3D(boxes, beam, ray, rayCross, nrays, ncrossings);//new spatial index
    xcurr = (id - 1) % nx;//update indices
    zcurr = (id - xcurr - 1)/nx;
  };
  return rayCross;
}

//Calculate the eta term in the CBET transfer equation: determines the sign of the exchange. Depends on difference in ray trajectories, frequency, and plasma velocity
double getEta(int beam_seed, int ray_seed, int cross_seed, int beam_pump, int ray_pump, int cross_pump, int ix, int iz)
{
  double mag1, mag2;//get the magnitudes of the rays, can calculate here from dkx and dkz instead of storing explicitly
  mag1 = vec3D(dkmag, beam_seed, ray_seed, cross_seed, nrays, ncrossings);
  mag2 = vec3D(dkmag, beam_pump, ray_pump, cross_pump, nrays, ncrossings);
  double kx_seed = vec3D(dkx, beam_seed, ray_seed, cross_seed, nrays, ncrossings);//get the x and y components of the ray vectors
  double kz_seed = vec3D(dkz, beam_seed, ray_seed, cross_seed, nrays, ncrossings);

  double kx_pump = vec3D(dkx, beam_pump, ray_pump, cross_pump, nrays, ncrossings);
  double kz_pump = vec3D(dkz, beam_pump, ray_pump, cross_pump, nrays, ncrossings);

  double machx = machnum[ix*nz+iz];//the mach number of the plasma velocity
  double neOverNc = eden[ix*nz + iz]/ncrit;//Interpolated neOverNc value


  double machz = 0.0;//TODO: Implement multidimensional plasma velocity;
  double omega1 = omega, omega2 = omega; //assuming uniform wavelength/frequency
  //find ion acoustic wave vector, difference of ray trajectories scaled to omega/(sqrt(1-neOverNc)*c)

  double iawVector[] = {(omega1*kx_seed - omega2*kx_pump)*sqrt(1-neOverNc)/c,(omega1*kz_seed - omega2*kz_pump)*sqrt(1-neOverNc)/c};
  double k_iaw = sqrt(iawVector[0]*iawVector[0] + iawVector[1]*iawVector[1]);//magnitude of the iaw vector
  double etaNumerator = omega1-omega2 - (iawVector[0]*machx + iawVector[1]*machz)*cs;//numerator of eta
  double etaDenominator = k_iaw*cs;//denominator of eta
  double eta = etaNumerator/etaDenominator;//find eta
  return eta;
};

//Get the actual interaction multiplier based on precalculated values
double getCouplingMultiplier(int beam, int ray, int cross, int ix, int iz, double cbetconst, double eta, double interactionMult)
{
  double neOverNc = eden[ix*nz + iz]/ncrit;//neOverNc taken from grid
  double param1 = cbetconst/(omega*(Te_eV/1e3 + 3 * Ti_eV/1e3/Z)); //Split coupling multiplier into discrete chuncks->easier debugging
  double param2 = neOverNc/iaw*iaw*iaw*eta;
  double param3 = (pow(eta*eta - 1, 2) + iaw*iaw*eta*eta);
  double param4 = interactionMult;
  double couplingMult = param1*param2/param3*param4;
  return couplingMult;
}
//limit the multiplier to ensure convergence. The limit on the multiplier increases over the course of the simulation, but the limit on the energy change will decrease to compensate
#define ABSLIM 1000
double limitMultiplier(double mult, int beam, int ray, int crossing, double medianDS, int iteration)
{

  double ds = vec3D(dkmag, beam, ray, crossing, nrays, ncrossings);//ray path length in cell
  double stepRatio = ds/medianDS;//ratio of path length to median
  double multPerStep = (1-mult)/stepRatio;//ratio of multiplier to step length
  double multiplierLim = maxIncr*iteration;//maximum allowable multiplier for iteration
  multPerStep = (multPerStep > multiplierLim) ? multiplierLim : multPerStep;//abs(mult) greater than multiplier limit, clamp with appropriate sign
  multPerStep = (multPerStep < -1*multiplierLim) ? -1*multiplierLim : multPerStep;
  mult = 1 - multPerStep*stepRatio;//restore mult value
  mult = (mult > ABSLIM) ? ABSLIM : mult;//make sure mult is never above absolute maximum
  mult = (mult < 1/ABSLIM) ? 1/ABSLIM : mult;
  return mult;
};

//limit the energy change allowed in a given cell in a given iteration
double limitEnergy(double multiplier_acc, double i0, double currMax, int beam, int ray, int crossing, double* maxChange)
{
  double i_prev = vec3D(i_b, beam, ray, crossing, nrays, ncrossings);//previous value in cell
  double i_curr = i0*multiplier_acc;//"current" value, unclamped update
  double fractionalChange = fabs(i_curr-i_prev)/i_prev;//the fractional change in energy from imposing the update as is
  //printf("(%e) ", i0);
  *maxChange = fmax(fractionalChange, *maxChange);//update the convergence check variable
  double correction = 0.0;
  if(fractionalChange > currMax)//If the fractional change is too large, clamp the value
  {
    int sign = (i_curr - i_prev > 0) ? 1 : -1;
    correction = 1 + currMax*sign;
    i_curr = i_prev*correction;
  }
  return i_curr;
}
//main CBET routine
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
          int raycross = getCrossing(q, ray_o, ix, iz);//get the crossing index of the other ray
          int islastq = !(vec3D(boxes, q, ray_o, raycross+1, nrays, ncrossings));//check if this is the last crossing for ray_o
          double interactionMultiplier = getInteractionMultiplier(q, ray_o, raycross, ix, iz, islast);//find ray_o's interaction multiplier
          double eta = getEta(i, j, k, q, ray_o, raycross, ix, iz);//calculate eta from the interacting rays
          double couplingMult = getCouplingMultiplier(i, j, k, ix, iz, cbetConst, eta, interactionMultiplier);//get the coupling multiplier
          couplingMult *= vec3D(dkmag, i,  j, k, nrays, ncrossings);//scale the coupling multiplier by the magnitude of the seed ray's path length
          double otherIntensity1 = vec3D(i_b, q, ray_o, raycross, nrays, ncrossings);
          double otherIntensity2 = vec3D(i_b, q, ray_o, raycross + !islastq, nrays, ncrossings);
          double avgIntensity = (otherIntensity2 + otherIntensity1)/2;//average the intensity of the other ray across the cell
          cbetSum += couplingMult*avgIntensity;//add to the accumulator
          vec3DW(spatialLog, i,  j, k, nrays, ncrossings, couplingMult);//logging variable
        }
        double mult = exp(-1*cbetSum);//exponentiate the sum
        mult = limitMultiplier(mult, i, j, k, medianDS, iteration);//limit the multiplier
        vec3DW(wMult, i, j, k, nrays, ncrossings, mult);//store the limited value
      }
    }
  }
  getchar();
}
*/

void updateIntensities(int iteration, double* convMax, double currMax)
{
  for(int i = 0; i < nbeams; i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      double i0 = vec3D(i_b, i, j, 0, nrays, ncrossings);
      double multAcc = 1.0;

      for(int k = 1; k < ncrossings; k++)
      {
        double mult = vec3D(wMult, i, j, k, nrays, ncrossings);
        double new_intensity = limitEnergy(multAcc, i0, currMax, i, j, k, convMax);
        multAcc *= mult;
        vec3DW(i_b, i, j, k, nrays, ncrossings, new_intensity);

        //vec3DW(spatialLog, i, j, k, nrays, ncrossings, new_intensity/1e14);

      }
    }
  }
}
#define CBETCONVERGENCE 0.9990
void cbet()
{
  spatialLog = new double[CROSS];
  initArrays();
  if(cudaCalc)
  {
    //launchCBETKernel();
    return;
  }
  int s = 0;
  if(s)
  {
    return;
  }
  wMult = new double[CROSS];//store the cumulative product of the normalized ray energies
  double* wMultOld = new double[CROSS];
  double* conv = new double[threads];
  for(int i = 0; i < CROSS; i++)
  {
    wMultOld[i] = 1.0;
    wMult[i] = 1.0;
  }
  double medianDS = median(dkmag, CROSS);
  double currmax = maxIncr;
  int raystore[nbeams*nx*nz]{0};
  for(int beam = 0; beam < nbeams; beam++)
  {
    for(int ix = 0; ix < nx; ix++)
    {
      for(int iz = 0; iz < nz; iz++)
      {
        int numpresent = vec3D(present, beam, ix, iz, nx, nz);
        if(!numpresent)
        {
          continue;
        }
        int ind = rand() % numpresent;
        raystore[(beam*nx + ix)*nz + iz] = vec4D(marked, beam, ix, iz, ind, nx, nz, numstored);
        //printf("%d %d %d\n", raystore[(beam*nx + ix)*nz + iz], );
      }
    }
  }
  auto startKernel = std::chrono::high_resolution_clock::now();
  int i;
  for(i = 1; i <= 500; i++)
  {
    double updateconv = 0.0;

    for(int j = 0; j < threads; j++)
    {
      conv[j] = 0;
    }
    getCBETGain(wMultOld, conv, spatialLog, medianDS, i, raystore);
    updateIntensities(i, &updateconv, currmax);
  //  printf("%f\n", updateconv);
    double convMax = 0.0;
    for(int i = 0; i < threads; i++)
    {
      convMax = fmax(convMax, conv[i]);
    }
    if(updateconv <= converge)
    {
      break;
    }
    double currmaxa = maxIncr*pow(CBETCONVERGENCE, i);
    double currmaxb = CBETCONVERGENCE*updateconv;
    currmax = fmin(currmaxa, currmaxb);
    printf("%f\n", currmax);
  }

    for(int i = 0; i < nbeams; i++)
    {
      for(int j = 0; j < nrays; j++)
      {
        for(int k = 0; k < ncrossings; k++)
        {
          int ix, iz;
          ix = vec3D(boxes, i, j, k, nrays, ncrossings) / nz;
          iz = vec3D(boxes, i, j, k, nrays, ncrossings) % nz;
          iz--;
          iz += nz*(iz == -1);
          if(!ix && !iz)
          {
            break;
          }
          //printf("%f ", vec3D(i_b, i, j, k, nrays, ncrossings)/1e14);

        }
       // printf("\n");
      }
     // printf("\n");
    }
  printf("%d\n", i);
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {
      ////printf("%f", vec2D(eden, i, j, nz)/ncrit);
    }
    ////printf("\n");
  }
  //double* fieldData = direct_map(i_b_new);
  delete [] conv;
    //elete [] wMult;
  delete [] wMultOld;
  auto stopKernel = std::chrono::high_resolution_clock::now();
  //*output << "CPUCBET " << threads <<" "<< nrays << " " << chrono::duration_cast<chrono::milliseconds>(stopKernel-startKernel).count() << std::endl;

}