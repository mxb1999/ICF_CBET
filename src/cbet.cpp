#include "CBET_Interface.hpp"
#define AVERAGEPUMP 0


//Different approach, break into unit-testable subroutines (inefficient, but should work)
#define NORM 1e14
double getInteractionMultiplier(int beam, int ray, int crossing, int ix, int iz, int islast, double neOverNc)
{
  //calculate the interaction multiplier given the ray areas, the spatial position, the ray index, and the beam index
  double area1, area2;
  //get the average area of the ray across the zone
  area1 = vec3D(areas, beam, ray, crossing+1, nrays, ncrossings);
  area2 = vec3D(areas, beam, ray, crossing, nrays, ncrossings);
  double areaAvg = (area1+area2)/2;
  neOverNc = eden[ix*nz + iz]/ncrit;//NOTE: This neOverNc value can be taken straight from the grid
  double neOverNcCorrected = fmin(neOverNc, 1.0);//clamp the maximum neOverNc value to
  double neTerm = sqrt(1-neOverNcCorrected);
  double epsilonEff = neTerm*neTerm;
  double interactionMult = 1/(areaAvg*neTerm)*1/sqrt(epsilonEff);

  return interactionMult;
}
double getConstant()
{
  double e = 4.8032e-10;
  me = 9.109384e-28;
  kb = 1.38065e-16;
  c=2.99792458e10;
  double cbetconst = (8*pi*1e7*NORM/c)*(e*e/(4*me*c*kb*1.1605e7))*1e-4;
  return cbetconst;
}
int* findInteractions(int ix, int iz, int otherbeam, int* size)
{
  int i, j;
  int size_local = vec3D(present, otherbeam, ix, iz, nx, nz);
  int* result = new int[size_local];
  *size = size_local;
  for(int q = 0; q < size_local; q++)
  {
    int ray = vec4D(marked, otherbeam, ix, iz, q, nx, nz, numstored) - 1;
    result[q] = ray;
  }
  return result;
}
int getCrossing(int beam, int ray, int ix, int iz)
{
  int xcurr, zcurr;
  int rayCross = 0;
  //printf("%d %d :: %d %d\n", ix, iz, xcurr, zcurr);
  //printf("\n%d:: %d %d\n", ray, ix, iz);
  while (xcurr != ix || zcurr != iz)
  {
    //printf("%d %d, ", xcurr, zcurr);
    //fflush(stdout);
    //getchar();
    rayCross += 1;
    int id = vec3D(boxes, beam, ray, rayCross, nrays, ncrossings);
    xcurr = (id - 1) % nx;
    zcurr = (id - xcurr - 1)/nx;

    //printf("\t%d %d :: %d %d\n", ix, iz, xcurr, zcurr);
  };

  return rayCross;
}
double getEta(int beam_seed, int ray_seed, int cross_seed, int beam_pump, int ray_pump, int cross_pump, int ix, int iz)
{
  double mag1, mag2;
  mag1 = vec3D(dkmag, beam_seed, ray_seed, cross_seed, nrays, ncrossings);
  mag2 = vec3D(dkmag, beam_pump, ray_pump, cross_pump, nrays, ncrossings);
  double kx_seed = vec3D(dkx, beam_seed, ray_seed, cross_seed, nrays, ncrossings);///mag1;
  double kz_seed = vec3D(dkz, beam_seed, ray_seed, cross_seed, nrays, ncrossings);///mag1;

  double kx_pump = vec3D(dkx, beam_pump, ray_pump, cross_pump, nrays, ncrossings);///mag2;
  double kz_pump = vec3D(dkz, beam_pump, ray_pump, cross_pump, nrays, ncrossings);///mag2;

  double machx = machnum[ix*nz+iz];
  double neOverNc = vec3D(neovernc, beam_pump, ray_pump, cross_pump, nrays, ncrossings);


  double machz = 0.0;//TODO: Implement multidimensional plasma velocity;
  double omega1 = omega, omega2 = omega;
  double iawVector[] = {(omega1*kx_seed - omega2*kx_pump)*sqrt(1-neOverNc)/c,(omega1*kz_seed - omega2*kz_pump)*sqrt(1-neOverNc)/c};
  double k_iaw = sqrt(iawVector[0]*iawVector[0] + iawVector[1]*iawVector[1]);
  double etaNumerator = omega1-omega2 - (iawVector[0]*machx + iawVector[1]*machz)*cs;
  double etaDenominator = k_iaw*cs;
  double eta = etaNumerator/etaDenominator;
  int ray_o = (beam_seed == 1) ? ray_seed + 16 + 1 : ray_seed+1;

  //printf("%f %f\n", (kx_seed), sqrt(1-neOverNc));

  return eta;
};
double getCouplingMultiplier(int beam, int ray, int cross, int ix, int iz, double cbetconst, double eta, double interactionMult, double neOverNc)
{
  neOverNc = eden[ix*nz + iz]/ncrit;
  double TeK = 2.0;
  double TiK = 1.0;
  Z = 3.100000000000000;
  omega = 5.366528681791605e15;
  double param1 = cbetconst/(omega*(TeK + 3 * TiK/Z));
  double param2 = neOverNc/iaw*iaw*iaw*eta;
  double param3 = (pow(eta*eta - 1, 2) + iaw*iaw*eta*eta);
  double param4 = interactionMult;
  double couplingMult = param1*param2/param3*param4;
  //printf("%f\n",couplingMult);
    //getchar();
  return couplingMult;
}
#define ABSLIM 1000
double limitMultiplier(double mult, int beam, int ray, int crossing, double medianDS, int iteration)
{
  double ds = vec3D(dkmag, beam, ray, crossing, nrays, ncrossings);
  double stepRatio = ds/medianDS;
  double multPerStep = (1-mult)/stepRatio;
  double multiplierLim = maxIncr*iteration;
  multPerStep = (multPerStep > multiplierLim) ? multiplierLim : multPerStep;
  multPerStep = (multPerStep < -1*multiplierLim) ? -1*multiplierLim : multPerStep;
  mult = 1 - multPerStep*stepRatio;
  mult = (mult > ABSLIM) ? ABSLIM : mult;
  mult = (mult < 1/ABSLIM) ? 1/ABSLIM : mult;

  return mult;
};
double limitEnergy(double multiplier_acc, double i0, double currMax, int beam, int ray, int crossing, double* maxChange)
{
  double i_prev = vec3D(i_b, beam, ray, crossing, nrays, ncrossings);
  double i_curr = i0*multiplier_acc;
  double fractionalChange = fabs(i_curr-i_prev)/i_prev;
  //printf("(%e) ", i0);
  *maxChange = fmax(fractionalChange, *maxChange);
  double correction = 0.0;

  if(fractionalChange > currMax)
  {
    int sign = (i_curr - i_prev > 0) ? 1 : -1;
    correction = 1 + currMax*sign;
    i_curr = i_prev*correction;
  }

  //printf("{%d %e}, ", crossing, (i_curr-i0*multiplier_acc)/1e14);
  return i_curr;
}
void getCBETGain(double* wMultOld, double* conv, double* logger, double medianDS, int iteration, int* raystore)
{
  double cbetConst = getConstant();

  int i, j, k, q, p;

  for(i = 0; i < nbeams; i++)
  {
    for(j = 0; j < nrays; j++)
    {
      int index = (i != 1) ? j + 1 : j + 360;
      //printf("%d: ", index);
      for(k = 0; k < ncrossings-1; k++)
      {
        int ix, iz;
        int id = vec3D(boxes, i, j, k, nrays, ncrossings);
        ix = (id - 1) % nx;
        iz = (id - ix - 1)/nx;
        //printf(":: %d %d %d| ", ix, iz, (iz)*nz + ix+1);
       // fflush(stdout);
        if(!id)
        {
          break;
        }

        int islast = 0;
        if(k != ncrossings - 1)
        {
          int ixnext, iznext;
          int idnext = vec3D(boxes, i, j, k+1, nrays, ncrossings);

          ixnext = (idnext - 1) % nx;
          iznext = (idnext - ixnext - 1)/nx;
          if(!ixnext || !iznext)
          {
            islast = 1;
          }
        }else
        {
          islast = 1;
        }
        double neOverNc = vec3D(neovernc, i, j, k, nrays, ncrossings);
        double cbetSum = 0.0;
        double ds = vec3D(dkmag, i, j, k, nrays, ncrossings);

        for(q = 0; q < nbeams; q++)
        {
          if(q == i)
          {
            continue;
          }
          int numInteractions;
          int* interactions = findInteractions(ix, iz, q, &numInteractions);
          if(!numInteractions)
          {
            continue;
          }
          int ray_o = raystore[(q*nx + ix)*nz + iz];
          //int other_index = vec3D(interactions_ML, i, j, k, nrays, ncrossings);
          //int ray_o = interactions[0];//other_index - nrays*q - 1;
          int raycross = getCrossing(q, ray_o, ix, iz);
          int islastq = !(vec3D(boxes, q, ray_o, raycross, nrays, ncrossings));
          double interactionMultiplier = getInteractionMultiplier(q, ray_o, raycross, ix, iz, islast, neOverNc);
          double area1 = vec3D(areas, q, ray_o, raycross+!islast, nrays, ncrossings);
          double area2 = vec3D(areas, q, ray_o, raycross, nrays, ncrossings);
          double eta = getEta(i, j, k, q, ray_o, raycross, ix, iz);
          double areaAvg = (area1+area2)/2;
          double couplingMult = getCouplingMultiplier(i, j, k, ix, iz, cbetConst, eta, interactionMultiplier, neOverNc);

          couplingMult *= vec3D(dkmag, i,  j, k, nrays, ncrossings);
          double otherIntensity1 = vec3D(i_b, q, ray_o, raycross, nrays, ncrossings);
          double otherIntensity2 = vec3D(i_b, q, ray_o, raycross + !islastq, nrays, ncrossings);
          double avgIntensity = (otherIntensity2 + otherIntensity1)/2;
          cbetSum += couplingMult*avgIntensity;

          vec3DW(spatialLog, i,  j, k, nrays, ncrossings, couplingMult);
          //printf("{%d, %f} ",o_index, eta);
          //if(i == 1 && j == 15 && k == 5)
          // {
            //printf("{%d %f %f}, ", ray_o, otherIntensity/1e14, couplingMult, cbetSum);
        //  }


        }

        double mult = exp(-1*cbetSum);

        //if(j + 1 + 16*i == 32)
        //  {
        //    printf("{%d %f}, ", k + 1, mult);
        //  }
        //printf("{%d %e}, ", k+1, mult);
        mult = limitMultiplier(mult, i, j, k, medianDS, iteration);
        vec3DW(wMult, i, j, k, nrays, ncrossings, mult);

      }
      if(i == 0 && index == 1)
      {
        printf("\n");
      }

    }
   // printf("\n");

  }
    //printf("\n");

}


void updateIntensities(int iteration, double* convMax, double currMax)
{
  for(int i = 0; i < nbeams; i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      //printf("%d: ", j + 1 + 16*i);
      double i0 = vec3D(i_b, i, j, 0, nrays, ncrossings);
      double multAcc = 1.0;
      //vec3DW(spatialLog, i, j, 0, nrays, ncrossings, i0/1e14);

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
    launchCBETKernel();
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
        raystore[(beam*nx + ix)*nz + iz] = vec4D(marked, beam, ix, iz, ind, nx, nz, numstored) - 1;
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