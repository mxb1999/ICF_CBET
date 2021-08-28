#include "CBET_Interface.hpp"
#define AVERAGEPUMP 0


//Different approach, break into unit-testable subroutines (inefficient, but should work)
#define NORM 1e15
double getInteractionMultiplier(int beam, int ray, int crossing, int ix, int iz, int islast)
{
  //calculate the interaction multiplier given the ray areas, the spatial position, the ray index, and the beam index
  double area1, area2;
  int offset = !islast;
  area1 = vec3D(areas, beam, ray, crossing+offset, nrays, ncrossings);
  area2 = vec3D(areas, beam, ray, crossing, nrays, ncrossings);
  double areaAvg = (area1+area2)/2;
  double neOverNc = eden[ix*nz+iz]/ncrit;
  double neTerm = sqrt(1-neOverNc);
  double epsilonEff = neTerm*neTerm;
  double interactionMult = 1/(areaAvg*neTerm)*1/sqrt(epsilonEff);


  return interactionMult;
}
double getConstant()
{
  double e = 4.8032e-10;
  me = 9.1094e-28;
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
  int id = vec3D(boxes, beam, ray, rayCross, nrays, ncrossings);
  xcurr = (id - 1) % nx;
  zcurr = (id - xcurr - 1)/nx;
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
  double ne = eden[ix*nz + iz];
  double neOverNc = ne/ncrit;
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
double getCouplingMultiplier(int beam, int ray, int cross, int ix, int iz, double cbetconst, double eta, double interactionMult)
{
  double ne = eden[ix*nz + iz];
  double neOverNc = ne/ncrit;
  double TeK = Te_eV/1e3;
  double TiK = Ti_eV/1e3;
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
  double fractionalChange = abs(i_curr-i_prev)/i_prev;
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
  int i, j, k, q, p;

  for(i = 0; i < nbeams; i++)
  {
    for(j = 0; j < nrays; j++)
    {
      int index = (i != 1) ? j + 1 : j + 17;
      //printf("%d: ", index);
      for(k = 0; k < ncrossings; k++)
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
          //int ray_index = raystore[(i*nx + ix)*nz + iz];
          int ray_o = raystore[(q*nx + ix)*nz + iz];
          int ray_oid = (q == 1) ? ray_o + 16 + 1 : ray_o+1;
          int raycross = getCrossing(q, ray_o, ix, iz);
          double interactionMultiplier = getInteractionMultiplier(q, ray_o, raycross, ix, iz, islast);
          double cbetConst = getConstant();
          double eta = getEta(i, j, k, q, ray_o, raycross, ix, iz);
          double couplingMult = getCouplingMultiplier(i, j, k, ix, iz, cbetConst, eta, interactionMultiplier);
          couplingMult *= vec3D(dkmag, i,  j, k, nrays, ncrossings);

          double otherIntensity = vec3D(i_b, q, ray_o, raycross, nrays, ncrossings);
          printf("{%d %f}, ", ray_o, otherIntensity/1e14);
          cbetSum += couplingMult*otherIntensity/NORM;
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
      //printf("\n\n");


      printf("\n\n");
    }
   // printf("\n");

  }
    printf("\n");

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
      for(int k = 1; k < ncrossings; k++)
      {
        double mult = vec3D(wMult, i, j, k, nrays, ncrossings);
        multAcc *= mult;
        double new_intensity = limitEnergy(multAcc, i0, currMax, i, j, k, convMax);
        vec3DW(i_b, i, j, k, nrays, ncrossings, new_intensity);
      }
    }
  }
  getchar();
}
void cbetGain(double* wMultOld, double* conv, double* logger, double medianDS, int iteration)
{
  me = 9.1094e-31;
  double constant1 = (pow(estat,2.0))/(4*(1.0e3*me)*c*omega*kb*Te*(1+3*Ti/(Z*Te)));
  double coullim1 = Ti*me/mi;
  double coullim2 = 10*pow(Z,2);
  int logswitch = 0;
  if(Te < coullim1)
  {
    logswitch = 1;
  }else if(coullim1 < Te && Te < coullim2)
  {
    logswitch = 2;
  }else if(coullim1 < coullim2 && coullim2 < Te)
  {
    logswitch = 3;
  }
  if(logswitch == 0)
  {
    printf("Error in inputs\n");
    exit(0);
  }
  for(int i = 0; i < nbeams;i++)
  {
    //#pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nrays;j++)
    {
      int thisThread = omp_get_thread_num();
      double i0 = vec3D(i_b, i,j,0,nrays,ncrossings);
      for(int m = 0; m < ncrossings;m++)
      {
        int ix = vec4D(boxes, i,j,m, 0, nrays, ncrossings, 2);
        int iz = vec4D(boxes, i,j,m, 1, nrays, ncrossings, 2);
        if(!ix || !iz)
        {
          break;
        }
        ix--;
        iz--;
        int icnt = vec3D(present, i, ix,iz,nx,nz);
        double omega1 = omega;
        double mag1 = vec3D(dkmag, i, j, m, nrays, ncrossings);
        double ne = vec2D(eden, ix,iz,nz);
        double epsilon = 1/sqrt(1.0-ne/ncrit);
        double coullog;//coulomb logarithm for electron-electron collisions
        if(logswitch == 1)
        {
          double mp_kg = 1.6726219e-27;//mass of proton in kg
          coullog = 30-log(sqrt(ne)/pow(Ti, 3/2)*pow(Z, 3/2)*mi_kg/mp_kg);
        }else if(logswitch == 2)
        {
          coullog = 23-log(sqrt(ne)/pow(Ti, 3/2)*Z);
        }else
        {
          coullog = 24-log(sqrt(ne)/Ti);
        }
        double vei = 4*sqrt(2*pi)*pow(estat,4)*pow(Z,2)*coullog/(3*sqrt(me)*pow(Te, 3/2));
        double L_aij = c*sqrt(epsilon)*ncrit/(vei*ne);
        double prev = exp(-1*mag1/L_aij);
        vec3DW(wMult,i,j,m, nrays, ncrossings, exp(-1*mag1/L_aij));
        for(int q = 0; q < nbeams;q++)
        {
          if(q == i)
          {
            continue;
          }
          if(!vec4D(marked,q,ix,iz,0, nx,nz,numstored))
          {
            continue;
          }
          int qcnt = vec3D(present, q, ix,iz,nx,nz);
          //int qcross[qcnt]{0};
          int lim = fmin(qcnt, icnt);
          double summation = 0.0;
          for(int l = 0; l < lim; l++)
          {
            int rayCross = 0;
            int r = vec4D(marked,q,ix,iz,l, nx,nz,numstored)-1;
            double multAcc = vec3D(wMultOld, q,r,0, nrays, ncrossings);
            for(int p = 0; p < ncrossings-1; p++)
            {
              int ox = vec3D(boxes, q,r,p, nrays, ncrossings)/ nx;
              int oz = vec3D(boxes, q,r,p, nrays, ncrossings)% nx;
              if(!ox && !oz)
              {
                break;
              }
              ox--;
              oz--;
              multAcc*=vec3D(wMultOld, q,r,p, nrays, ncrossings);
              if(ox == ix && oz == iz)
              {
                rayCross = p;
                break;
              }
            }
            double mag2 = vec3D(dkmag, q, r, rayCross, nrays, ncrossings);
            if(mag2 < 1e-10)
            {
              continue;
            }
            double kmag = (omega/c)*sqrt(epsilon);
            double kx_seed = kmag * vec3D(dkx, i, j, m, nrays, ncrossings) / mag1;
            double kx_pump = kmag * vec3D(dkx, q, r, rayCross, nrays, ncrossings) / mag2;

            double kz_seed = kmag * vec3D(dkz, i, j, m, nrays, ncrossings) / mag1;
            double kz_pump = kmag * vec3D(dkz, q, r, rayCross, nrays, ncrossings) / mag2;
            double kiaw = sqrt(pow(kx_pump-kx_seed,2.0)+pow(kz_pump-kz_seed,2.0));

            double ws = kiaw*cs;
            double pumpI;
            if(AVERAGEPUMP)
            {
              pumpI = (vec3D(i_b, q, r, rayCross, nrays, ncrossings) + vec3D(i_b, q, r, rayCross + 1, nrays, ncrossings))/2;
              pumpI/=2;
            }else
            {
              pumpI = vec3D(i_b, q, r, rayCross, nrays, ncrossings);
            }
            double neOverNc = eden[ix*nz+iz]/ncrit;
            double e = 4.8032e-10;
            me = 9.1094e-28;
            double cbetconst = (8*pi*1e7*1e14/c)*(e*e/(4*me*c*kb*1.1605e7))*1e-4;
            double omega2 = omega;
            int offset = (rayCross != ncrossings);
            double averageArea = (vec3D(areas, q, r, rayCross+offset, nrays, ncrossings) + vec3D(areas, q, r, rayCross, nrays, ncrossings))/2;
            double epseff = 1/sqrt(1-neOverNc);
            double interactionMult = 1/(averageArea*sqrt(1-neOverNc))/sqrt(epseff);
            double eta = ((omega2-omega1)-(kx_seed-kx_pump)*vec2D(u_flow,ix,iz,nz))/(ws+1.0e-10);
            double couplingCoeff = cbetconst*1/omega*1/(Te_eV/1e3 + 3*Ti_eV/(1e3*Z))*neOverNc/iaw*iaw*iaw*eta/(pow(eta*eta-1, 2) + iaw*iaw*eta*eta)*interactionMult;
            double efield2 = 8.*pi*1.0e7*pumpI/c;
            double P = (pow(iaw,2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw),2)*pow(eta,2));
            double gain1 = constant1*efield2*(ne/ncrit)*(1/iaw)*P;               //L^-1 from Russ's paper
            double newEnergy1Mult = multAcc*mag1*gain1/sqrt(epsilon)/lim;
            summation += newEnergy1Mult;
          }
          double mult = exp(-1*mag1/L_aij + summation);
          double temp = mult;
          double dsratio = mag1/medianDS;
          double mult_per_step = (1-mult)/dsratio;
          double mult_limit = iteration*maxIncr;
          mult_per_step += (mult_per_step < -1*mult_limit)*-1*(mult_per_step + mult_limit) + (mult_per_step > mult_limit)*(-1*mult_per_step + mult_limit);
          mult = 1 - mult_per_step*dsratio;
          mult += (mult > 1000)*(1000.0-mult) + (mult < 1.0/1000)*(1.0/1000-mult);
          double prevVal = vec3D(wMultOld, i, j, m, nrays, ncrossings);
          conv[thisThread] = fmax(conv[thisThread], abs(mult-prevVal)/prevVal);
          //printf("(%f, %f) ", mult_per_step, mult);
          vec3DW(wMult, i, j, m, nrays, ncrossings, mult);
          //if(mult > 1.2 || mult < 0.8)
        }
      }
    //printf("\n");
    }
  }
  //getchar();
}

#define CBETCONVERGENCE 0.9990
void cbetUpdate(double* wMultOld, double medianDS, int iteration, double currmax, double* conv)
{

  for(int i = 0; i < nbeams;i++)
  {
    //#pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nrays; j++)
    {
      double multAcc = 1;
      double i0 = vec3D(i_b, i,j,0, nrays, ncrossings);
      for(int m = 1; m < ncrossings;m++)
      {
        double mult = vec3D(wMult, i,j,m, nrays, ncrossings);
        multAcc *= mult;
        int ix = vec3D(boxes, i,j,m, nrays, ncrossings) / nx;
        int iz = vec3D(boxes, i,j,m, nrays, ncrossings) % nx;
        if(!ix || !iz)
        {
          break;
        }
        iz--;
        double iprev = vec3D(i_b, i,j,m, nrays, ncrossings);
        double icurr = i0*multAcc;
        double temp = icurr;
        double frac_change = abs((icurr-iprev)/iprev);
        if(frac_change > currmax)
        {
          int sign = (icurr-iprev > 0) ? 1 : -1;
          double correction = (1+currmax*sign);
          icurr = iprev*correction;
        }
        *conv = fmax(frac_change, *conv);
        //printf("%f ", currmax);
        vec3DW(i_b, i,j,m, nrays, ncrossings, icurr);
        vec3DW(wMultOld, i,j,m, nrays, ncrossings, mult);
      }
      //printf("\n");
    }
    //printf("\n");
  }
  //getchar();
}
void cbetUpdateFinal()
{
  for(int i = 0; i < nbeams;i++)
  {
    //#pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nrays; j++)
    {
      double i0 = vec3D(i_b, i,j,0, nrays, ncrossings);
      for(int m = 1; m < ncrossings;m++)
      {
        double mult = vec3D(wMult, i,j,m, nrays, ncrossings);
        i0 *= mult;
        vec3DW(i_b, i,j,m, nrays, ncrossings, i0);
      }
    }
  }
}

//note: may be difficult to generalize to 3D with initial areas
//calculate the amplification of the ray amplitude due to CBET
double* cbet_amplification(double* multiplier, double* initial_areas, double* amplitudes)
{
    double* amplification = (double*)malloc(sizeof(double)*CROSS);
    int i, j, k;
    for(i = 0; i < nbeams; i++)
    {
        for(j = 0; j < nrays; j++)
        {
            double init_area = vec2D(initial_areas, i, j, nrays);
            for(k = 0; k < ncrossings; k++)
            {
                double imult = vec3D(multiplier, i, j, k, nrays, ncrossings);
                double iamp = vec3D(amplitudes, i, j, k, nrays, ncrossings)/1e14;
                double val = (imult-1)*init_area*iamp;
                vec3DW(amplification, i, j, k, nrays, ncrossings, val);
            }
        }
    }
    return amplification;
}

void cbet()
{
  double spatialLog[100];
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
  int raystore[nbeams*nx*nz];
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
  for(i = 1; i <= 100; i++)
  {
    double updateconv = 0.0;

    for(int j = 0; j < threads; j++)
    {
      conv[j] = 0;
    }
    getCBETGain(wMultOld, conv, spatialLog, medianDS, i, raystore);
    updateIntensities(i+1, &updateconv, currmax);
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
          printf("%f ", vec3D(i_b, i, j, k, nrays, ncrossings)/1e14);

        }
        printf("\n");
      }
      printf("\n");
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