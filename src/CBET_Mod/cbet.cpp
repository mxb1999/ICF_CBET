#include "CBET_Interface.hpp"

void cbetGain()
{
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
  printf("REe");
  for(int i = 0; i < nbeams;i++)
  {
    for(int j = 0; j < nrays;j++)
    {
      for(int m = 0; m < ncrossings;m++)
      {
        int ix = vec4D(boxes, i,j,m, 0, nrays, ncrossings, 2);
        int iz = vec4D(boxes, i,j,m, 0, nrays, ncrossings, 2);
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
        double epsilon = 1.0-ne/ncrit;
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
          for(int l = 0; l < lim; l++)
          {
            int rayCross = 0;
            int r = vec4D(marked,q,ix,iz,l, nx,nz,numstored)-1;
            for(int p = 0; p < ncrossings; p++)
            {
              int ox = vec4D(boxes, q,r,p, 0, nrays, ncrossings, 2);
              int oz = vec4D(boxes, q,r,p, 1, nrays, ncrossings, 2);
              if(!ox || !oz)
              {
                break;
              }
              ox--;
              oz--;
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
            double cumProd = vec3D(wMult, q, r, rayCross, nrays, ncrossings);

            
            double kmag = (omega/c)*sqrt(epsilon);
            double kx1 = kmag * vec3D(dkx, i, j, m, nrays, ncrossings) / mag1;
            double kx2 = kmag * vec3D(dkx, q, r, rayCross, nrays, ncrossings) / mag2;
            
            double kz1 = kmag * vec3D(dkz, i, j, m, nrays, ncrossings) / mag1;
            double kz2 = kmag * vec3D(dkz, q, r, rayCross, nrays, ncrossings) / mag2;
            double kiaw = sqrt(pow(kx2-kx1,2.0)+pow(kz2-kz1,2.0));
            
            double ws = kiaw*cs;
            
            double omega2 = omega;
            double eta = ((omega2-omega1)-(kx2-kx1)*vec2D(u_flow,ix,iz,nz))/(ws+1.0e-10);
            double efield2 = sqrt(8.*pi*1.0e7*vec3D(i_b, q, r, rayCross, nrays, ncrossings)/c);   
            double P = (pow(iaw,2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw),2)*pow(eta,2));  
            double gain1 = constant1*pow(efield2,2)*(ne/ncrit)*(1/iaw)*P/icnt;               //L^-1 from Russ's paper
            int sign = (i < q) ? -1 : 1;
            double oldEnergy2 = vec3D(wMult, q,r,rayCross,nrays, ncrossings);
            double newEnergy1Mult = exp(oldEnergy2*mag1*gain1/sqrt(epsilon));
            newEnergy1Mult = fmin(newEnergy1Mult, 1.1);
            if(newEnergy1Mult < 10)
            {
              vec3DM(wMult, i, j, m, nrays, ncrossings,newEnergy1Mult);
              
            }
          }
            double curr = vec3D(wMult, i, j, m, nrays, ncrossings);
          
          
        }
      }
    }
  }
}
void cbetUpdate()
{
  for(int i = 0; i < nbeams;i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      double i0 = vec3D(i_b, i,j,0, nrays, ncrossings);
      for(int m = 1; m < ncrossings;m++)
      {
        double mult = vec3D(wMult, i,j,m, nrays, ncrossings);
        vec3DW(i_b, i,j,m, nrays, ncrossings, i0*mult);
        i0 = vec3D(i_b, i,j,m, nrays, ncrossings);

      }
    }
  }
}
void cbetUpdateFinal()
{
  for(int i = 0; i < nbeams;i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      double i0 = vec3D(i_b_new, i,j,0, nrays, ncrossings);
      for(int m = 1; m < ncrossings;m++)
      {
        double mult = vec3D(wMult, i,j,m, nrays, ncrossings);
        vec3DW(i_b_new, i,j,m, nrays, ncrossings, i0*mult);
        if(i == 1 && j == 19)
        {
          printf("beam %d %d %d mult %e %e\n", i,j,m, i0, mult);

        }
        i0 = vec3D(i_b_new, i,j,m, nrays, ncrossings);
      }
    }
  }
}
void cbet()
{
  initArrays();
  if(cudaCalc)
  {
    launchCBETKernel();
  }else{
    wMult = new double[CROSS];//store the cumulative product of the normalized ray energies
    for(int i = 0; i < CROSS; i++)
    {
      wMult[i] = 1.0;
    }
    for(int i = 0; i < 10; i++)
    {
      printf("%d\n",i);
      cbetGain();
      cbetUpdate();
    }
    cbetUpdateFinal();
  }
}