#include "implSim.h"
#include "customMath.h"
#include <iomanip>
#include <cstring>
using namespace std;
//initializing arrays and constants for the CBET subroutines
void initializeArr()
{
  i_b = new double**[nbeams];
  for(int i = 0; i < nbeams; i++)
  {
    i_b[i] = new double*[nx];
    for(int j = 0; j < nx; j++)
    {
        i_b[i][j] = new double[nz]{0.0};
    }
  }
  cs = 1e2*sqrt(ec*(Z*Te_eV+3.0*Ti_eV)/mi_kg);	// acoustic wave speed, approx. 4e7 cm/s in this example
  for(int m = 0; m < nbeams; m++)
  {
    for(int i = 0; i < nx;i++)
    {
      for(int j = 0; j < nz;j++)
      {
        i_b[m][i][j] = edep[m][i][j];
        i_b_new[m][i][j] = edep[m][i][j];
        i_b_prev[m][i][j] = 0;
      //  i_b_prev[m][i][j] = edep[m][i][j] + ((m)%2 != 0)*(-1)*edep[m][i][j]*0.2 + ((m)%2 == 0)*edep[m][i][j]*0.2;
        //
        W[m][i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);
        W_new[m][i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);
        W_init[m][i][j] =   W[m][i][j];
        u_flow[i][j] = machnum[i][j]*cs;
      }

    }
  }
  //cout << scientific;
  for(int i = 0; i < nbeams;i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings-1; m++)
      {
        dkx[i][j][m] = crossesx[i][j][m+1]-crossesx[i][j][m];
        dkz[i][j][m] = crossesz[i][j][m+1]-crossesz[i][j][m];
        dkmag[i][j][m] = sqrt(pow(dkx[i][j][m],2.0)+pow(dkz[i][j][m],2.0));
      }
    }
  }
}


void calculateMult(int ix, int iz, int i, int j, int m, int* numrays, int** marker, int** markerCopy)
{
  double constant1 = (pow(estat,2.0))/(4*(1.0e3*me)*c*omega*kb*Te*(1+3*Ti/(Z*Te)));
  int n2limit = fmin(present[i][ix][iz],numrays[i+1]);
  for ( int n2 = 0; n2 < n2limit; n2++)
  {
    double ne = eden[ix][iz];
    double epsilon = 1.0-ne/ncrit;
    double kmag = (omega/c)*sqrt(epsilon);         // magnitude of wavevector
    double l_d = sqrt(e0*kb*Te/(pow(ec,2)*ne)); //Debye length of the plasma
    double N_d = ne*4*pi/3*pow(l_d,3); //number of electrons in the debye sphere
    double coulLog = log(9*N_d/Z); //Coulomb logarithm of the plasma
    double kx1 = kmag * dkx[i][j][m] / (dkmag[i][j][m] + 1.0e-10);
    double kx2 = kmag * dkx[i+1][marker[i+1][n2]][markerCopy[i+1][n2]] / (dkmag[i+1][marker[i+1][n2]][markerCopy[i+1][n2]] + 1.0e-10);
    double kz1 = kmag * dkz[i][j][m]/(dkmag[i][j][m]+1.0e-10);
    double kz2 = kmag * dkz[i+1][marker[i+1][n2]][markerCopy[i+1][n2]] / (dkmag[i+1][marker[i+1][n2]][markerCopy[i+1][n2]] + 1.0e-10);
    double kiaw = sqrt(pow(kx2-kx1,2.0)+pow(kz2-kz1,2.0));
  // magnitude of the difference between the two vectors
  //  double vei = 4*sqrt(2*pi)*pow(ec,4)*pow(Z,2),
    double ws = kiaw*cs;            // acoustic frequency, cs is a constant
    double omega1= omega;  // laser frequency difference. To start, just zero.
    double omega2 = omega;
    double eta = ((omega2-omega1)-(kx2-kx1)*u_flow[ix][iz])/(ws+1.0e-10);

   double efield1 = sqrt(8.*pi*1.0e7*i_b[i][ix][iz]/c);             // initial electric field of rayi;lnlni46
   double efield2 = sqrt(8.*pi*1.0e7*i_b[i+1][ix][iz]/c);             // initial electric field of ray

   double P = (pow((iaw),2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw),2)*pow(eta,2));         // From Russ's paper
   double gain1 = constant1*pow(efield2,2)*(ne/ncrit)*(1/iaw)*P;               //L^-1 from Russ's paper
   double gain2 = constant1*pow(efield1,2)*(ne/ncrit)*(1/iaw)*P;
            //L^-1 from Russ's paper

                    // new energy of crossing (PROBE) ray (beam 2)
    //updating arrays based upon the calculated magnitude changes
    if (dkmag[i+1][marker[i+1][n2]][markerCopy[i+1][n2]] >= 1.0*dx)
    {
      #pragma omp atomic write
      W_new[i+1][ix][iz] =W[i+1][ix][iz]*exp(-1*W[i][ix][iz]*dkmag[i+1][marker[i+1][n2]][markerCopy[i+1][n2]]*gain2/sqrt(epsilon));
      #pragma omp atomic write
      W_new[i][ix][iz] =W[i][ix][iz]*exp(1*W[i+1][ix][iz]*dkmag[i][j][m]*gain1/sqrt(epsilon));
    }
  }
}
//updating CBET Intensities
void update_CBETPar()
{
  //storage array to allow for data level parallelization
  //cout << "Updating CBET Intensities" << endl;

  double track1 = 0.0;
  double track2 = 0.0;
  //perform calculation for each crossing of each ray calculated for the beams in launch rays
  for(int i = 0; i < nbeams-1;i++)
  {
    #pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nrays; j++)
    {
      int endCheck = 0;
      for(int m = 0; m < ncrossings; m++)
      {
        if(boxes[i][j][m][0] == 0 || boxes[i][j][m][1] == 0)
        {
          break;
        }
        int ix = boxes[i][j][m][0] - 1;
        int iz = boxes[i][j][m][1] - 1;
        //if two beams intersect
        if(ix != -1 && iz != -1 && intersections[ix][iz] != 0)
        {


          //cout << "maxDev 1: " <<maxDev << " " << m << " " << j <<  endl;
          double tempStore = maxDev;
          //find nonzero locations of rays
          int numrays[nbeams]{0};
          int** marker = new int*[nbeams];
          int** markerCopy = new int*[nbeams];
          for(int q = 0; q < numstored; q++)
          {
            if(marked[ix*nz+iz][q*nbeams + 0] != 0)
            {
              numrays[0]++;
            }
            if(marked[ix*nz+iz][q*nbeams + 1] != 0)
            {
              numrays[1]++;
            }
          }
          marker[0] = new int[numrays[0]];
          marker[1] = new int[numrays[1]];
          markerCopy[0] = new int[numrays[0]];
          markerCopy[1] = new int[numrays[1]];
          int cnt1 = 0;
          int cnt2 = 0;
          for(int q = 0; q < numstored; q++)
          {
            if(marked[ix*nz+iz][q*nbeams + i + 1] != 0)
            {
              marker[1][cnt1] = marked[ix*nz+iz][q*nbeams + i + 1] - 1;
              markerCopy[1][cnt1] = marked[ix*nz+iz][q*nbeams + i + 1] - 1;
              cnt1++;
            }
            if(marked[ix*nz+iz][q*nbeams + i] != 0)
            {
              marker[0][cnt2] = marked[ix*nz+iz][q*nbeams + i] - 1;
              markerCopy[0][cnt2] = marked[ix*nz+iz][q*nbeams + i] - 1;
              cnt2++;
            }
          }

          int ray1num = 0;
          for(int r = 0; r < numrays[i]; r++)
          {
            if(marker[i][r] == j)
            {
              ray1num = r;
              break;
            }
          }
          //ix2 is
      //    //cout << markerCopy[1][0] << endl;
          for(int n = 0; n < numrays[i+1]; n++)
          {
            for(int q = 0; q < ncrossings;q++)
            {
              int ix2 = boxes[i+1][marker[i+1][n]][q][0] - 1;
              int iz2 = boxes[i+1][marker[i+1][n]][q][1] - 1;
              if(ix == ix2 && iz == iz2)
              {
                markerCopy[i+1][n] = q;
                break;
              }
            }
          }
          calculateMult(ix, iz, i, j, m, numrays, marker, markerCopy);
          //calculate the additive change in CBET field to be applied
          //cout << scientific;
          double additive_change[2];
          additive_change[0] = (-1.0 * (1.0 - (W_new[i][ix][iz]/W[i][ix][iz])) * abs(i_b[i][ix][iz] - i_b_prev[i][ix][iz])*(i_b[i][ix][iz] != 0));
          additive_change[1] = (-1.0 * (1.0 - (W_new[i+1][ix][iz]/W[i+1][ix][iz])) *  abs(i_b[i+1][ix][iz] - i_b_prev[i+1][ix][iz]));
          int kill = 0;
          if(i_b_new[i][ix][iz] < abs(additive_change[0]) && additive_change[0] < 0 || i_b_new[i][ix][iz] <= 0)
          {
            additive_change[1] = i_b_new[i][ix][iz] * (i_b_new[i][ix][iz] > 0);
            #pragma omp atomic write
            i_b_new[i][ix][iz] = 0;
            additive_change[0] = 0;
            kill = 1;
          }else
          {
            #pragma omp atomic update
            i_b_new[i][ix][iz] += additive_change[0];
            #pragma omp atomic update
            track1+=additive_change[0];
          }
          tempStore = fmax(abs(additive_change[0]/(i_b[i][ix][iz]+(i_b_new[i][ix][iz] <= 1e-10))), tempStore);
          tempStore = fmax(abs(additive_change[1]/(i_b[i+1][ix][iz]+(i_b_new[i+1][ix][iz] <= 1e-10))), tempStore);
        //  #pragma omp atomic update
      //    orderplot1[ix][iz] += 1*kill;//*(orderplot1[ix][iz] < 1);
        //  #pragma omp atomic write
    //      orderplot1[ix][iz] = 1;//*(orderplot1[ix][iz] < 1)
          #pragma omp atomic update
          i_b_new[i+1][ix][iz] += additive_change[1];
          #pragma omp atomic update
          track2+=additive_change[1];
          int iter;
          int n2 = fmin(ray1num, numrays[i+1] - 1);

          double x_prev = x[ix];
          double z_prev = z[iz];
          for(int l = m+1; l < ncrossings;l++)
          {
            int ix_next_i = boxes[i][j][l][0] - 1;
            int iz_next_i = boxes[i][j][l][1] - 1;
            if(ix_next_i != -1 && iz_next_i != -1)
            {
              double x_curr = x[ix_next_i];
              double z_curr = z[iz_next_i];
              if(x_curr != x_prev || z_curr != z_prev)
              {
                double add = additive_change[0]*present[i][ix][iz]/present[i][ix_next_i][iz_next_i];
                double comp = add/(i_b_new[i][ix_next_i][iz_next_i]+(i_b_new[i][ix_next_i][iz_next_i] < 1e-10));
                tempStore = fmax(tempStore, comp);
              //  cout << "3: " << tempStore << endl;
            //  cout << "Beam 1 Downstream Original: "<< i_b_new[i][ix_next_i][iz_next_i] << endl;
              //  #pragma omp atomic update
            //    orderplot1[ix_next_i][iz_next_i] += 1*kill;//*(orderplot1[ix_next_i][iz_next_i] < 1);
;
                #pragma omp atomic update
                i_b_new[i][ix_next_i][iz_next_i] += add + kill*(-1)*(i_b_new[i][ix_next_i][iz_next_i] + add);
                #pragma omp atomic update
                i_b[i][ix_next_i][iz_next_i] += kill*(-1)*i_b[i][ix_next_i][iz_next_i];
                #pragma omp atomic write
                orderplot1[ix_next_i][iz_next_i] = 1;//*(orderplot1[ix_next_other][iz_next_other] < 1);
              //  cout << "Downstream for Beam 1: " << add<<endl;
            //    cout << "Beam 1 Downstream Final: "<< i_b_new[i][ix_next_i][iz_next_i] << endl;

              //  #pragma omp atomic update
                //track1+=add;
              }
              x_prev = x_curr;
              z_prev = z_curr;
            }
          }
          for(int l = markerCopy[i+1][n2] + 1; l < ncrossings;l++)
          {
            int ix_next_other = boxes[i+1][marker[i+1][n2]][l][0] - 1;
            int iz_next_other = boxes[i+1][marker[i+1][n2]][l][1] - 1;
            ////cout << i+1 << " " << n2 << " "<< marker[i+1][n2] <<" " << l  << " " << endl;

            if(ix_next_other != -1 && iz_next_other != -1)
            {
              double x_curr = x[ix_next_other];
              double z_curr = z[iz_next_other];
              if(x_curr != x_prev || z_curr != z_prev)
              {
                double add = additive_change[1]*present[i][ix][iz]/present[i+1][ix_next_other][iz_next_other];
                double comp = add/(i_b_new[i+1][ix_next_other][iz_next_other]+(i_b_new[i+1][ix_next_other][iz_next_other] < 1e-10));
                //cout << ix_next_other << " " << iz_next_other << " " << additive_change[1] * (double)present[i][ix][iz]/present[i+1][ix_next_other][iz_next_other] << endl;
                tempStore = fmax(tempStore, comp);
            //    cout << "4: " << tempStore << endl;
                //cout << "Beam 2 Downstream Original: "<< i_b_new[i+1][ix_next_other][iz_next_other] << endl;


                #pragma omp atomic update
                i_b_new[i+1][ix_next_other][iz_next_other] += add;

              //  cout << "Downstream for Beam 2: " << add<<endl;
              //  cout << "Beam 2 Downstream Final: "<< i_b_new[i+1][ix_next_other][iz_next_other] << endl;

                ////cout << "Downstream raise for beam " << i+1 << " and Ray " << markerCopy[i+1][n2] << " at " << ix_next_other << " " << iz_next_other << endl;
              //  #pragma omp atomic update
              //  track2 += add;
              }
              x_prev = x_curr;
              z_prev = z_curr;
            }
          }
        //  cout << "Beam 1 Final: " << i_b_new[i][ix][iz] << endl;
    //      cout << "Beam 2 Final: " << i_b_new[i+1][ix][iz] << endl;
        //  cout << endl;
          #pragma omp atomic write
          maxDev = tempStore;
        //  cout << "maxDev 2: " <<maxDev << " " << m << " " << j << endl;
        }
      }

      if ( j % 20 == 0 )
      {
        ////cout << "     ..."<<(int)(100.*(1.0-(double)(j/double(1*nrays))))<<"%  remaining..."<<endl;
      }
    }
    cout << "maxDev i: " <<maxDev << " "  << i <<  endl;

  }
  //cout << "Finished Update loop" << endl;
  double sum1 = 0;
  double sum2 = 0;
  for(int i = 0; i < nz; i++)
  {
    for(int j = 0; j < nx; j++)
    {
      sum1+=i_b_new[0][j][i]+i_b_new[1][j][i];
      sum2+=edep[0][j][i]+edep[1][j][i];
    }
  }
  cout << "ENERGY CONSERVATION: " << sum2 << "  "<<sum1 << " " << 100*(sum1-sum2)/sum2<< endl;
  cout << "______TRACK 1: " << track1 << endl;
    cout << "______TRACK 2: " << track2 << endl;
for(int i = 0; i < nbeams; i++)
{
  #pragma omp parallel for num_threads(threads)
  for(int j = 0; j < nx; j++)
  {
    for(int m = 0; m < nz; m++)
    {
      W_init[i][j][m] = W_new[i][j][m];
      W[i][j][m] = W_new[i][j][m];
      i_b_prev[i][j][m] = i_b[i][j][m];
      i_b[i][j][m] = i_b_new[i][j][m];
    }
  }
}
}
void cbet()
{
  i_b_prev = new double**[nbeams];
  iter = 0;

  orderplot1 = new double*[nrays];
  for(int i = 0; i < nrays; i++)
  {
    orderplot1[i] = new double[nrays]{0};
  }
  i_b_new = new double**[nbeams];
  double*** icopy = new double**[nbeams];
  for(int i = 0; i < nbeams; i++)
  {
    icopy[i] = new double*[nx];
    i_b_prev[i] = new double*[nx];
    i_b_new[i] = new double*[nx];
    for(int j = 0; j < nx; j++)
    {
      icopy[i][j] = new double[nz];
      i_b_prev[i][j] = new double[nz]{0.0};
      i_b_new[i][j] = new double[nz]{0.0};
    }
  }
  int cntMax = 0;
  //initializeArr();
//  gain_CBETSeq();
//  update_CBETSeq();
int rayorder[nrays];
for(int i = 0; i < nrays; i++)
{
  rayorder[i] = i;
}
//cout << "nRays = " << nrays << endl;
maxDev = 0;
double prev = maxDev;
initializeArr();

auto start = chrono::high_resolution_clock::now();
//gain_CBETPar();
auto inter = chrono::high_resolution_clock::now();
update_CBETPar();
//cout << "MaxDev: " << maxDev << endl;
int cnt = 0;
cout << "MaxDev Preloop: " << maxDev << endl;
/*
for(int i = 0; i < nrays; i++)
{
  int r1 = ((int)rand())%nrays;
  int r2 = ((int)rand())%nrays;
  int temp = rayorder[r1];
  rayorder[r1] = rayorder[r2];
  rayorder[r2] = temp;
}
*/
 while(abs(maxDev) > converge)
  {
    iter++;
    prev = maxDev;
    maxDev = 0;
    //gain_CBETPar();
    update_CBETPar();
    cout << "MaxDev Postupdate: " << maxDev << endl;
    cnt++;
  }
  cout << "Number of iterations: " << iter << endl;

  i_b1Error = new double*[nx];
  i_b2Error = new double*[nx];
  /*
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {
      cout << i_b_new[0][i][j] << endl;
    }
  }
  */
  i_bplot = new double*[nz];
  i_b_newplot = new double*[nx];
  edepplot = new double*[nz];
  edenplot = new double*[nz];
  for(int i = 0; i < nz;i++)
  {
    i_b1Error[i] = new double[nx];
    i_b2Error[i] = new double[nx];
    i_bplot[i] = new double[nx];
    i_b_newplot[i] = new double[nz];
    edepplot[i] = new double[nx]{0.0};
    edenplot[i] = new double[nx]{0.0};
  }
  for(int i = 0; i < nz;i++)
  {
    for(int j = 0; j < nx;j++)
    {
      i_b1Error[j][i] = abs(icopy[0][j][i]-i_b_new[0][j][i])/((i_b_new[0][j][i] == 0) + i_b_new[0][j][i]);
      i_b2Error[j][i] = abs(icopy[1][j][i]-i_b_new[1][j][i])/((i_b_new[1][j][i] == 0) + i_b_new[1][j][i]);

      edenplot[j][i] = eden[j][i]/ncrit;
      i_bplot[j][i] = 8.53e-10*sqrt(i_b[0][j][i]+i_b[1][j][i]+1.0e-10)*(1.053/3.0);
      i_b_newplot[j][i] = 8.53e-10*sqrt(fmax(1.0e-10,i_b_new[0][j][i])+fmax(1.0e-10,i_b_new[1][j][i]))*(1.053/3.0);
      for(int m = 0;m<nbeams;m++)
      {
        edepplot[j][i]+=edep[m][j][i];
      }
    }
  }
  //cout << nrays << "  " << ncrossings << endl;
}
