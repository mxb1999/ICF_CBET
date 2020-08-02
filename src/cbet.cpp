#include "include/implSim.h"
#include "include/customMath.h"
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
        //
        W[m][i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);
        W_new[m][i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);
        W_init[m][i][j] =   W[m][i][j];
        u_flow[i][j] = machnum[i][j]*cs;
      }

    }
  }
  cout << scientific;
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
void gain_CBETSeq()
{
  double constant1 = (pow(estat,2.0))/(4*(1.0e3*me)*c*omega*kb*Te*(1+3*Ti/(Z*Te)));
  cout << "Calculating CBET gains" << endl;
  initializeArr();
  //performs the CBET calculation for each crossing of each ray calculated for the beams
  for(int i = 0; i < nbeams-1;i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings;m++)
      {
        //marked array has dimensions nx nz nrays nbeams
        int ix = boxes[i][j][m][0] - 1;
        int iz = boxes[i][j][m][1] - 1;
          if(ix >= 0 && iz >= 0 && intersections[ix][iz] != 0)
          {
            //given marked stores the trajectory of a given ray, the following loop gets the rays that pass through ix iz
            int numrays[nbeams];
            int* marker[nbeams];
            int* markerCopy[nbeams];
            for(int p = 0; p < nbeams; p++)
            {
              numrays[p] = 0;
              for(int q = 0; q < nrays; q++)
              {
                if(marked[ix*nz+iz][q*nbeams + p] != 0)
                {
                  numrays[p]++;
                }
              }
              if(numrays[p] > 0)
              {
                marker[p] = new int[numrays[p]];
                markerCopy[p] = new int[numrays[p]];
              }else
              {
                marker[p] = new int[1]{0};
                markerCopy[p] = new int[1]{0};
              }
              for(int q = 0; q < nrays; q++)
              {
                if(marked[ix*nz+iz][q*nbeams + p] != 0)
                {
                  marker[p][q] = marked[ix*nz+iz][q*nbeams + p] - 1;
                  markerCopy[p][q] = marked[ix*nz+iz][q*nbeams + p] - 1;
                }
              }
            }


          for(int r = 0; r < numrays[i];r++)
          {
            if(marker[i][r] == j)
            {
              ray1num = r;
              break;
            }
          }
          //
          for(int l = i+1; l < nbeams; l++)
          {
            for(int n = 0; n < numrays[l]; n++)
            {
              for(int q = 0; q <ncrossings; q++)
              {
               int ix2 = boxes[l][marker[l][n]][q][0] - 1;
                int iz2 = boxes[l][marker[l][n]][q][1] - 1;
                if ( ix == ix2 && iz == iz2 )
                {
                  markerCopy[l][n] = q;
                //  break;
                }

              }

            }
            int n2limit = fmin(present[i][ix][iz],numrays[l]);
            for ( int n2 = 0; n2 < n2limit; n2++)
            {
              double ne = eden[ix][iz];
              double epsilon = 1.0-ne/ncrit;
              double kmag = (omega/c)*sqrt(epsilon);         // magnitude of wavevector

              double kx1 = kmag * dkx[i][j][m] / (dkmag[i][j][m] + 1.0e-10);
              double kx2 = kmag * dkx[l][marker[l][n2]][markerCopy[l][n2]] / (dkmag[l][marker[l][n2]][markerCopy[l][n2]] + 1.0e-10);
              double kz1 = kmag * dkz[i][j][m]/(dkmag[i][j][m]+1.0e-10);
              double kz2 = kmag * dkz[l][marker[l][n2]][markerCopy[l][n2]] / (dkmag[l][marker[l][n2]][markerCopy[l][n2]] + 1.0e-10);
              double kiaw = sqrt(pow(kx2-kx1,2.0)+pow(kz2-kz1,2.0));
            // magnitude of the difference between the two vectors

              double ws = kiaw*cs;            // acoustic frequency, cs is a constant
              double omega1= omega;  // laser frequency difference. To start, just zero.
              double omega2 = omega;
              double eta = ((omega2-omega1)-(kx2-kx1)*u_flow[ix][iz])/(ws+1.0e-10);

             double efield1 = sqrt(8.*pi*1.0e7*i_b[i][ix][iz]/c);             // initial electric field of rayi;lnlni46
            // double efield2 = sqrt(8.*pi*1.0e7*(*(i_b2[ix]+iz))/c);             // initial electric field of ray

             double P = (pow((iaw),2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw),2)*pow(eta,2));         // From Russ's paper
             //double gain1 = constant1*pow(efield2,2)*(ne/ncrit)*(1/iaw)*P;               //L^-1 from Russ's paper
             double gain2 = constant1*pow(efield1,2)*(ne/ncrit)*(1/iaw)*P;
                      //L^-1 from Russ's paper

                              // new energy of crossing (PROBE) ray (beam 2)
              //updating arrays based upon the calculated magnitude changes
              if (dkmag[l][marker[l][n2]][markerCopy[l][n2]] >= 1.0*dx)
              {
                //#pragma omp atomic write
                W_storage[l][j][ix][iz] =W[l][ix][iz]*exp(-1*W[i][ix][iz]*dkmag[l][marker[l][n2]][markerCopy[l][n2]]*gain2/sqrt(epsilon));
                //#pragma omp atomic write
                W_storage[i][j][ix][iz] =W[i][ix][iz]*exp(1*W[l][ix][iz]*dkmag[i][j][m]*gain2/sqrt(epsilon));
              }
            }
          }
          for(int p = 0; p < nbeams;p++)
          {
            delete [] marker[p];
            delete [] markerCopy[p];
          }
        }
      }
      if ( j % 20 == 0 )
      {
      //  cout << "     ..."<<(int)(100.*(1.0-(double)(j/double(1*nrays))))<<"%  remaining..."<<endl;
      }
    }
  }

  for(int q = 0; q < nbeams; q++)
  {
    //#pragma omp parallel for num_threads(threads)
    for(int i = 0; i < nx; i++)
    {
      for(int j = 0; j < nz; j++)
      {
        for(int m = nrays - 1; m >= 0; m--)
        {
          if(W_storage[q][m][i][j] != 0)
          {
            W_new[q][i][j] = W_storage[q][m][i][j];
            break;
          }
        }
      }
    }
  }

}

//updating CBET Intensities: Original sequential version
void update_CBETSeq()
{
  //storage array to allow for data level parallelization
  cout << "Updating CBET Intensities" << endl;

  //perform calculation for each crossing of each ray calculated for the beams in launch rays
  for(int i = 0; i < nbeams-1;i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings; m++)
      {
        if(boxes[i][j][m][0] == 0 || boxes[i][j][m][1] == 0)
        {
          break;
        }
        int ix = boxes[i][j][m][0] - 1;
        int iz = boxes[i][j][m][1] - 1;
        //if two beams intersect

        if(intersections[ix][iz] != 0)
        {

          double tempStore = maxDev;
          //find nonzero locations of rays
          int numrays[nbeams];
          int** marker = new int*[nbeams];
          int** markerCopy = new int*[nbeams];
          for(int p = 0; p < nbeams; p++)
          {
            numrays[p] = 0;
            for(int q = 0; q < nrays; q++)
            {
              if(marked[ix*nz+iz][q*nbeams + p] != 0)
              {
                numrays[p]++;
              }
            }
            marker[p] = new int[numrays[p]];
            markerCopy[p] = new int[numrays[p]];
            for(int q = 0; q < nrays; q++)
            {
              if(marked[ix*nz+iz][q*nbeams + p] != 0)
              {
                marker[p][q] = marked[ix*nz+iz][q*nbeams + p] - 1;
                markerCopy[p][q] = marked[ix*nz+iz][q*nbeams + p] - 1;
              }
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
          for(int k = i+1; k < nbeams; k++)
          {
            if(k != i)
            {
              for(int n = 0; n < numrays[k]; n++)
              {
                for(int q = 0; q < ncrossings;q++)
                {
                  int ix2 = boxes[k][marker[k][n]][q][0] - 1;
                  int iz2 = boxes[k][marker[k][n]][q][1] - 1;
                  if(ix == ix2 && iz == iz2)
                  {
                    markerCopy[k][n] = q;
                    break;
                  }
                }
              }
            }
          }

          //calculate the additive change in CBET field to be applied
          cout << scientific;
          double additive_change[nbeams];
          for(int q = 0; q < nbeams; q++)
          {
            additive_change[q] = (-1.0 * (1.0 - (W_new[q][ix][iz]/W_init[q][ix][iz])) * i_b[q][ix][iz]);
            while(abs(additive_change[q]/i_b[q][ix][iz]) > maxInc)
            {
              additive_change[q] *= maxInc;
            }
            tempStore = fmax(additive_change[q]/i_b_new[q][ix][iz], tempStore);
            i_b_new[q][ix][iz] += additive_change[q];

            int iter;
            int n2 = fmin(ray1num, numrays[q] - 1);
            if(q == i)
            {
              iter = m+1;
            }else
            {
              iter = markerCopy[q][n2] + 1;
            }
            double x_prev = x[ix];
            double z_prev = z[iz];
            for(int l = iter; l < ncrossings;l++)
            {
              int ix_next;
              int iz_next;
              if(i == q)
              {
                ix_next = boxes[q][j][l][0] - 1;
                iz_next = boxes[q][j][l][1] - 1;
              }else
              {
                ix_next = boxes[q][marker[q][n2]][l][0] - 1;
                iz_next = boxes[q][marker[q][n2]][l][1] - 1;
              }
              if(ix_next == -1 || iz_next == -1)
              {
                break;
              }else
              {
                double x_curr = x[ix_next];
                double z_curr = z[iz_next];
                if(x_curr != x_prev || z_curr != z_prev)
                {
                  int loc;
                  if(i == q)
                  {
                    loc = i;
                  }else
                  {
                    loc = q;
                  }
                i_b_new[q][ix_next][iz_next] += additive_change[q] * ((double)present[i][ix][iz])/present[loc][ix_next][iz_next];
                }
                x_prev = x_curr;
                z_prev = z_curr;
              }
            }
          }
          maxDev = tempStore;

        }
      }
      if ( j % 20 == 0 )
      {
        //cout << "     ..."<<(int)(100.*(1.0-(double)(j/double(1*nrays))))<<"%  remaining..."<<endl;
      }

    }
  }
  cout << "Finished Update loop" << endl;

}



















double update_Change(double additive_change, double additive_change_old, int iteration, int ix, int iz)
{
  double temp = additive_change;
  additive_change = 1 - additive_change;
  if(additive_change <= 0)
  {
    additive_change = 1e-5 + (additive_change - 1e-5)*(additive_change == 0);
    if(iteration > 10)
    {
      cout << "Setting frac_change equal to min" << endl;
    }
  }
  double stepSizeRatio = 1; //In Russ' paper, ratio of grid ds to median ds, for rays not
  //changePerStep = (1-additive_change) /
  return 0.0;
}
void gain_CBETPar()
{
  double constant1 = (pow(estat,2.0))/(4*(1.0e3*me)*c*omega*kb*Te*(1+3*Ti/(Z*Te)));
  cout << "Calculating CBET gains" << endl;

  //performs the CBET calculation for each crossing of each ray calculated for the beams
  for(int i = 0; i < nbeams-1;i++)
  {
    //#pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings;m++)
      {
        //marked array has dimensions nx nz nrays nbeams
        int ix = boxes[i][j][m][0] - 1;
        int iz = boxes[i][j][m][1] - 1;

          if(ix >= 0 && iz >= 0 && intersections[ix][iz] != 0)
          {
            int numrays[nbeams];
            int* marker[nbeams];
            int* markerCopy[nbeams];
            for(int p = 0; p < nbeams; p++)
            {
              numrays[p] = 0;
              for(int q = 0; q < nrays; q++)
              {
                if(marked[ix*nz+iz][q*nbeams + p] != 0)
                {
                  numrays[p]++;
                }
              }
              if(numrays[p] > 0)
              {
                marker[p] = new int[numrays[p]];
                markerCopy[p] = new int[numrays[p]];
              }else
              {
                marker[p] = new int[1]{0};
                markerCopy[p] = new int[1]{0};
              }
              for(int q = 0; q < nrays; q++)
              {
                if(marked[ix*nz+iz][q*nbeams + p] != 0)
                {
                  marker[p][q] = marked[ix*nz+iz][q*nbeams + p] - 1;
                  markerCopy[p][q] = marked[ix*nz+iz][q*nbeams + p] - 1;
                }
              }
            }


          for(int r = 0; r < numrays[i];r++)
          {
            if(marker[i][r] == j)
            {
              ray1num = r;
              break;
            }
          }
          for(int l = i+1; l < nbeams; l++)
          {
            for(int n = 0; n < numrays[l]; n++)
            {
              for(int q = 0; q <ncrossings; q++)
              {
               int ix2 = boxes[l][marker[l][n]][q][0] - 1;
                int iz2 = boxes[l][marker[l][n]][q][1] - 1;
                if ( ix == ix2 && iz == iz2 )
                {
                  markerCopy[l][n] = q;
                  break;
                }

              }

            }
            int n2limit = fmin(present[i][ix][iz],numrays[l]);
            for ( int n2 = 0; n2 < n2limit; n2++)
            {
              double ne = eden[ix][iz];
              double epsilon = 1.0-ne/ncrit;
              double kmag = (omega/c)*sqrt(epsilon);         // magnitude of wavevector
              double l_d = sqrt(e0*kb*Te/(pow(ec,2)*ne)); //Debye length of the plasma
              double N_d = ne*4*pi/3*pow(l_d,3); //number of electrons in the debye sphere
              double coulLog = log(9*N_d/Z); //Coulomb logarithm of the plasma
              double kx1 = kmag * dkx[i][j][m] / (dkmag[i][j][m] + 1.0e-10);
              double kx2 = kmag * dkx[l][marker[l][n2]][markerCopy[l][n2]] / (dkmag[l][marker[l][n2]][markerCopy[l][n2]] + 1.0e-10);
              double kz1 = kmag * dkz[i][j][m]/(dkmag[i][j][m]+1.0e-10);
              double kz2 = kmag * dkz[l][marker[l][n2]][markerCopy[l][n2]] / (dkmag[l][marker[l][n2]][markerCopy[l][n2]] + 1.0e-10);
              double kiaw = sqrt(pow(kx2-kx1,2.0)+pow(kz2-kz1,2.0));
            // magnitude of the difference between the two vectors
            //  double vei = 4*sqrt(2*pi)*pow(ec,4)*pow(Z,2),
              double ws = kiaw*cs;            // acoustic frequency, cs is a constant
              double omega1= omega;  // laser frequency difference. To start, just zero.
              double omega2 = omega;
              double eta = ((omega2-omega1)-(kx2-kx1)*u_flow[ix][iz])/(ws+1.0e-10);


             double efield1 = sqrt(8.*pi*1.0e7*i_b[i][ix][iz]/c);             // initial electric field of rayi
             //double efield2 = sqrt(8.*pi*1.0e7*i_b[l][ix][iz]/c);             // initial electric field of ray

             double P = (pow((iaw),2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw),2)*pow(eta,2));         // From Russ's paper
             //double gain1 = constant1*pow(efield2,2)*(ne/ncrit)*(1/iaw)*P;               //L^-1 from Russ's paper
             double gain2 = constant1*pow(efield1,2)*(ne/ncrit)*(1/iaw)*P;
                      //L^-1 from Russ's paper
             double L_a = 0;
             //cout << L_a << endl;
                              // new energy of crossing (PROBE) ray (beam 2)
              //updating arrays based upon the calculated magnitude changes
              if (dkmag[l][marker[l][n2]][markerCopy[l][n2]] >= 1.0*dx)
              {
                #pragma omp atomic write
                W_storage[l][j][ix][iz] =W[l][ix][iz]*exp(-1*W[i][ix][iz]*dkmag[l][marker[l][n2]][markerCopy[l][n2]]*gain2/sqrt(epsilon));
                #pragma omp atomic write
                W_storage[i][j][ix][iz] =W[i][ix][iz]*exp(1*W[l][ix][iz]*dkmag[i][j][m]*gain2/sqrt(epsilon));
              }
            }
          }
          for(int p = 0; p < nbeams;p++)
          {
            delete [] marker[p];
            delete [] markerCopy[p];
          }
        }
      }

      if ( j % 20 == 0 )
      {
      //  cout << "     ..."<<(int)(100.*(1.0-(double)(j/double(1*nrays))))<<"%  remaining..."<<endl;
      }
    }
  }

  for(int q = 0; q < nbeams; q++)
  {
    #pragma omp parallel for num_threads(threads)
    for(int i = 0; i < nx; i++)
    {
      for(int j = 0; j < nz; j++)
      {
        for(int m = nrays - 1; m >= 0; m--)
        {
          if(W_storage[q][m][i][j] != 0)
          {
            W_new[q][i][j] = W_storage[q][m][i][j];
            break;
          }
        }
      }
    }
  }
}

//updating CBET Intensities
void update_CBETPar()
{
  //storage array to allow for data level parallelization
  cout << "Updating CBET Intensities" << endl;
  int checker = 0;

  //perform calculation for each crossing of each ray calculated for the beams in launch rays
  for(int i = 0; i < nbeams-1;i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings; m++)
      {
        if(boxes[i][j][m][0] == 0 || boxes[i][j][m][1] == 0)
        {
          break;
        }
        int ix = boxes[i][j][m][0] - 1;
        int iz = boxes[i][j][m][1] - 1;
        //if two beams intersect
        if(intersections[ix][iz] != 0)
        {

          double tempStore = maxDev;
          //find nonzero locations of rays
          int numrays[nbeams]{0};
          int** marker = new int*[nbeams];
          int** markerCopy = new int*[nbeams];
          for(int q = 0; q < nrays; q++)
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
          for(int q = 0; q < nrays; q++)
          {
            if(marked[ix*nz+iz][q*nbeams + i + 1] != 0)
            {
              marker[1][q] = marked[ix*nz+iz][q*nbeams + i + 1] - 1;
              markerCopy[1][q] = marked[ix*nz+iz][q*nbeams + i + 1] - 1;
            }
            if(marked[ix*nz+iz][q*nbeams + i] != 0)
            {
              marker[0][q] = marked[ix*nz+iz][q*nbeams + i] - 1;
              markerCopy[0][q] = marked[ix*nz+iz][q*nbeams + i] - 1;
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
          for(int n = 0; n < numrays[i+1]; n++)
          {
            for(int q = 0; q < ncrossings;q++)
            {
              int ix2 = boxes[i+1][marker[i+1][n]][q][0] - 1;
              int iz2 = boxes[i+1][marker[i+1][n]][q][1] - 1;
              if(ix == ix2 && iz == iz2)
              {
                markerCopy[i+1][n] = q;
              }
            }
          }


          //calculate the additive change in CBET field to be applied
          cout << scientific;
          double additive_change[2];
          additive_change[0] = (-1.0 * (1.0 - (W_new[i][ix][iz]/W_init[i][ix][iz])) * abs(i_b[i][ix][iz] - i_b_prev[i][ix][iz]));
          additive_change[1] = (-1.0 * (1.0 - (W_new[i+1][ix][iz]/W_init[i+1][ix][iz])) *  abs(i_b[i+1][ix][iz] - i_b_prev[i+1][ix][iz]));
          tempStore = fmax(abs(additive_change[0]/i_b_new[i][ix][iz]), tempStore);
          if(ix == 101 && iz == 100)
          {
            convergeplot[count] = checker*i_b_new[i][ix][iz];
            count++;
            checker = 1;
          }
          i_b_new[i][ix][iz] += additive_change[0];
          cout << additive_change[0] << endl;
          i_b_new[i+1][ix][iz] += additive_change[1];
          //cout << "Initial lower for beam " << i << " and Ray " << j << " at " << ix << " " << iz << endl;
          //cout << "Initial raise for beam " << i+1  << " at " << ix << " " << iz << endl;
          counter[ix][iz]++;
          int iter;
          int n2 = fmin(ray1num, numrays[i] - 1);
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
                if(ix_next_i == 101 && iz_next_i == 100)
                {
                  convergeplot[count] = checker*i_b_new[i][ix_next_i][iz_next_i] ;
                  count++;
                  checker = 1;
                }
                counter[ix_next_i][iz_next_i]++;
                tempStore = fmax(tempStore, abs(additive_change[0] * ((double)present[i][ix][iz])/present[i][ix_next_i][iz_next_i]/i_b_new[i][ix_next_i][iz_next_i]));
                i_b_new[i][ix_next_i][iz_next_i] += additive_change[0] * (double)present[i][ix][iz]/present[i][ix_next_i][iz_next_i];
              //  cout << "Downstream lower for beam " << i << " and Ray " << j << " at " << ix_next_i << " " << iz_next_i << endl;
              }
              x_prev = x_curr;
              z_prev = z_curr;
            }
          }
          for(int l = markerCopy[i+1][n2] + 1; l < ncrossings;l++)
          {
            int ix_next_other = boxes[i+1][markerCopy[i+1][n2]][l][0] - 1;

            int iz_next_other = boxes[i+1][markerCopy[i+1][n2]][l][1] - 1;
            if(ix_next_other != -1 && iz_next_other != -1)
            {
              double x_curr = x[ix_next_other];
              double z_curr = z[iz_next_other];
              if(x_curr != x_prev || z_curr != z_prev)
              {

                if(ix_next_other == 101 && iz_next_other == 100)
                {
                  convergeplot[count] = checker*i_b_new[i+1][ix_next_other][iz_next_other];
                  count++;
                  checker = 1;
                }
                counter[ix_next_other][iz_next_other]++;
                tempStore = fmax(tempStore, abs(additive_change[1] * ((double)present[i][ix][iz])/present[i+1][ix_next_other][iz_next_other]/i_b_new[i+1][ix_next_other][iz_next_other]));
                i_b_new[i+1][ix_next_other][iz_next_other] += additive_change[1] * (double)present[i][ix][iz]/present[i+1][ix_next_other][iz_next_other]/i_b_new[i+1][ix_next_other][iz_next_other];
                //cout << "Downstream raise for beam " << i+1 << " and Ray " << markerCopy[i+1][n2] << " at " << ix_next_other << " " << iz_next_other << endl;
              }
              x_prev = x_curr;
              z_prev = z_curr;
            }
          }

          #pragma omp atomic write
          maxDev = tempStore;

        }
      }
      if ( j % 20 == 0 )
      {
        //cout << "     ..."<<(int)(100.*(1.0-(double)(j/double(1*nrays))))<<"%  remaining..."<<endl;
      }

    }
  }
  cout << "Finished Update loop" << endl;
for(int i = 0; i < nbeams; i++)
{
  for(int j = 0; j < nx; j++)
  {
    for(int m = 0; m < nz; m++)
    {
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
  convergeplot = new double[20016]{0};
  counter = new int*[nx];
  for(int i = 0; i < nx; i++)
  {
    counter[i] = new int[nz]{0};
  }
  orderplot1 = new double*[nrays];
  orderplot2 = new double*[nrays];
  for(int i = 0; i < nrays; i++)
  {
    orderplot1[i] = new double[nrays]{0};
    orderplot2[i] = new double[nrays]{0};
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
  initializeArr();
  gain_CBETSeq();
  update_CBETSeq();
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {

      icopy[0][i][j] = i_b_new[0][i][j];
      icopy[1][i][j] = i_b_new[1][i][j];

    }
  }
  cout << "nRays = " << nrays << endl;
  maxDev = 0;
  maxInc = maxIncrement;
  double prev = maxDev;
  initializeArr();
  auto start = chrono::high_resolution_clock::now();
  gain_CBETPar();
  auto inter = chrono::high_resolution_clock::now();
  update_CBETPar();
  cout << "MaxDev: " << maxDev << endl;
  int cnt = 0;

  while(abs(maxDev) > converge && cnt < 50)
  {
  //maxInc *= exp(-1*pow(cnt,2));
    iter++;
    prev = maxDev;
    maxDev = 0;
    cout << "_________________Loop________________" << endl;
    gain_CBETPar();
    update_CBETPar();
    cout << "Count: " << count << endl;
    cout << "MaxDev: " << maxDev << endl;
    cout << "MaxInc: " << maxInc << endl;
    cnt++;
  }

  cout << "Reeee" << endl;
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {
      if(cntMax < counter[i][j])
      {
        cntMax = counter[i][j];
      //  cout << i << " || " << j << endl;
      }
      if(i_b_new[0][i][j] !=icopy[0][i][j])
      {
      cout << i_b_new[0][i][j] << " " << icopy[0][i][j]<< " " << edep[0][i][j] << " " << i << " " << j << endl;
      }
      if(i_b_new[1][i][j] != 0)
      {
        //cout << i_b_new[1][i][j] << endl;
      }

    }
  }
  cout << cntMax << endl;
  auto parallel = chrono::high_resolution_clock::now();
  i_b1Error = new double*[nx];
  i_b2Error = new double*[nx];

  auto seqgain = chrono::high_resolution_clock::now();

  //gain_CBETSeq();
  auto sequpdate = chrono::high_resolution_clock::now();
//  update_CBETSeq();
  auto finish = chrono::high_resolution_clock::now();
  cout << "Gain Parallel CBET: " << chrono::duration_cast<chrono::milliseconds>(inter-start).count() << " ms" << endl;
  cout << "Update Parallel CBET: " << chrono::duration_cast<chrono::milliseconds>(parallel-inter).count() << " ms" << endl;
  cout << endl;
  cout << "Gain Sequential CBET: " << chrono::duration_cast<chrono::milliseconds>(sequpdate-seqgain).count() << " ms" << endl;
  cout << "Update Sequential CBET: " << chrono::duration_cast<chrono::milliseconds>(finish-sequpdate).count() << " ms" << endl;
  i_bplot = new double*[nx];
  i_b_newplot = new double*[nx];
  edepplot = new double*[nx];
  edenplot = new double*[nx];
  for(int i = 0; i < nx;i++)
  {
    i_b1Error[i] = new double[nz];
    i_b2Error[i] = new double[nz];
    i_bplot[i] = new double[nz];
    i_b_newplot[i] = new double[nz];
    edepplot[i] = new double[nz]{0.0};
    edenplot[i] = new double[nz]{0.0};
    for(int j = 0; j < nz;j++)
    {
      if(icopy[0][i][j] != 0)
      {
        i_b1Error[i][j] = icopy[0][i][j]-i_b_new[0][i][j];
      }
      if(icopy[1][i][j] != 0)
      {
        i_b2Error[i][j] = icopy[1][i][j]-i_b_new[1][i][j];
      }
      edenplot[i][j] = eden[i][j]/ncrit;
      i_bplot[i][j] = 8.53e-10*sqrt(i_b[0][i][j]+i_b[1][i][j]+1.0e-10)*(1.053/3.0);
      i_b_newplot[i][j] = 8.53e-10*sqrt(fmax(1.0e-10,i_b_new[0][i][j])+fmax(1.0e-10,i_b_new[1][i][j]))*(1.053/3.0);


      for(int m = 0;m<nbeams;m++)
      {
        edepplot[i][j]+=edep[m][i][j];
      }
    }
  }
  cout << nrays << "  " << ncrossings << endl;
}
