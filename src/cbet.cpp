#include "implSim.hpp"
#include "customMath.hpp"
#include <iomanip>
#include <cstring>
using namespace std;
//initializing arrays and constants for the CBET subroutines
void initializeArr()
{
  //allocate space for the intensity and CBET multiplier arrays
  i_b = new double**[nbeams];//Initial intensity
  i_b_prev = new double**[nbeams];//Updated intensity
  i_b_new = new double**[nbeams];//Stores previous iteration intensity, initially 0
  W_new = new double**[nbeams];
  W = new double**[nbeams];
  for(int i = 0; i < nbeams; i++)
  {
    W[i] = new double*[nx];
    W_new[i] = new double*[nx];
    i_b_prev[i] = new double*[nx];
    i_b_new[i] = new double*[nx];
    i_b[i] = new double*[nx];
    for(int j = 0; j < nx; j++)
    {
      W[i][j] = new double[nx];
      W_new[i][j] = new double[nx];
      i_b[i][j] = new double[nz];
      i_b_new[i][j] = new double[nz];
      i_b_prev[i][j] = new double[nz]{0.0};
    }
  }
  cs = 1e2*sqrt(ec*(Z*Te_eV+3.0*Ti_eV)/mi_kg);	// acoustic wave speed, approx. 4e7 cm/s in this example
  //Initialize CBET array values
  for(int m = 0; m < nbeams; m++)
  {
    #pragma omp parallel for num_threads(threads)
    for(int i = 0; i < nx;i++)
    {
      for(int j = 0; j < nz;j++)
      {
        i_b[m][i][j] = edep[m][i][j];
        i_b_new[m][i][j] = edep[m][i][j];
        W[m][i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);//
        W_new[m][i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);//
        u_flow[i][j] = machnum[i][j]*cs;//
      }

    }
  }
  //cout << scientific;
  for(int i = 0; i < nbeams;i++)
  {
    #pragma omp parallel for num_threads(threads)
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
//Calculates the update multipler for CBET using th
void calculateMult(int ix, int iz, int i, int j, int m, int* numrays, int** marker, int** markerCopy)
{
  double constant1 = (pow(estat,2.0))/(4*(1.0e3*me)*c*omega*kb*Te*(1+3*Ti/(Z*Te)));

  //Iterate over the rays present in the zone (stops once beyond lowest ray count)
  int n2limit = fmin(present[i][ix][iz],numrays[i+1]);
  for ( int n2 = 0; n2 < n2limit; n2++)
  {
    double ne = eden[ix][iz];
    double epsilon = 1.0-ne/ncrit;
    double kmag = (omega/c)*sqrt(epsilon);         // magnitude of wavevector
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

    //Variable intensity in electric field causes change in multipliers
   double efield1 = sqrt(8.*pi*1.0e7*i_b[i][ix][iz]/c);             // electric field of a ray of beam 1
   double efield2 = sqrt(8.*pi*1.0e7*i_b[i+1][ix][iz]/c);             // electric field of ray of beam 2

   double P = (pow((iaw),2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw),2)*pow(eta,2));         // From Russ's paper
   double gain1 = constant1*pow(efield2,2)*(ne/ncrit)*(1/iaw)*P;               //L^-1 from Russ's paper
   double gain2 = constant1*pow(efield1,2)*(ne/ncrit)*(1/iaw)*P;
            //L^-1 from Russ's paper

                    // new energy of crossing (PROBE) ray (beam 2)
    //updating arrays based upon the calculated magnitude changes
    //Exponential gain function, important to note
    if (dkmag[i+1][marker[i+1][n2]][markerCopy[i+1][n2]] >= 1.0*dx)
    {
      //beam 2 CBET multiplier
      #pragma omp atomic write
      W_new[i+1][ix][iz] =W[i+1][ix][iz]*exp(-1*W[i][ix][iz]*dkmag[i+1][marker[i+1][n2]][markerCopy[i+1][n2]]*gain2/sqrt(epsilon));
      //beam 1 CBET multiplier
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
  //perform calculation for each crossing of each ray calculated for the beams in launch rays
  for(int i = 0; i < nbeams-1;i++)
  {
    #pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings; m++)
      {
        //If no ray is present at crossing m
        if(boxes[i][j][m][0] == 0 || boxes[i][j][m][1] == 0)
        {
          break;
        }
        int ix = boxes[i][j][m][0] - 1;
        int iz = boxes[i][j][m][1] - 1;
        //if there is a ray intersection at (ix, iz) and (ix, iz) is a valid coordinate
        if(ix != -1 && iz != -1 && intersections[ix][iz] != 0)
        {
          //First have to find all of the rays present from both beams in the coordinate (To generalize for N beams, encapsulate in a for loop)
          double tempStore = maxDev;
          int numrays[nbeams]{0};
          int** marker = new int*[nbeams];
          int** markerCopy = new int*[nbeams];
          //find the number of beams
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
          //allocate space for marker arrays
          marker[0] = new int[numrays[0]];
          marker[1] = new int[numrays[1]];
          markerCopy[0] = new int[numrays[0]];
          markerCopy[1] = new int[numrays[1]];
          int cnt1 = 0;
          int cnt2 = 0;
          //iterate over marked at (ix, iz) to find all incident rays
          for(int q = 0; q < numstored; q++)
          {
            //beam 1 rays
            if(marked[ix*nz+iz][q*nbeams + i] != 0)
            {
              marker[0][cnt1] = marked[ix*nz+iz][q*nbeams + i] - 1;
              markerCopy[0][cnt1] = marked[ix*nz+iz][q*nbeams + i] - 1;
              cnt1++;
            }
            //beam 2 rays
            if(marked[ix*nz+iz][q*nbeams + i + 1] != 0)
            {
              marker[1][cnt2] = marked[ix*nz+iz][q*nbeams + i + 1] - 1;
              markerCopy[1][cnt2] = marked[ix*nz+iz][q*nbeams + i + 1] - 1;
              cnt2++;
            }
          }
          //locate ray j's index within the zone
          int ray1num = 0;
          for(int r = 0; r < numrays[i]; r++)
          {
            if(marker[i][r] == j)
            {
              ray1num = r;
              break;
            }
          }

          //find a ray from beam 2 that intersects at (ix, iz)
          for(int n = 0; n < numrays[i+1]; n++)
          {
            for(int q = 0; q < ncrossings;q++)
            {
              int ix2 = boxes[i+1][marker[i+1][n]][q][0] - 1;
              int iz2 = boxes[i+1][marker[i+1][n]][q][1] - 1;
              if(ix == ix2 && iz == iz2)//if q has a valid crossing in (ix, iz)
              {
                markerCopy[i+1][n] = q;
                break;
              }
            }
          }
          calculateMult(ix, iz, i, j, m, numrays, marker, markerCopy);//calculate CBET multipliers using ray parameters

          //calculate the additive change in CBET field to be applied
          double additive_change[2];
          additive_change[0] = (-1.0 * (1.0 - (W_new[i][ix][iz]/W[i][ix][iz])) * abs(i_b[i][ix][iz] - i_b_prev[i][ix][iz])*(i_b[i][ix][iz] != 0));//change to be applied to beam 1 rays
          additive_change[1] = (-1.0 * (1.0 - (W_new[i+1][ix][iz]/W[i+1][ix][iz])) *  abs(i_b[i+1][ix][iz] - i_b_prev[i+1][ix][iz]));//change to be applied to beam 2 rays
          int kill = 0;
          if(i_b_new[i][ix][iz] < abs(additive_change[0]) && additive_change[0] < 0 || i_b_new[i][ix][iz] <= 0)//if beam 1's intensity is/will be below zero, kill the rays in this zone and any zones that they will pass through in the future, only happens at high intensities
          {
            additive_change[1] = i_b_new[i][ix][iz] * (i_b_new[i][ix][iz] > 0);//if beam 1 still has energy, that will be transferred to beam 2 to conserve energy
            #pragma omp atomic write
            i_b_new[i][ix][iz] = 0;
            additive_change[0] = 0;
            kill = 1;//stores that beam 1 has been killed in this sector
          }else
          {
            #pragma omp atomic update
            i_b_new[i][ix][iz] += additive_change[0];//standard case
          }
          tempStore = fmax(abs(additive_change[0]/(i_b[i][ix][iz]+(i_b_new[i][ix][iz] <= 1e-10))), tempStore);//variable used to store convergence
          tempStore = fmax(abs(additive_change[1]/(i_b[i+1][ix][iz]+(i_b_new[i+1][ix][iz] <= 1e-10))), tempStore);//variable used to store convergence
          #pragma omp atomic update
          i_b_new[i+1][ix][iz] += additive_change[1];//apply update to beam 1


          double x_prev = x[ix];
          double z_prev = z[iz];
          //for every future crossing of ray j
          for(int l = m+1; l < ncrossings;l++)
          {
            int ix_next_i = boxes[i][j][l][0] - 1;
            int iz_next_i = boxes[i][j][l][1] - 1;
            //if l is a valid crossing
            if(ix_next_i != -1 && iz_next_i != -1)
            {
              double x_curr = x[ix_next_i];
              double z_curr = z[iz_next_i];
              //if l is not a previously updated crossing (within this loop)
              if(x_curr != x_prev || z_curr != z_prev)
              {
                double add = additive_change[0];//*present[i][ix][iz]/present[i][ix_next_i][iz_next_i]; // Change to be applied
                double comp = add/(i_b_new[i][ix_next_i][iz_next_i]+(i_b_new[i][ix_next_i][iz_next_i] < 1e-10));//Variable simply for the sake of clarity when updating tempstore
                tempStore = fmax(tempStore, comp);//store convergence variable
                #pragma omp atomic update
                i_b_new[i][ix_next_i][iz_next_i] += add + kill*(-1)*(i_b_new[i][ix_next_i][iz_next_i] + add);//apply update, set intensity to 0 if kill == 1
                #pragma omp atomic update
                i_b[i][ix_next_i][iz_next_i] += kill*(-1)*i_b[i][ix_next_i][iz_next_i];//set initial intensity to 0 if kill == 1`
              }
              x_prev = x_curr;
              z_prev = z_curr;
            }else
            {
              break;
            }
          }

          int n2 = fmin(ray1num, numrays[i+1] - 1);//either stores the current ray in beam 1, or the last ray in beam 2, whichever is older, make sure no segfault
          //for every future crossing of ray n2, update intensity
          for(int l = markerCopy[i+1][n2] + 1; l < ncrossings;l++)
          {
            //find next grid coordinates
            int ix_next_other = boxes[i+1][marker[i+1][n2]][l][0] - 1;
            int iz_next_other = boxes[i+1][marker[i+1][n2]][l][1] - 1;
            //if it is a valid crossing
            if(ix_next_other != -1 && iz_next_other != -1)
            {
              double x_curr = x[ix_next_other];
              double z_curr = z[iz_next_other];
              //if this coordinate has not been updated in this loop
              if(x_curr != x_prev || z_curr != z_prev)
              {

                double add = additive_change[1];//*present[i][ix][iz]/present[i+1][ix_next_other][iz_next_other]; //change to be applied
                double comp = add/(i_b_new[i+1][ix_next_other][iz_next_other]+(i_b_new[i+1][ix_next_other][iz_next_other] < 1e-10)); // convergence variable
                tempStore = fmax(tempStore, comp);
                #pragma omp atomic update
                i_b_new[i+1][ix_next_other][iz_next_other] += add;//update downstream intensity
              }
              //Set up variables for next iteration
              x_prev = x_curr;
              z_prev = z_curr;
            }else
            {
              break;
            }
          }
          #pragma omp atomic write
          maxDev = tempStore;
        }
      }
    }
  }
  //Print energy conservation value
  if(printCBETDiagnostics)
  {
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
  }
  //Update arrays for the next iteration
  for(int i = 0; i < nbeams; i++)
  {
    #pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nx; j++)
    {
      for(int m = 0; m < nz; m++)
      {
        W[i][j][m] = W_new[i][j][m];
        i_b_prev[i][j][m] = i_b[i][j][m];
        i_b[i][j][m] = i_b_new[i][j][m];
      }
    }
  }
}
void cbet()
{

  //Used to track ray updates for beam 1, can be switched in the CBET function
  orderplot1 = new double*[nrays];
  for(int i = 0; i < nrays; i++)
  {
    orderplot1[i] = new double[nrays]{0};
  }
  //Initialize the arrays:
  auto start = chrono::high_resolution_clock::now();
  maxDev = 0;
  initializeArr();

  auto inter = chrono::high_resolution_clock::now();
  double prev = maxDev;
  //Initial CBET Update
  update_CBETPar();
  int cnt = 0;
  //CBET Iteration Loop
  while(abs(maxDev) > converge && cnt < maxIterations)
  {
    prev = maxDev;
    maxDev = 0;
    update_CBETPar();
    cnt++;
  }
  cout << maxDev << endl;
  auto end = chrono::high_resolution_clock::now();
  if(printUpdates)
  {
    cout << "Number of iterations: " << cnt << endl;
  }
  if(printSpecificTimings)
  {
    cout << "CBET Initialization Time: " << chrono::duration_cast<chrono::milliseconds>(inter-start).count() << " ms" << endl;
    cout << "CBET Calculation Time: " << chrono::duration_cast<chrono::milliseconds>(end-inter).count() << " ms" << endl;
  }
}
