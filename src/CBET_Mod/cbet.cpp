#include "CBET_Interface.hpp"
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
        i_b[m][i][j] = edep[m][i][j];///present[m][i][j];
        i_b_new[m][i][j] = edep[m][i][j];///present[m][i][j];
        W[m][i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);//
        W_new[m][i][j] = sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);//
        u_flow[i][j] = machnum[i][j]*cs;//
      }

    }
  }
  //cout << scientific;
  if(!switchvar)
  {
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
/*  for(int i = 0; i < ncrossings;i++)
  {
    if(dkmag[0][0][i] >= 1e-10)
    {
      cout << dkmag[0][0][i] << "\n";

    }
  }
*/
}

void updateMult()
{
  //Outer summation over "sheets"/beams
    //inner summation over interactions with other ray k at location l in that beam
    //Sum G_jmkl from Russ' paper, function of intensity and validity of interaction
      //Multiply G_jmkl by product of energy change ratio for ray k along path so far
  double constant1 = (pow(estat,2.0))/(4*(1.0e3*me)*c*omega*kb*Te*(1+3*Ti/(Z*Te)));

  for(int i = 0;i < nbeams-1;i++)
  {
    
    for(int j = 0; j < nrays;j++)//for each ray
    {
      for(int m = 0; m < ncrossings;m++)
      {
        //identify all rays that j interacts with
        int boxx = *vec4D(boxes, i,j,m, 0, nrays, ncrossings,2);
        int boxz = *vec4D(boxes, i,j,m, 1, nrays, ncrossings,2);
        if(!boxx || !boxz)//no more valid locations
        {
          break;
        }
        boxx--;
        boxz--;
        if(intersections[boxx][boxz])//if there is a valid intersection
        {
          double ne = eden[boxx][boxz];
          double epsilon = 1.0-ne/ncrit;
          double kmag = (omega/c)*sqrt(epsilon);         // magnitude of wavevector
          int rayloc = 0;
          double sum1 = 0;//outer sum, over each beam
          double s_jm = dkmag[i][j][m];
          double L = dx/100*ncrit/eden[boxx][boxz];
          for(int n = 0; n < nbeams;n++)//for every higher order beam (only consider those with higher indices to prevent counting twice)
          {
            if(n == i)
            {
              if(i == nbeams-1)
                break;
              n++;
            }
            double sum2 = 0;//inner sum, over each interaction
            int cnt = 0;
            int icnt = 0;
            for(int l = 0; l < numstored;l++)
            {
              if(*vec4D(marked, i,boxx,boxz, l, nrays, nx,nz))
              {
                icnt++;
              }
              if(*vec4D(marked, n,boxx,boxz, l, nrays, nx,nz))
              {
                cnt++;
              }
            }
          //printf("Count: %d %d :: %d %d\n", icnt, cnt, present[i][boxx][boxz],present[n][boxx][boxz]);
          for(int l = 0; l < cnt;l++)//for every crossing 
          {
            
            int initx = *vec4D(boxes, n,j,0, 0, nrays, ncrossings,2)-1;
            int initz = *vec4D(boxes, n,j,0, 1, nrays, ncrossings,2)-1;
            int kcross = 0;
            double prod = 1;//inner product, each change in energy for kth beam
            double xi = 0.9*(omega*L/c);
            double dS_ratio = (double)present[n][initx][initz]/present[n][boxx][boxz];
            //printf("Ratio: %e\n", dS_ratio);
            double delS = fmin(1/sqrt(epsilon)*dS_ratio, xi*sqrt(ne/ncrit));
            for(int q = 1; q < ncrossings;q++)
            {
              int kboxx = *vec4D(boxes, n,l,q, 0, nrays, ncrossings,2);//location of kth beam
              int kboxz = *vec4D(boxes, n,l,q, 1, nrays, ncrossings,2);
              
              if(!kboxx || !kboxz)//if invalid crossing
              {
                break;
              }
            
              kboxx--;
              kboxz--;
              
              if(kboxx == boxx && kboxz == boxz)//if crossing q is at this grid location
              {
                kcross = q;
                break;
              }
              int kboxx_next = *vec4D(boxes, i,l,q+1, 0, nrays, ncrossings,2);//next crossing
              int kboxz_next = *vec4D(boxes, i,l,q+1, 1, nrays, ncrossings,2);
              if(!kboxx_next || !kboxz_next)
              {
                break;
              }
              kboxx_next--;
              kboxz_next--;
              if(!present[n][kboxx][kboxz])
              {
                //printf("None Present: %d %d :: ([%e,%e], [%e,%e])\n", kboxx, kboxz, kboxx*dx + xmin,(kboxx+1)*dx + xmin,kboxz*dz + xmin,(kboxz+1)*dz + xmin);
              }
              double this_E = W[n][kboxx][kboxz]/present[n][kboxx][kboxz];
              double next_E = W[n][kboxx_next][kboxz_next]/present[n][kboxx_next][kboxz_next];
              double ratio = next_E/this_E;
              prod *= ratio;
              
            }
            double epsilon = 1.0-ne/ncrit;
            double kmag = (omega/c)*sqrt(epsilon);         // magnitude of wavevector
            double kx1 = kmag * dkx[i][j][m] / (dkmag[i][j][m] + 1.0e-10);
            double kx2 = kmag * dkx[n][l][kcross] / (dkmag[i+1][l][kcross] + 1.0e-10);
            double kz1 = kmag * dkz[i][j][m]/(dkmag[i][j][m]+1.0e-10);
            double kz2 = kmag * dkz[n][l][kcross] / (dkmag[i+1][l][kcross] + 1.0e-10);
            double kiaw = sqrt(pow(kx2-kx1,2.0)+pow(kz2-kz1,2.0));
          // magnitude of the difference between the two vectors
          //  double vei = 4*sqrt(2*pi)*pow(ec,4)*pow(Z,2),
            double ws = kiaw*cs;            // acoustic frequency, cs is a constant
            double omega1= omega;  // laser frequency difference. To start, just zero.
            double omega2 = omega;
            double eta = ((omega2-omega1)-(kx2-kx1)*u_flow[boxx][boxz])/(ws+1.0e-10);
            double efield2 = sqrt(8.*pi*1.0e7*i_b[n][boxx][boxz]/c);             // electric field of a ray of beam 1
            double P = (pow((iaw),2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw),2)*pow(eta,2));         // From Russ's paper
            double L_jmkl = constant1*pow(efield2,2)*ne/ncrit*(1/iaw)*P;
            double G = dkmag[i][j][m]*delS/(sqrt(epsilon)*L_jmkl*cnt);
            sum2+=(G*prod);
            mult[(n*nrays+l)*ncrossings+kcross] = exp(sum2);
            
          }
          short check = (sum1 != sum1);

          sum1+=sum2;

          
          //printf("sum2: %e\n", sum2);
        }
        double La = c*sqrt(epsilon)*ncrit/(iaw*omega*ne);
        double A_jm = s_jm*La;
        mult[(i*nrays+j)*ncrossings+m] = exp(sum1-A_jm);
        
       // printf("Ray %d at Crossing %d: %e\n", j,m, mult[(i*nrays+j)*ncrossings+m]);
       // printf("\t|%e|%e|%e|%e|\n", sum1, s_jm, A_jm, La);
      }
        
      }
    }
  }
}
void updateRays()
{
  double convergeStore[threads]{0};
  for(int i = 0; i < nbeams;i++)
  {
    for(int j = 0; j < nrays;j++)
    {
      int numT = omp_get_thread_num();
      int initx = *vec4D(boxes, i,j,0,0,nrays,ncrossings,2);
      int initz = *vec4D(boxes, i,j,0,1,nrays,ncrossings,2);
      double prevE = sqrt(8*pi*i_b[i][initx-1][initz-1]/c);
      for(int m = 1; m < ncrossings;m++)
      {
        int currx = *vec4D(boxes, i,j,m,0,nrays,ncrossings,2);
        int currz = *vec4D(boxes, i,j,m,1,nrays,ncrossings,2);
        if(!currx || !currz)
        {
          break;
        }
        currx--;
        currz--;
        double currE = mult[(i*nrays+j)*ncrossings + m]*prevE;
       
       // printf("%d at %d: Curr E: %e, Prev E: %e\n",j,m, currE, prevE);
       // printf("\tNew Intensity: %e %e\n", pow(currE, 2)/(8*pi)*c, mult[(i*nrays+j)*ncrossings + m]);
        i_b_new[i][currx][currz] = pow(currE, 2)/(8*pi)*c;
        double temp = convergeStore[numT];
        double comp = sqrt(8*pi*i_b[i][currx][currz]/c);
        double denominator = (abs(comp) < 1e-10) ? fmax(comp,1) : comp;
        convergeStore[numT] = fmax(currE/denominator-1, convergeStore[numT]);
        printf("%e %e\n", currE, sqrt(8*pi*i_b[i][currx][currz]/c));
        if(temp != convergeStore[numT])
        {
          printf("Converge: %e %e %d %d %d %e %e\n",currE, sqrt(8*pi*i_b[i][currx][currz]/c), i, j, m, mult[(i*nrays+j)*ncrossings + m], convergeStore[numT]);
        }
        prevE=currE;
      }
    }
  }
  for(int i = 0; i < threads; i++)
  {
    maxDev = fmax(convergeStore[i], maxDev);
  }
  
  for(int i = 0; i < nbeams; i++)
  {
    for(int j = 0; j < nx; j++)
    {
      for(int m = 0; m < nz; m++)
      {
        i_b[i][j][m] = i_b_new[i][j][m];
      }
    }
  }
  
  for(int i = 0; i < nbeams; i++)
  {
    for(int j = 0; j < nx; j++)
    {
      for(int m = 0; m < nz; m++)
      {
       // if(i_b[i][j][m] != i_b_new[i][j][m])
        //{
          //rintf("What the fuck\n");
       // }
      }
    }
  }
}
//Calculates the update multipler for CBET using th
void calculateMult(int ix, int iz, int i, int j, int m, int* numrays, int* marker, int* markerCopy)
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
    double kx2 = kmag * dkx[i+1][marker[(i+1)*ncrossings+n2]][markerCopy[(i+1)*ncrossings+n2]] / (dkmag[i+1][marker[(i+1)*ncrossings+n2]][markerCopy[(i+1)*ncrossings+n2]] + 1.0e-10);
    double kz1 = kmag * dkz[i][j][m]/(dkmag[i][j][m]+1.0e-10);
    double kz2 = kmag * dkz[i+1][marker[(i+1)*ncrossings+n2]][markerCopy[(i+1)*ncrossings+n2]] / (dkmag[i+1][marker[(i+1)*ncrossings+n2]][markerCopy[(i+1)*ncrossings+n2]] + 1.0e-10);
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
   double s_jm = dkmag[i+1][marker[(i+1)*ncrossings+n2]][markerCopy[(i+1)*ncrossings+n2]];

   double P = (pow((iaw),2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw),2)*pow(eta,2));         // From Russ's paper
   double gain1 = constant1*pow(efield2,2)*(ne/ncrit)*(1/iaw)*P;               //L^-1 from Russ's paper
   double sum1 = 0;
   //double gain2 = constant1*pow(efield1,2)*(ne/ncrit)*(1/iaw)*P;
            //L^-1 from Russ's paper

                    // new energy of crossing (PROBE) ray (beam 2)
    //updating arrays based upon the calculated magnitude changes
    //Exponential gain function, important to note

    //if (dkmag[i+1][marker[(i+1)*ncrossings+n2]][markerCopy[(i+1)*ncrossings+n2]] )//  <= sqrt(pow(dx,2) + pow(dz,2)))
   // {
      //if(j == 0)
      //{
              //#pragma omp atomic write

        gain1arr[ix][iz] += gain1;
        //gain2arr[ix][iz] += gain2;
        mag[ix][iz]+=s_jm;
       // printf("Crossing: %d w/ mult %e, orig %e, mag %e, gain %e\n",m,exp(-1*W[i][ix][iz]*dkmag[i+1][marker[(i+1)*ncrossings+n2]][markerCopy[(i+1)*ncrossings+n2]]*gain2/sqrt(epsilon)),
        //  W[i][ix][iz], dkmag[i+1][marker[(i+1)*ncrossings+n2]][markerCopy[(i+1)*ncrossings+n2]], gain2);
     // }
      //beam 2 CBET multiplier
      #pragma omp atomic write
      W_new[i+1][ix][iz] =W[i+1][ix][iz]*exp(-1*W[i][ix][iz]*dkmag[i+1][marker[(i+1)*ncrossings+n2]][markerCopy[(i+1)*ncrossings+n2]]*gain1/sqrt(epsilon));
      //beam 1 CBET multiplier
      #pragma omp atomic write
      W_new[i][ix][iz] =W[i][ix][iz]*exp(1*W[i+1][ix][iz]*dkmag[i][j][m]*gain1/sqrt(epsilon));
 //   }
  }
}
//updating CBET Intensities
void update_CBETPar()
{
  //storage array to allow for data level parallelization
  //cout << "Updating CBET Intensities" << endl;
  double convergeArr[threads]{0.0};

  //perform calculation for each crossing of each ray calculated for the beams in launch rays
  for(int i = 0; i < nbeams-1;i++)
  {
    #pragma omp parallel for num_threads(threads)
    for(int j = 0; j < nrays; j++)
    {
      for(int m = 0; m < ncrossings; m++)
      {
        //get the location of the mth crossing for ray j
        int boxx = *vec4D(boxes,i,j,m,0, nrays, ncrossings, 2);//x coor
        int boxz = *vec4D(boxes,i,j,m,1, nrays, ncrossings, 2);//z coor
        if(boxx == 0 || boxz == 0)//no more crossings for ray j
        {
          break;
        }
        int ix = boxx - 1;
        int iz = boxz - 1;
        //if there is a ray intersection at (ix, iz) and (ix, iz) is a valid coordinate
        if(ix != -1 && iz != -1 && intersections[ix][iz] != 0)
        {

          //First have to find all of the rays present from both beams in the coordinate (To generalize for N beams, encapsulate in a for loop)
          double tempStore = maxDev;

          int* marker = new int[nbeams*ncrossings];//saves time, assuming the majority of rays will have close to ncrossing crossings
          int* markerCopy = new int[nbeams*ncrossings];
          //printf("Ray %d at crossing %d: (%d, %d)\n",j, m, ix,iz);
          //find the number of beams
          //allocate space for marker arrays
          int cnt1 = 0;
          int cnt2 = 0;
          //iterate over marked at (ix, iz) to find all incident rays
          int ray1num = 0;
          for(int q = 0; q < ncrossings; q++)
          {
            int temp1 = *vec4D(marked,i,ix,iz,q,nx,nz,ncrossings);
            if(temp1)
            {
              marker[cnt1] = temp1-1;
              markerCopy[cnt1] = temp1-1;
              if(temp1-1 == j)
              {
                ray1num = cnt1;//find the index of ray j in the marker array
              }
              cnt1++;
            }
            temp1 = *vec4D(marked,i+1,ix,iz,q,nx,nz,ncrossings);
            if(temp1)
            {
              marker[ncrossings+cnt2] = temp1-1;
              markerCopy[ncrossings+cnt2] = temp1-1;
              cnt2++;
            }
          }
          int numrays[2] = {cnt1, cnt2};
          //locate ray j's index within the zone
         /* for(int r = 0; r < numrays[i]; r++)
          {

            if(marker[i*ncrossings+r] == j)
            {

              ray1num = r;
              break;
            }
          }
*/
          //find a ray from beam 2 that intersects at (ix, iz)
          for(int n = 0; n < numrays[i+1]; n++)
          {
            for(int q = 0; q < ncrossings;q++)
            {
              int ix2 = *vec4D(boxes,i+1,marker[(i+1)*ncrossings+n],q,0, nrays, ncrossings, 2) - 1;
              int iz2 = *vec4D(boxes,i+1,marker[(i+1)*ncrossings+n],q,1, nrays, ncrossings, 2)  - 1;
              if(ix == ix2 && iz == iz2)//if q has a valid crossing in (ix, iz)
              {
                markerCopy[(i+1)*ncrossings+n] = q;
                break;
              }
            }
          }

          calculateMult(ix, iz, i, j, m, numrays, marker, markerCopy);//calculate CBET multipliers using ray parameters

          //calculate the additive change in CBET field to be applied
          double additive_change[2];
          additive_change[0] = (-1.0 * (1.0 - (W_new[i][ix][iz]/W[i][ix][iz])) * abs(i_b[i][ix][iz]));// - i_b_prev[i][ix][iz])*(i_b[i][ix][iz] != 0));//change to be applied to beam 1 rays
          additive_change[1] = (-1.0 * (1.0 - (W_new[i+1][ix][iz]/W[i+1][ix][iz])) *  abs(i_b[i+1][ix][iz]));// - i_b_prev[i+1][ix][iz]));//change to be applied to beam 2 rays
          int kill = 0;
          if(i_b_new[i][ix][iz] < abs(additive_change[0]) && additive_change[0] < 0 || abs(i_b[i][ix][iz]) <= 1e-10)//if beam 1's intensity is/will be below zero, kill the rays in this zone and any zones that they will pass through in the future, only happens at high intensities
          {
            additive_change[1] = i_b_new[i][ix][iz] * (abs(i_b_new[i][ix][iz]) <= 1e-10);//if beam 1 still has energy, that will be transferred to beam 2 to conserve energy
            #pragma omp atomic write
            i_b_new[i][ix][iz] = 0;
            kill = 1;//stores that beam 1 has been killed in this sector
          }else
          {
            #pragma omp atomic update
            i_b_new[i][ix][iz] += additive_change[0];//standard case
          }
          if(abs(i_b_new[i][ix][iz]) < 1e-10 && additive_change[1] > 1)
          {
            printf("%e %e\n",i_b_new[i][ix][iz], additive_change[1]);
          }
          tempStore = fmax(abs(additive_change[0]/(i_b[i][ix][iz]+(abs(i_b_new[i][ix][iz]) <= 1e-10))), tempStore);//variable used to store convergence
          tempStore = fmax(abs(additive_change[1]/(i_b[i+1][ix][iz]+(abs(i_b_new[i+1][ix][iz]) <= 1e-10))), tempStore);//variable used to store convergence
          if(tempStore  > 1e10)
          {
       //     printf("TempStore 1 %e\n", tempStore);
          }
          #pragma omp atomic update
          i_b_new[i+1][ix][iz] += additive_change[1];//apply update to beam 1

          if(isnan(i_b_new[i][ix][iz]))
          {
          //  printf("Old: %f || Change %f || Kill Status: %d\n", i_b[i][ix][iz], additive_change[0], kill);
          //  printf("W_new[i][ix][iz]: %f || W[i][ix][iz] %f || Prev: %d\n", W_new[i][ix][iz],W[i][ix][iz], i_b_prev[i][ix][iz]);
          }
          int x_prev = -1;//x[ix];
          int z_prev = -1;//z[iz];
          //for every future crossing of ray j
          for(int l = m+1; l < ncrossings;l++)
          {
            int ix_next_i = *vec4D(boxes,i,j,l,0, nrays, ncrossings, 2) - 1;
            int iz_next_i = *vec4D(boxes,i,j,l,1, nrays, ncrossings, 2) - 1;
            //if l is a valid crossing
            if(ix_next_i != -1 && iz_next_i != -1)
            {
              double x_curr = x[ix_next_i];
              double z_curr = z[iz_next_i];
              //if l is not a previously updated crossing (within this loop)
              if(ix_next_i != x_prev || iz_next_i != z_prev)//x_curr != x_prev && z_curr != z_prev)
              {
                double add = additive_change[0]*present[i][ix][iz]/present[i][ix_next_i][iz_next_i]; // Change to be applied
                double comp = add/(i_b_new[i][ix_next_i][iz_next_i]+(abs(i_b_new[i][ix_next_i][iz_next_i]) < 1e-10));//Variable simply for the sake of clarity when updating tempstore
                if(tempStore == comp && comp > 1e10)
                {
           //       printf("TempStore 2 %e || %e : %e || %e\n", tempStore, i_b_new[i][ix_next_i][iz_next_i], i_b_new[i+1][ix_next_i][iz_next_i], add);
                }
                int killspot;
                if(abs(add) > abs(i_b_new[i][ix_next_i][iz_next_i]) && add < 0 || abs(i_b_new[i][ix][iz]) < 1e-10 && add < 0)
                {
                  //printf("\tKilled Crossing %d: (%d, %d)\n",l, ix_next_i,iz_next_i);
                  #pragma omp atomic write
                  i_b_new[i][ix_next_i][iz_next_i] = 0;
                  //#pragma omp atomic write
                  //i_b[i][ix_next_i][iz_next_i] = 0;//set initial intensity to 0 if kill == 1`
                }else{
                //  printf("\tAdjusted Crossing %d: (%d, %d) %e\n",l, ix_next_i,iz_next_i, add);
                  #pragma omp atomic update
                  i_b_new[i][ix_next_i][iz_next_i] += add;
                  tempStore = fmax(tempStore, comp);//store convergence variable

                }
              }
              x_prev = ix_next_i;//x_curr;
              z_prev = iz_next_i;//z_curr;
            }else
            {
              break;
            }
          }
         // printf("\tOther Ray Intersections:\n");
          int n2 = fmin(ray1num, numrays[i+1] - 1);//either stores the current ray in beam 1, or the last ray in beam 2, whichever is older, make sure no segfault
          //for every future crossing of ray n2, update intensity
          for(int l = markerCopy[(i+1)*ncrossings+n2] + 1; l < ncrossings;l++)
          {
            //find next grid coordinates

            int ix_next_other = *vec4D(boxes,i+1,marker[(i+1)*ncrossings+n2],l,0, nrays, ncrossings, 2) - 1;//[i+1] [marker[i+1][n2]] [l] [0]
            int iz_next_other = *vec4D(boxes,i+1,marker[(i+1)*ncrossings+n2],l,1, nrays, ncrossings, 2) - 1;//[i+1] [marker[i+1][n2]] [l] [1]
            //if it is a valid crossing
            if(ix_next_other != -1 && iz_next_other != -1)
            {
              double x_curr = x[ix_next_other];
              double z_curr = z[iz_next_other];
              //if this coordinate has not been updated in this loop
              if(ix_next_other != x_prev || iz_next_other != z_prev)//x_curr != x_prev && z_curr != z_prev)
              {

                double add = additive_change[1]*present[i][ix][iz]/present[i+1][ix_next_other][iz_next_other]; //change to be applied
                double comp = add/(i_b_new[i+1][ix_next_other][iz_next_other]+(i_b_new[i+1][ix_next_other][iz_next_other] < 1e-10)); // convergence variable
                tempStore = fmax(tempStore, comp);
                if(tempStore == comp && comp > 1e10)
                {
             //     printf("TempStore 3 %e\n", tempStore);
                }
                //printf("\tCrossing %d w/ ray %d: (%d, %d)\n",l, marker[(i+1)*ncrossings+n2],ix_next_other,iz_next_other);

                #pragma omp atomic update
                i_b_new[i+1][ix_next_other][iz_next_other] += add;//update downstream intensity
              }
              //Set up variables for next iteration
              x_prev = ix_next_other;//x_curr;
              z_prev = iz_next_other;//z_curr;
            }else
            {
              break;
            }
          }
          delete [] marker;
          delete [] markerCopy;
          convergeArr[omp_get_thread_num()] = max(abs(tempStore), convergeArr[omp_get_thread_num()]);
        }
      }
    }
  }
  for(int i = 0; i < threads; i++)
  {
    maxDev = fmax(maxDev, convergeArr[i]);
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
  //update_CBETPar();
  int cnt = 0;
  //CBET Iteration Loop
  do{
    prev = maxDev;
    maxDev = 0;
    //updateMult();
    //updateRays();
    update_CBETPar();
    cnt++;
    if(printCBETDiagnostics)
    {
      cout << "Iteration: " << cnt << "  Converging: " << maxDev << endl;
    }
  }
  while(abs(maxDev) > converge && cnt < 100 && iterate);// && cnt < maxIterations)
  
  if(printCBETDiagnostics)
  {
    cout << "Iteration: " << cnt << "  Converging: " << maxDev << endl;
  }
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
