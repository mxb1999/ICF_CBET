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
        
        if(temp != convergeStore[numT])
        {
         // printf("Converge: %e %e %d %d %d %e %e\n",currE, sqrt(8*pi*i_b[i][currx][currz]/c), i, j, m, mult[(i*nrays+j)*ncrossings + m], convergeStore[numT]);
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
    updateMult();
    updateRays();
    //update_CBETPar();
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
