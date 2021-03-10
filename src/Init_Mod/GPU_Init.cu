#include "cuda_var_init.hpp"
/*__device__ double mi_kg_cu = 10230.0*me_cu;	   // Mass of ion in kg
__device__ double mi_cu = 10230*(1.0e3*me_cu);          // Mass of ion in g
__device__ double uray_mult_cu = intensity_cu*(courant_mult_cu)*pow(double(rays_per_zone_cu),-1.0); //multiplier which determines intensity deposited in a given zone
__device__ double dz_cu = (zmax_cu-zmin_cu)/(nz_cu-1);
__device__ double dx_cu = (xmax_cu-xmin_cu)/(nx_cu-1);
__device__ int nrays_cu= int(rays_per_zone_cu*(beam_max_z_cu-beam_min_z_cu)/dz_cu)+0;//number of rays per beam
__device__ double dt_cu=courant_mult_cu*fmin(dx_cu,dz_cu)/c;//time stepping
__device__ int nt_cu=int(pow(courant_mult_cu,-1.0)*fmax(nx_cu,nz_cu)*2.0)+1;//number of time steps to track for a given ray
__device__ int numstored_cu = nx_cu*6;//number of rays stored per grid zone
__device__ int ncrossings_cu = nx_cu * 3;//max number of ray crossings that can be stored per ray
__device__ double freq_cu = c_cu/lambda_cu;		// frequency of light, in Hz
__device__ double omega_cu = 2*pi_cu*freq_cu;	// frequency of light, in rad/s
__device__ double ncrit_cu = 1e-6*(pow(omega_cu,2.0)*me_cu*e0_cu/pow(ec_cu,2.0));*/

void gpuInit()
{
  cudaMallocManaged(&intersections, sizeof(int)*GRID);
  //intersections = new int[GRID]{0}; //nx nz
  cudaMallocManaged(&marked, sizeof(int)*GRID*RAYS);
  //marked = new int[GRID*nbeams*numstored]{0}; //nx nz nrays nbeams
  cudaMallocManaged(&dedendx, sizeof(int)*GRID);
  //dedendx = new double[GRID]; //nx nz
  cudaMallocManaged(&dedendz, sizeof(int)*GRID);
  //dedendz = new double[GRID]; //nx nz
  cudaMallocManaged(&x, sizeof(int)*nx);
  //x = new double[nx]{0.0}; //nx nz
  cudaMallocManaged(&z, sizeof(int)*nz);
  //z = new double[nz]{0.0}; //nx nz
  cudaMallocManaged(&eden, sizeof(int)*GRID);
  //eden = new double[GRID]; //nx nz
  cudaMallocManaged(&edep, sizeof(int)*nbeams * (nx+2) * (nz+2));
  //edep = new double[nbeams * nx+2 * nz+2]{0.0}; //nx+2 nz+2 nrays
  cudaMallocManaged(&present, sizeof(int)*GRID);
  //present = new int[nbeams*GRID]; //nx nz nbeams
  cudaMallocManaged(&machnum, sizeof(int)*GRID);
  //machnum = new double[GRID]; //nx nz
  cudaMallocManaged(&boxes, sizeof(int)*RAYS*ncrossings);
  //boxes = new int[nbeams*nrays*ncrossings*2]{0}; //nbeams nrays ncrossings 2
  cudaMallocManaged(&u_flow, sizeof(int)*GRID);
  //u_flow = new double[GRID]; //nx nz
  cudaMallocManaged(&dkx, sizeof(int)*RAYS*ncrossings);
  //dkx = new double[nbeams*nrays*ncrossings]; //nbeams nrays 2
  cudaMallocManaged(&dkz, sizeof(int)*RAYS*ncrossings);
  //dkz = new double[nbeams*nrays*ncrossings]; //nbeams nrays 2
  cudaMallocManaged(&dkmag, sizeof(int)*RAYS*ncrossings);
  //dkmag = new double[nbeams*nrays*ncrossings]; //nbeams nrays 2
  cudaMallocManaged(&W, sizeof(int)*GRID*nbeams);
  //W = new double[nbeams*GRID];//nx nz
  cudaMallocManaged(&W_new, sizeof(int)*GRID*nbeams);
  //W_new = new double[nbeams*GRID];//nx nz
  cudaMallocManaged(&wpe, sizeof(int)*GRID);
  //wpe = new double[GRID]; //nx nz
  cudaMallocManaged(&crossesz, sizeof(int)*RAYS*ncrossings);
  //crossesz = new double[nbeams*nrays*ncrossings]; //nbeams nrays ncrossings
  cudaMallocManaged(&crossesx, sizeof(int)*RAYS*ncrossings);
  //crossesx = new double[nbeams*nrays*ncrossings]; //nbeams nrays ncrossings
  cudaMallocManaged(&ints, sizeof(int)*RAYS*nrays*);
  //ints = new int[nbeams*nrays*ncrossings*numstored]; //nbeams nrays ncrossings
  cudaMallocManaged(&raypath, sizeof(int)*GRID);
  //raypath = new int[GRID];
  cudaMallocManaged(&mult, sizeof(int)*RAYS*ncrossings);
  //mult = new double[nbeams*nrays*ncrossings]{1.0};
  for(int i = 0; i < RAYS*ncrossings;i++)
  {
    mult[i] = 1;
  }
  gain2arr = new double[GRID];
  gain1arr = new double[GRID];
  mag = new double[GRID];

  auto check1 = chrono::high_resolution_clock::now();
  //Calculating the initial energy density, wpe, and machnum values
  span(x, xmin, xmax, nx);
  span(z, zmin, zmax, nz);

  for(int i = 0; i < nx;i++)
  {
    for(int j = 0; j < nz;j++)
    {
      //eden[i][j] = temp[i];
      double temp = fmax(0.0,((0.3*ncrit-0.1*ncrit)/(xmax-xmin))*(x[i]-xmin)+(0.1*ncrit));
      vec2DW(eden,i,j,nz,temp);

      temp = sqrt(vec2D(eden,i,j,nz)*1e6*pow(ec,2.0)/(me*e0));
      vec2DW(wpe,i,j,nz, temp);

      temp = fmax(0.0,(((-0.4)-(-2.4))/(xmax-xmin))*(x[i]-xmin))+(-2.4);
      vec2DW(machnum,i,j,nz, temp);
    }
  }
  for(int i = 0; i < nx-1; i++)
  {
    for(int j = 0; j < nz-1; j++)
    {
      double tempx = (vec2D(eden,i+1,j,nz)-vec2D(eden,i,j,nz))/(x[i+1]-x[i]);
      double tempz = (vec2D(eden,i,j+1,nz)-vec2D(eden,i,j,nz))/(z[j+1]-z[j]);
      vec2DW(dedendx,i,j,nz, tempx);
      vec2DW(dedendz,i,j,nz, tempz);
    }
  }
  for(int i = 0; i < max(nx,nz);i++)
  {
    if(i < nx)
    {
      double temp = vec2D(dedendz,i,nz-2,nz);
      vec2DW(dedendz,i,nz-1,nz,temp);

    }
    if(i < nx-1)
    {
      double tempx = (vec2D(eden,i+1,nz-1,nz)-vec2D(eden,i,nz-1,nz))/(x[i+1]-x[i]);
      vec2DW(dedendx,i,nz-1,nz, tempx);
    }
    if(i < nz-1)
    {
      double tempz = (vec2D(eden,nx-1,i+1,nz)-vec2D(eden,nx-1,i,nz))/(z[i+1]-z[i]);
      vec2DW(dedendz,nx-1,i,nz, tempz);
    }
    if(i < nz)
    {

      double temp = vec2D(dedendx,nx-2,i,nz);
      vec2DW(dedendx,nx-1,i,nz,temp);
    }
  }
}   