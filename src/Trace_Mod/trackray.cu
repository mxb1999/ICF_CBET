#include <stdio.h>
#include <stdlib.h>
#include "Trace_interface.hpp"
#include "io_interface.hpp"
#include "cuda_help.hpp"
#include <stdarg.h>
#define MAX(i,j) ((i > j) ? i : j)
#define MIN(i,j) ((i < j) ? i : j)
__device__ double* edepcuComp;

//initializing necessary arrays for the calculations
__global__ void rayLaunchKernel(TrackConst val, TrackArrs arrs,rayinit* rays_cu, int* raypath)
{

  double dx_cu = val.dx_cu;
  double dz_cu = val.dz_cu;
  int nx_cu = val.nx_cu;
  int nz_cu = val.nz_cu;

  double xmax_cu = val.xmax_cu;
  double xmin_cu = val.xmin_cu;
  double zmax_cu = val.zmax_cu;
  double zmin_cu = val.zmin_cu;
  int nrays_cu = val.nrays_cu;
  int nbeams_cu = val.nbeams_cu;
  int ncrossings_cu = val.ncrossings_cu;
  int numstored_cu = val.numstored_cu;
  int nt_cu = val.nt_cu;
  double dt_cu = val.dt_cu;
  double omega_cu = val.omega_cu;
  double c_cu = val.c_cu;
  double ncrit_cu = val.ncrit_cu;
  
  
  double* dedendx_cu = arrs.dedendx_cu;
  double* dedendz_cu = arrs.dedendz_cu;
  double* x_cu = arrs.x_cu;
  double* z_cu = arrs.z_cu;
  double* crossesx_cu = arrs.crossesx_cu;
  double* crossesz_cu = arrs.crossesz_cu;
  double* edep_cu = arrs.edep_cu;
  double* wpe_cu = arrs.wpe_cu;
  int* boxes_cu = arrs.boxes_cu;

  int raynum = threadIdx.x;
  int beam = blockIdx.x;


  if(beam >= nbeams_cu || raynum >= nrays_cu)
  {
    return;
  }
  rayinit thisInit = *((rayinit*)rays_cu + beam*nrays_cu+ raynum);
  //double xinit = thisInit.xinit;
  //double zinit = thisInit.zinit;
  double kxinit = thisInit.kxinit;
  double kzinit = thisInit.kzinit;
  //Launch_Ray_XZ Array Declaration
  int thisx = 0;
  int thisz = 0;
  int thisx_0 = 0;
  int thisz_0 = 0;
  //double uray;
  //int //markingxprev;// = new double[nt]{0.0};
  //int //markingzprev;// = new double[nt]{0.0};
  //int //markingx;// = new double[nt]{0.0};
  //int //markingz;
  double myx = thisInit.xinit;
  double myz = thisInit.zinit;
  double myvx;
  double myvz;
  double mykx;
  double mykz;
  double myvxprev;
  double myvzprev;
  //determining the initial x grid index within the desired range to track the beam


  //determining the initial z grid index within the desired range to track the beam
  for(int i = 0;i < nz_cu;i++)
  {
    if(thisInit.zinit - z_cu[i] <= ((0.5+1.0e-10)*dz_cu + 1e-11) && thisInit.zinit - z_cu[i] >= -1*((0.5+1.0e-10)*dz_cu + 1e-11) )
    {

      thisz_0=i;
      break;
    }
  }

  for(int i = 0;i < nx_cu;i++)
  {

    double initX = thisInit.xinit;
    double currX = x_cu[i];
    if(initX - currX <= ((0.5+1.0e-10)*dx_cu + 1e-11) && initX - currX >= -1*((0.5+1.0e-10)*dx_cu + 1e-11) )
    {
      thisx_0=i;
      break;
    }
  }
  
  //determining the velocity characteristics of the ray based upon its initial position
  double k = sqrt((pow(omega_cu,2.0)-pow(vec2D_cu(wpe_cu,thisx_0, thisz_0, nz_cu),2.0))/pow(c_cu,2.0));
  double knorm = sqrt(pow(kxinit,2.0)+pow(kzinit,2.0));
  mykx=(kxinit/knorm)*k;			// Normalized value for the ray's initial k_x
  mykz=(kzinit/knorm)*k;			// Normalized value for the ray's initial k_z
  myvxprev = pow(c_cu,2.0)*mykx/omega_cu;                   // v_group, group velocity (dw/dk) from D(k,w).
  myvzprev =  pow(c_cu,2.0)*mykz/omega_cu;
  ////markingxprev = thisx_0;
  ////markingzprev = thisz_0;
  //__________Time Stepping__________
    int numcrossing = 0;

    //looping through time intervals
    for(int i = 1; i < nt_cu;i++)
    {
      if(raynum == 20 && beam == 1)
      {
        raypath[thisx*nz_cu + thisz] = 1;
      }
      double forcex = pow(c_cu,2.0)/(2.0*ncrit_cu)*vec2D_cu(dedendx_cu,thisx_0,thisz_0, nz_cu);
      double forcez = pow(c_cu,2.0)/(2.0*ncrit_cu)*vec2D_cu(dedendz_cu,thisx_0,thisz_0, nz_cu);


      myvx = myvxprev - forcex*dt_cu;
      myvz = myvzprev - forcez*dt_cu;
      myx += myvx*dt_cu;
      myz += myvz*dt_cu;

      int search_index_x = 1;
      int search_index_z = 1;
      int thisx_m = MAX(0, thisx_0-search_index_x );
      int thisx_p = MIN(nx_cu-1, thisx_0+search_index_x);
      int thisz_m = MAX(0, thisz_0-search_index_z);
      int thisz_p = MIN(nz_cu-1, thisz_0+search_index_z);
      //determining the current x index of the ray
      for(int j = thisx_m; j <= thisx_p;j++)
      {
        if ( myx - x_cu[j] <= ((0.5+1.0e-10)*dx_cu + 1e-12) && myx - x_cu[j] >= -1*((0.5+1.0e-10)*dx_cu + 1e-12))
        {
          thisx = j;
          break;
        }
      }
      if(raynum == 0 && beam == 0)
      {

       // printf("check 2 %d\n",i);
      }
      //determining the current z index of the ray
      for(int j = thisz_m; j <= thisz_p; j++)
      {
        if (myz - z_cu[j] <= ((0.5+1.0e-10)*dz_cu + 1e-12) && myz - z_cu[j] >= -1*((0.5+1.0e-10)*dz_cu + 1e-12))
        {

          thisz = j;
          break;
         }
      }
      //double linez[2]={thisInit.zinit, myz};
      //double linex[2]={thisInit.xinit, myx};
      int lastx = 10000;
      int lastz = 10000;
      //iterating through the selected portions of the x spatial tracking arrays
      //boxes_cu stores the spatial locations of each crossing of each ray
      //Marked = trajectory of a single ray, boxes_cu = coordinates of each ray intersection
      
      for(int j = thisx_m; j <= thisx_p;j++)
      {
        double currx = x_cu[j]-dx_cu/2;//crossing into 
        //if the ray is currently between within the desired caustic zone for a crossing
        if((myx > currx && thisInit.xinit <= (currx + 1e-10)) || (myx < currx && thisInit.xinit >= (currx- 1e-10)))
        {
          double m = (myz - thisInit.zinit)/(myx-thisInit.xinit);
          double b = myz - myx*m;
          double crossx = m*currx+b;
          //if the ray has moved since last update
          if(abs(crossx-lastz)>1.0e-20)
          {
            vec3DW_cu(crossesx_cu, beam, raynum, numcrossing, nrays_cu, ncrossings_cu,currx);
            vec3DW_cu(crossesz_cu, beam, raynum, numcrossing, nrays_cu, ncrossings_cu,crossx);
            //if ray is still within the grid
            if(myx < (xmax_cu+dx_cu/2 + 1e-10) && myx > (xmin_cu-dx_cu/2 - 1e-10))
            {
              //printf("(%f, %f), ([%f,%f], [%f, %f]) \n", myx, myz, dx_cu*thisx+ xmin_cu, dx_cu*thisx+dx_cu + xmin_cu, dz_cu*thisz+ zmin_cu, dz_cu*thisz+dz_cu+ zmin_cu);
              vec4DW_cu(boxes_cu, beam,raynum,numcrossing,0, nrays_cu, ncrossings_cu, 2, thisx+1);//[beam][raynum][numcrossing][0]
              vec4DW_cu(boxes_cu, beam,raynum,numcrossing,1, nrays_cu, ncrossings_cu, 2, thisz+1);//[beam][raynum][numcrossing][1]
              //vec4DW_cu(marked_cu, beam,raynum,thisx,thisz, nrays_cu, nx_cu, nz_cu, 1);
            }
            lastx = currx;
            numcrossing += 1;
            break;
          }
        }
      }
      if(raynum == 0 && beam == 0)
      {

        //printf("check 3 %d\n",i);
      }
      //iterating through the selected portions of the z spatial tracking arrays
      //Same idea as previous loop, but for the Z coordinate instead of x
        for(int j = thisz_m; j <= thisz_p;j++)//for [thisz_m, thisz_p] previous z locations, iterate through spatial locations centered on thisz
        {
          double currz = z_cu[j]-dz_cu/2;//center of the jth zone
          //printf("REE\n");
          if((myz > (currz) && thisInit.zinit < (currz + 1e-10)) || (myz < (currz) && thisInit.zinit > (currz - 1e-10)))//if myz is approximately equal to a zone crossing
          {

            double m = (myx - thisInit.xinit)/(myz-thisInit.zinit);
            double b = myx - myz*m;
            double crossz = m*currz+b;
            if(abs(crossz-lastx) > 1.0e-20)
            {
              vec3DW_cu(crossesx_cu, beam, raynum, numcrossing, nrays_cu, ncrossings_cu,crossz);
              vec3DW_cu(crossesz_cu, beam, raynum, numcrossing, nrays_cu, ncrossings_cu,currz);
              if(myz < (zmax_cu+dz_cu/2 +1e-10) && myz > (zmin_cu-dz_cu/2-1e-10))
              {
               vec4DW_cu(boxes_cu, beam,raynum,numcrossing,0, nrays_cu, ncrossings_cu, 2, thisx+1);//[beam][raynum][numcrossing][0]
               vec4DW_cu(boxes_cu, beam,raynum,numcrossing,1, nrays_cu, ncrossings_cu, 2, thisz+1);//[beam][raynum][numcrossing][1]
              // vec4DW_cu(marked_cu, beam,raynum,thisx,thisz, nrays_cu, nx_cu, nz_cu, 1);
              }
              lastz = currz;

              numcrossing += 1;
              break;
            }
          }
        }
        //Sets the "initial conditions" for the next iteration
        thisx_0 = thisx;
        thisz_0 = thisz;
        //markingx = thisx;
        //markingz = thisz;
        //Deposit energy due to the incident ray
  	    double increment = thisInit.urayinit;
        double xp = (myx - (x_cu[thisx]+dx_cu/2.0))/dx_cu;
        double zp = (myz - (z_cu[thisz]+dz_cu/2.0))/dz_cu;
        int xadd = (xp >= 0) ? 1 : -1;
        int zadd = (zp >= 0) ? 1 : -1;
        double dl = zp * zadd;
        double dm = xp * xadd;
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x+1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z+1)
        double a4 = dl*dm;			// red 		: (x+1, z+1)
        if((xadd != xadd) && (zadd != zadd))
        {

          break;
        }
        //printf("edep: %d %d %d %d %d\n",thisx+1, thisz+1, beam, raynum, RAYS);
        vec4DI_cu(edep_cu, beam, raynum, thisx+1, thisz+1, nrays_cu, nx_cu+2,nz_cu+2, a1*increment);	// blue
        vec4DI_cu(edep_cu, beam, raynum, thisx+xadd+1, thisz+1, nrays_cu, nx_cu+2,nz_cu+2, a2*increment);// green
        vec4DI_cu(edep_cu, beam, raynum, thisx+1, thisz+zadd+1, nrays_cu, nx_cu+2,nz_cu+2, a3*increment);// yellow
        vec4DI_cu(edep_cu, beam, raynum, thisx+xadd+1, thisz+zadd+1, nrays_cu, nx_cu+2,nz_cu+2, a4*increment);	// red
        thisInit.xinit = myx;
        thisInit.zinit = myz;
        myvxprev = myvx;
        myvzprev = myvz;
        ////markingxprev = //markingx;
        ////markingzprev = //markingz;
      if ( (myx < (xmin_cu-(dx_cu/2.0))) || (myx > (xmax_cu+(dx_cu/2.0))))
      {

        break;                  // "breaks" out of the i loop once the if condition is satisfied
      } else if ( (myz < (zmin_cu-(dz_cu/2.0))) || (myz > (zmax_cu+(dz_cu/2.0)))){
           // the "|" means "or" (symbol above the return key)

        break;
    }

  }
  //delete [] mytime;
  //delete [] nuei;
  //delete [] amplitude_norm;
  
}
template <typename T>
__device__ void flattenArray(T* in, T* out,int targetDim, int dims...)//function to compress edep_cu
{
  
  va_list args;
  va_start(args,dims);
  int* dimlist = new int[dims];
  int bprod=1;
  int bnum = 1;
  int aprod=1;
  bool check = false;
  for(int i = 1; i < dims; i++)
  {
    int dim = va_arg(args, int);
    if(!check)
    {
      bprod *= dim;
    }
    if(check)
    {
      aprod *= dim;
    }
    if(targetDim == i)
    {
      check = true;
    }
    dimlist[i] = dim;
  }
  int index = blockDim.x*blockIdx.x+threadIdx.x;
  //printf("%d %d\n", index, prod);
  int mainDim = dimlist[targetDim];
  bnum = bprod/mainDim;
  if(index >=  mainDim)
  {
    return;
  }
  //printf("%d\n", mainDim);

  for(int i = 0; i < bprod;i++)
  {
    for(int j = 0; j < aprod;j++)
    {
      out[(i+index)*aprod+j] += in[(i+index)*aprod+j];

    }
    //printf("%e\n", in[index + i]);
  }
  delete dimlist;
}

__global__ void locateInts(TrackConst val, TrackArrs arr, double* edep_flat)
{

  int nz_cu = val.nz_cu;
  int nx_cu = val.nx_cu;
  int nrays_cu = val.nrays_cu;
  int ncrossings_cu = val.ncrossings_cu;
  int numstored_cu = val.numstored_cu;
  int nbeams_cu = val.nbeams_cu;
  int* ints_cu = arr.ints_cu;
  int* boxes_cu = arr.boxes_cu;
  double* edep_cu = arr.edep_cu;
  int raynum = blockDim.x*blockIdx.x+threadIdx.x;
  
  int ix = raynum / nz_cu;
  int iz = raynum % nz_cu;
  if(ix < nx_cu + 2 && iz < nz_cu + 2 && raynum < nx_cu*nz_cu)
  {
    printf("%d\n", ix, iz, raynum);
    double acc0 = 0.0;
    double acc1 = 0.0;
    for(int i = 0; i < nrays_cu;i++)
    {
      acc0 += vec4D_cu(edep_cu, 0,i,ix,iz, nrays_cu, nx_cu+2,nz_cu+2);
      acc1 += vec4D_cu(edep_cu, 1,i,ix,iz, nrays_cu, nx_cu+2,nz_cu+2);
    }
    vec3DW_cu(edep_flat, 0,ix,iz,nx_cu+2,nz_cu+2,acc0);
    vec3DW_cu(edep_flat, 1,ix,iz,nx_cu+2,nz_cu+2,acc1);
  }

  //flattenEdep(edep_cu, edep_flat, 1, 4, nbeams_cu, nrays_cu, (nx_cu+2), (nz_cu+2));
  //collapseEdep_CU();
  //Locate all intersections
  if(raynum >= nrays_cu*nbeams_cu)
  {
    return;
  }
  int beam = raynum/nrays_cu;
  //for each crossing of raynum
  for(int i = 0; i < ncrossings_cu; i++)
  {
    //spatial coordinates
    int rx = vec4D_cu(boxes_cu, beam,raynum,i,0, nrays_cu, ncrossings_cu, 2);
    int rz = vec4D_cu(boxes_cu, beam,raynum,i,0, nrays_cu, ncrossings_cu, 2);
    int cnt = 0;
    if(!rx || !rz)
      break;
    rx--;
    rz--;
    for(int n = i+1; n < nbeams_cu;n++)
    {
      for(int j = 0; j < nrays_cu; j++)
      {
        for(int q = 0; q < ncrossings_cu; q++)
        {
          int cx = vec4D_cu(boxes_cu, n,j,q,0, nrays_cu, ncrossings_cu, 2);
          int cz = vec4D_cu(boxes_cu, n,j,q,0, nrays_cu, ncrossings_cu, 2);
          if(cx == rx && cz == rz)
          {
            vec4DW_cu(ints_cu, beam, raynum, i,cnt, nrays_cu, ncrossings_cu, numstored_cu, j+1);
            cnt++;
          }
          if(cnt >= numstored_cu)
          {
            break;
          }
        }
        if(cnt >= numstored_cu)
        {
          break;
        }
      }
      if(cnt >= numstored_cu)
      {
        break;
      }
    }
  }
}


#define B 2
void LaunchCUDARays(rayinit* rays)
{ 

  TrackConst constVals = *deviceTrackConst(0);
  TrackArrs constArrs = *deviceTrackArrs(0);
  rayLaunchKernel<<<2,nrays>>>(constVals,constArrs,rays,raypath);
  cudaDeviceSynchronize();
  cudaMallocManaged(&edep_flat, sizeof(double)*nbeams*(nx+2)*(nz+2));
  locateInts<<<4,1024>>>(constVals,constArrs,edep_flat);
  cudaDeviceSynchronize();
  for(int i = 0; i < nrays; i++)
  {
    for(int j = 0; j < nx+2;j++)
    {
      for(int m = 0; m < nz+2;m++)
      {
        double edep41 = vec4D(edep, 0,i,j,m,nrays, nx+2, nz+2);
        double edep42 = vec4D(edep, 1,i,j,m,nrays, nx+2, nz+2);
        double edep31 = vec3D(edep_flat, 0,j,m,nx+2, nz+2);
        double edep32 = vec3D(edep_flat, 1,j,m,nx+2, nz+2);
        if(abs(edep31) > 0 || abs(edep32) > 0)
        {
          //printf("3d %e %e\n", edep31, edep32); 
        }
        if(abs(edep41) > 0 || abs(edep42) > 0)
        {
          //printf("4d %e %e\n", edep41, edep42); 
        }
      }
    }
  }  
  cudaFree(edep);//free arrays no longer needed
  cudaFree(dedendx);
  cudaFree(dedendz);
}
