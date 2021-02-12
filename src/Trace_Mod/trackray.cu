#include <stdio.h>
#include <stdlib.h>
#include "Trace_interface.hpp"
#include "io_interface.hpp"
#include "cuda_help.hpp"
#include "cuda_var_init.hpp"
#define MAX(i,j) ((i > j) ? i : j)
#define MIN(i,j) ((i < j) ? i : j)

//initializing necessary arrays for the calculation
__global__ void rayLaunchKernel()
{
  int raynum = blockDim.x*blockIdx.x+threadIdx.x;
  if(raynum >= nrays_cu)
  {
    return;
  }
  rayinit thisInit = *((rayinit*)rays_cu + raynum);
  //double xinit = thisInit.xinit;
  //double zinit = thisInit.zinit;
  double kxinit = thisInit.kxinit;
  double kzinit = thisInit.kzinit;
  int beam = thisInit.beam;
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
  double myx;
  double myz;
  double myvx;
  double myvz;
  double mykx;
  double mykz;
  double myvxprev;
  double myvzprev;

  //determining the initial x grid index within the desired range to track the beam
  for(int i = 0;i < nx_cu;i++)
  {
    if(thisInit.xinit - x_cu[i] <= ((0.5+1.0e-10)*dx_cu + 1e-11) && thisInit.xinit - x_cu[i] >= -1*((0.5+1.0e-10)*dx_cu + 1e-11) )
    {
      thisx_0=i;
      break;
    }
  }

  //determining the initial z grid index within the desired range to track the beam
  for(int i = 0;i < nz_cu;i++)
  {
    if(thisInit.zinit - z_cu[i] <= ((0.5+1.0e-10)*dz_cu + 1e-11) && thisInit.zinit - z_cu[i] >= -1*((0.5+1.0e-10)*dz_cu + 1e-11) )
    {

      thisz_0=i;
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
      double forcex = pow(c_cu,2.0)/(2.0*ncrit_cu)*vec2D_cu(dedendx_cu,thisx_0,thisz_0, nz_cu);
      double forcez = pow(c_cu,2.0)/(2.0*ncrit_cu)*vec2D_cu(dedendz_cu,thisx_0,thisz_0, nz_cu);
      myvz = myvzprev - forcex*dt_cu;
      myvx = myvxprev - forcez*dt_cu;
      myx = thisInit.xinit + myvx*dt_cu;
      myz = thisInit.zinit + myvz*dt_cu;
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
        double currx = x_cu[j];//-dx_cu/2;//crossing into 
        //if the ray is currently between within the desired caustic zone for a crossing
        if((myx > currx && thisInit.xinit <= (currx + 1e-10)) || (myx < currx && thisInit.xinit >= (currx- 1e-10)))
        {
          double m = (myz - thisInit.zinit)/(myx-thisInit.xinit);
          double b = myz - myx*m;
          double crossx = m*currx+b;
          //if the ray has moved since last update
          if(abs(crossx-lastz)>1.0e-20)
          {
            vec3DW_cu(crossx_cu, beam, raynum, numcrossing, nrays_cu, ncrossings_cu,currx);
            vec3DW_cu(crossz_cu, beam, raynum, numcrossing, nrays_cu, ncrossings_cu,crossx);
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

      //iterating through the selected portions of the z spatial tracking arrays
      //Same idea as previous loop, but for the Z coordinate instead of x
        for(int j = thisz_m; j <= thisz_p;j++)//for [thisz_m, thisz_p] previous z locations, iterate through spatial locations centered on thisz
        {
          double currz = z_cu[j];//-dz_cu/2;//center of the jth zone
          //printf("REE\n");
          if((myz > (currz) && thisInit.zinit < (currz + 1e-10)) || (myz < (currz) && thisInit.zinit > (currz - 1e-10)))//if myz is approximately equal to a zone crossing
          {

            double m = (myx - thisInit.xinit)/(myz-thisInit.zinit);
            double b = myx - myz*m;
            double crossz = m*currz+b;
            if(abs(crossz-lastx) > 1.0e-20)
            {
              vec3DW_cu(crossx_cu, beam, raynum, numcrossing, nrays_cu, ncrossings_cu,crossz);
              vec3DW_cu(crossz_cu, beam, raynum, numcrossing, nrays_cu, ncrossings_cu,currz);
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
          printf("NaN Error: Error in Grid Interpolation \n");
          break;
        }
        vec4DI_cu(edep_cu, beam, raynum, thisx+1, thisz+1, nrays_cu, nx_cu+2,nz_cu+2, a1*increment);	// blue
        vec4DI_cu(edep_cu, beam, raynum, thisx+xadd+1, thisz+1, nrays_cu, nx_cu+2,nz_cu+2, a2*increment);// green
        vec4DI_cu(edep_cu, beam, raynum, thisx+1, thisz+zadd+1, nrays_cu, nx_cu+2,nz_cu+2, a3*increment);// yellow
        vec4DI_cu(edep_cu, beam, raynum, thisx+xadd+1, thisz+zadd+1, nrays_cu, nx_cu+2,nz_cu+2, a4*increment);	// red
        thisInit.xinit = myx;
        myvxprev = myvx;
        thisInit.zinit = myz;
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
__device__ void collapseEdep_CU(double* edepcu, double* edepcuComp, int beam)//function to compress edep
{
  int indexX = blockDim.x*blockIdx.x+threadIdx.x;
  int indexZ = blockDim.y*blockIdx.y+threadIdx.y;
  if(indexX >= nx_cu || indexZ >= nz_cu)
  {
    return;
  }
  double accumulator = 0.0;
  for(int i = 0; i < nrays_cu;i++)
  {
    accumulator += vec4D_cu(edep_cu, beam, i, indexX, indexZ, nrays_cu, nx_cu+2, nz_cu+2); 
  }
  vec3DW_cu(edepcuComp, beam, indexX, indexZ, nx_cu+2, nz_cu+2, accumulator);

}
__global__ void locateInts(int* ints, int beam,double* edepcu, double* edepcuComp)
{
  collapseEdep_CU(edep_cu, edepcuComp, beam);
  //Locate all intersections
  int raynum = blockDim.x*blockIdx.x+threadIdx.x;
  if(raynum >= nrays_cu)
  {
    return;
  }
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
    for(int n = 0; n < nbeams_cu;n++)
    {
      for(int j = 0; j < nrays_cu; j++)
      {
        if(vec4D_cu(marked_cu, n,j,rx,rz,nrays_cu, nx_cu, nz_cu))
        {
          vec4DW_cu(ints, beam, raynum, i, cnt++, nrays_cu, ncrossings_cu, numstored_cu,1);
        }
        if(cnt >= numstored_cu)
        {
          break;
        }
      }
    }
  }
}
#define T 512
#define B 1024
void LaunchCUDARays(GConfig gpu,rayinit* rays)
{ int n = 0;
  hostP[0] = rays;
  sizes[0] = nbeams*nrays*sizeof(rayinit);
  for(int* i = (int*)hostP; i != NULL; i++)
  {
    n++;
  }
  printf("n %d\n",n);
  gpu.deviceDataTransfer(hostP, devP, sizes, n, 0);
  rayLaunchKernel<<<T,B>>>();
}
