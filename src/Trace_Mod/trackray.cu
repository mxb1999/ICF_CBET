#include <stdio.h>
#include <stdlib.h>
#include "Trace_interface.hpp"
#include "io_interface.hpp"
 



//initializing necessary arrays for the calculation
__global__ void rayLaunch(rayinit* init, double* edepcu, double* force, double* crossx_cu, double* crossz_cu, int nrays)
{
  int raynum = blockDim.x*blockIdx.x+threadIdx.x;
  if(raynum >= nrays)
  {
    return;
  }
  rayinit thisInit = init[raynum];
  double xinit = thisInit.xinit;
  double zinit = thisInit.zinit;
  double kxinit = thisInit.kxinit;
  double kzinit = thisInit.kzinit;
  int beam = thisInit.beam;
  double edeninit = thisInit.edeninit;
  double wpeInit = thisInit.wpeInit;
  //Launch_Ray_XZ Array Declaration
  int thisx = 0;
  int thisz = 0;
  int thisx_0 = 0;
  int thisz_0 = 0;
  double uray;
  int markingxprev;// = new double[nt]{0.0};
  int markingzprev;// = new double[nt]{0.0};
  int markingx;// = new double[nt]{0.0};
  int markingz;
  double myx;
  double myz;
  double myvx;
  double myvz;
  double mykx;
  double mykz;
  double myvxprev;
  double myvzprev;

  //determining the initial x grid index within the desired range to track the beam
  for(int i = 0;i < nx;i++)
  {
    if(thisInit.xinit - x_cu[i] <= ((0.5+1.0e-10)*dx + 1e-11) && thisInit.xinit - x_cu[i] >= -1*((0.5+1.0e-10)*dx + 1e-11) )
    {
      thisx_0=i;
      break;
    }
  }

  //determining the initial z grid index within the desired range to track the beam
  for(int i = 0;i < nz;i++)
  {
    if(thisInit.zinit - z_cu[i] <= ((0.5+1.0e-10)*dz + 1e-11) && thisInit.zinit - z_cu[i] >= -1*((0.5+1.0e-10)*dz + 1e-11) )
    {

      thisz_0=i;
      break;
    }
  }
  //determining the velocity characteristics of the ray based upon its initial position
  double k = sqrt((pow(omega,2.0)-pow(wpeInit,2.0))/pow(c,2.0));
  double knorm = sqrt(pow(kx_init,2.0)+pow(kz_init,2.0));
  mykx=(kx_init/knorm)*k;			// Normalized value for the ray's initial k_x
  mykz=(kz_init/knorm)*k;			// Normalized value for the ray's initial k_z
  myvxprev = pow(c,2.0)*mykx/omega;                   // v_group, group velocity (dw/dk) from D(k,w).
  myvzprev =  pow(c,2.0)*mykz/omega;
  markingxprev = thisx_0;
  markingzprev = thisz_0;
  //__________Time Stepping__________
    int numcrossing = 0;
    //looping through time intervals
    for(int i = 1; i < nt;i++)
    {
      myvz = myvzprev - pow(c,2.0)/(2.0*ncrit)*vec3D(force, thisx_0, thisz_0,0, nz, 2)*dt;
      myvx = myvxprev - pow(c,2.0)/(2.0*ncrit)*vec3D(force, thisx_0, thisz_0,1, nz, 2)*dt;
      myx = thisInit.xinit + myvx*dt;
      myz = thisInit.zinit + myvz*dt;
      int search_index_x = 1;
      int search_index_z = 1;
      int thisx_m = fmax(0, thisx_0-search_index_x );
      int thisx_p = fmin(nx-1, thisx_0+search_index_x);
      int thisz_m = fmax(0, thisz_0-search_index_z);
      int thisz_p = fmin(nz-1, thisz_0+search_index_z);
      //determining the current x index of the ray
      for(int j = thisx_m; j <= thisx_p;j++)
      {
        if ( myx - x_cu[j] <= ((0.5+1.0e-10)*dx + 1e-12) && myx - x_cu[j] >= -1*((0.5+1.0e-10)*dx + 1e-12))
        {
          thisx = j;
          break;
        }
      }

      //determining the current z index of the ray
      for(int j = thisz_m; j <= thisz_p; j++)
      {
        if (myz - z_cu[j] <= ((0.5+1.0e-10)*dz + 1e-12) && myz - z_cu[j] >= -1*((0.5+1.0e-10)*dz + 1e-12))
        {

          thisz = j;
          break;
         }
      }
      double linez[2]={thisInit.zinit, myz};
      double linex[2]={thisInit.xinit, myx};
      int lastx = 10000;
      int lastz = 10000;
      //iterating through the selected portions of the x spatial tracking arrays
      //boxes_cu stores the spatial locations of each crossing of each ray
      //Marked = trajectory of a single ray, boxes_cu = coordinates of each ray intersection
      for(int j = thisx_m; j <= thisx_p;j++)
      {
        double currx = x_cu[j];//-dx/2;//crossing into 
        //if the ray is currently between within the desired caustic zone for a crossing
        if((myx > currx && thisInit.xinit <= (currx + 1e-10)) || (myx < currx && thisInit.xinit >= (currx- 1e-10)))
        {
          double m = (myz - thisInit.zinit)/(myx-thisInit.xinit);
          double b = myz - myx*m;
          double crossx = m*currx+b;
          //if the ray has moved since last update
          if(abs(crossx-lastz)>1.0e-20)
          {
            vec3DW(crossx_cu, beam, raynum, numcrossing, nrays, ncrossings,currx);
            vec3DW(crossz_cu, beam, raynum, numcrossing, nrays, ncrossings,crossx);
            //if ray is still within the grid
            if(myx < (xmax+dx/2 + 1e-10) && myx > (xmin-dx/2 - 1e-10))
            {
              //printf("(%f, %f), ([%f,%f], [%f, %f]) \n", myx, myz, dx*thisx+ xmin, dx*thisx+dx + xmin, dz*thisz+ zmin, dz*thisz+dz+ zmin);
              vec4DW(boxes_cu, beam,raynum,numcrossing,0, nrays, ncrossings, 2, thisx+1);//[beam][raynum][numcrossing][0]
              vec4DW(boxes_cu, beam,raynum,numcrossing,1, nrays, ncrossings, 2, thisz+1);//[beam][raynum][numcrossing][1]
              vec4DW(marked_cu, beam,raynum,thisx,thisz, nrays, nx, nz, 1);
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
          double currz = z_cu[j];//-dz/2;//center of the jth zone
          //printf("REE\n");
          if((myz > (currz) && thisInit.zinit < (currz + 1e-10)) || (myz < (currz) && thisInit.zinit > (currz - 1e-10)))//if myz is approximately equal to a zone crossing
          {

            double m = (myx - thisInit.xinit)/(myz-thisInit.zinit);
            double b = myx - myz*m;
            double crossz = m*currz+b;
            if(abs(crossz-lastx) > 1.0e-20)
            {
              vec3DW(crossz_cu, beam, raynum, numcrossing, nrays, ncrossings,crossx);
            vec3DW(crossx_cu, beam, raynum, numcrossing, nrays, ncrossings,currz);
              if(myz < (zmax+dz/2 +1e-10) && myz > (zmin-dz/2-1e-10))
              {
               vec4DW(boxes_cu, beam,raynum,numcrossing,0, nrays, ncrossings, 2, thisx+1);//[beam][raynum][numcrossing][0]
               vec4DW(boxes_cu, beam,raynum,numcrossing,1, nrays, ncrossings, 2, thisz+1);//[beam][raynum][numcrossing][1]
               vec4DW(marked_cu, beam,raynum,thisx,thisz, nrays, nx, nz, 1);
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
        markingx = thisx;
        markingz = thisz;
        //Deposit energy due to the incident ray
  	    double increment = thisInit.urayinit;
        double xp = (myx - (x[thisx]+dx/2.0))/dx;
        double zp = (myz - (z[thisz]+dz/2.0))/dz;
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
        vec3DI(edepcu, beam, numray, thisx+1, thisz+1, nrays, nx,nz, a1*increment);	// blue
        vec3DI(edepcu, beam, numray, thisx+xadd+1, thisz+1, nrays, nx,nz, a2*increment);// green
        vec3DI(edepcu, beam, numray, thisx+1, thisz+zadd+1, nrays, nx,nz, a3*increment);// yellow
        vec3DI(edepcu, beam, numray, thisx+xadd+1, thisz+zadd+1, nrays, nx,nz, a4*increment);	// red
        thisInit.xinit = myx;
        myvxprev = myvx;
        thisInit.zinit = myz;
        myvzprev = myvz;
        markingxprev = markingx;
        markingzprev = markingz;
        
      if ( (myx < (xmin-(dx/2.0))) || (myx > (xmax+(dx/2.0))))
      {
        break;                  // "breaks" out of the i loop once the if condition is satisfied
      } else if ( (myz < (zmin-(dz/2.0))) || (myz > (zmax+(dz/2.0)))){
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
  if(indexX >= nx || indexZ >= nz)
  {
    return;
  }
  double accumulator = 0.0;
  for(int i = 0; i < nrays;i++)
  {
    accumulator += vec4D(edepcu, beam, i, indexX, indexZ, nrays, nx+2, nz+2); 
  }
  vec3DW(edepcuComp, beam, indexX, indexZ, nx+2, nz+2);

}
__global__ void locateInts(int* ints, int beam,double* edepcu, double* edepcuComp)
{
  collapseEdep_CU(edepcu, edepcuComp);
  //Locate all intersections
  int raynum = blockDim.x*blockIdx.x+threadIdx.x;
  if(raynum >= nrays)
  {
    return;
  }
  //for each crossing of raynum
  for(int i = 0; i < ncrossings; i++)
  {
    //spatial coordinates
    int rx = vec4D(boxes_cu, beam,raynum,i,0, nrays, ncrossings, 2);
    int rz = vec4D(boxes_cu, beam,raynum,i,0, nrays, ncrossings, 2);
    int cnt = 0;
    if(!rx || !rz)
      break;
    rx--;
    rz--;
    for(int n = 0; n < nbeams;n++)
    {
      for(int j = 0; j < nrays; j++)
      {
        if(vec4D(marked_cu, n,j,rx,rz))
        {
          vec4D(ints, beam, raynum, i, cnt++, nrays, ncrossings, numstored);
        }
        if(cnt >= numstored)
        {
          break;
        }
      }
    }
  }
}
