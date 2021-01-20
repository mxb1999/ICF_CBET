#include <stdio.h>
#include <stdlib.h>
#include "Trace_interface.h"
#include "io_interface.h"
__device__
struct rayinit//composite struct, compactly store necessary initial information
{
    double xinit;
    double zinit;
    double kxinit;
    double kzinit;
    double beam;
    double urayinit;
    double edeninit;
};
__global__ void track(rayinit* init, double* deposit, double* force, double* dkx, double* dkz)
{
    int raynum = blockDim.x*blockIdx.x+threadIdx.x;
    
    double xloc = xinit[index];//initial x position
    double zloc = zinit[index];//initial z position

    int currX = (int)((xloc-xmin)/(dx));//initial x cell location
    int currZ = (int)((zloc-zmin)/(dz));//initial z cell location
    int relx = ((xloc-xmin)-dx*currX > dx/2);
    int relz = ((zloc-zmin)-dz*currZ > dz/2);
    currX += relx;
    currZ += relz;
    //printf("(%d, %d)\n", currX, currZ);
  
    //initialize ray variables
    double wpeInit = sqrt(eden[currX][currZ]*1e6*pow(ec,2.0)/(me*e0));
    double k = sqrt((pow(omega,2.0)-pow(wpeInit,2.0))/pow(c,2.0));
    double knorm = sqrt(pow(kxinit,2.0)+pow(kzinit,2.0));
    double kx = k*kxinit/knorm;
    double kz = k*kzinit/knorm;
    double rayV[] = {kx*pow(c,2)/omega, kz*pow(c,2)/omega};//ray velocity
    //storage variables to enable intelligent updates
    int prevX = -1;
    int prevZ = -1;
    double deltaX = 0;
    double deltaZ = 0;
    int numcrossing = 0;
    double mag;
    while((currX < nx && currZ < nz) && (currX >= 0 && currZ >= 0))
    {
      
      //If the ray enters a new cell
      int cond = (prevX != currX || prevZ != currZ);
      if(prevX != currX || prevZ != currZ)
      {
        //printf("XCross %d, ZCross %d:(%f, %f), ([%f,%f], [%f, %f]) \n",(prevX != currX), prevZ != currZ, xloc, zloc, dx*currX+ xmin, dx*currX+dx + xmin, dz*currZ+ zmin, dz*currZ+dz+ zmin);
        //deltaX -= (abs(xloc-xmin)-currX*dx);//remove overshoot distance
       // deltaZ -= (abs(zloc-zmin)-currZ*dz);
        //update boxes at each crossing into a new zone
        int* boxx = vec4D(boxes, beam,raynum,numcrossing,0, nrays, ncrossings, 2);//[beam][raynum][numcrossing][0]
        int* boxz = vec4D(boxes, beam,raynum,numcrossing,1, nrays, ncrossings, 2);//[beam][raynum][numcrossing][1]
        *boxx = currX + 1;
        *boxz = currZ + 1;//[beam][raynum][numcrossing][1]
        //printf("%d\n",marked[(currX*nz + currZ)*nbeams + beam].size());
        dkx[(beam*nrays+raynum)*ncrossings+numcrossing] = deltaX;//calculate dkmag later, reduce main loop bloat
        dkz[beam][raynum][numcrossing] = deltaZ;
        numcrossing++;
        if(raynum == 0 && beam == 1)
        {
          raypath[currX][currZ]++;
        }
        deltaX = 0;
        deltaZ = 0;
      }
      //printf("(%d, %d)\n",currX, currZ);
      rayV[0] -= pow(c,2.0)/(2.0*ncrit)*dedendx[currX][currZ]*dt;
      rayV[1] -= pow(c,2.0)/(2.0*ncrit)*dedendz[currX][currZ]*dt;
      deltaX += rayV[0]*dt;
      deltaZ += rayV[1]*dt;
      //Update position variables for comparison in the next iteration
      prevX = currX;
      prevZ = currZ;
      xloc += rayV[0]*dt;
      zloc += rayV[1]*dt;
      //if(raynum > nrays)
        //printf("XF <%e,%e>, %d\n", rayV[0], rayV[1], raynum);
      double xp = (xloc - (x[currX]+dx/2.0))/dx;
      double zp = (zloc - (z[currZ]+dz/2.0))/dz;
      //uray = urayinit;
      double increment = urayinit;
  
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
  
      #pragma omp atomic update
      edep[beam][currX+1][currZ+1] += a1*increment;	// blue
      #pragma omp atomic update
      edep[beam][currX+xadd+1][currZ+1] += a2*increment;	// green
      #pragma omp atomic update
      edep[beam][currX+1][currZ+zadd+1] += a3*increment;	// yellow
      #pragma omp atomic update
      edep[beam][currX+xadd+1][currZ+zadd+1] += a4*increment;	// red
  
      double xabs = xloc-xmin;
      double zabs = zloc-zmin;
      currX = (int)(xabs/dx);
      currZ = (int)(zabs/dz);
      relx = (xabs-dx*currX > dx/2);
      relz = (zabs-dz*currZ > dz/2);
      currX += relx;
      currZ += relz;
  
    }
}


int cudaRayTrace(gconfig* gpu)
{

};
