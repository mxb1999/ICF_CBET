#include "include/implSim.h"
#include "include/customMath.h"
using namespace std;
//initializing necessary arrays for the calculation
void rayLaunch(double x_init, double z_init, double kx_init, double kz_init, double urayinit, int raynum)
{
  //Launch_Ray_XZ Array Declaration
  int thisx = 0;
  int thisz = 0;
  int thisx_0 = 0;
  int thisz_0 = 0;
  int thisx_00 = 0;
  int thisz_00 = 0;
  double* uray = new double[nt]{1.0};
  double* mytime = new double[nt]{0.0};
  span(mytime, dt, nt*dt, nt);
  double* amplitude_norm= new double[nt]{0.0};
  double* markingx = new double[nt]{0.0};
  double* markingz = new double[nt]{0.0};
  double* myx = new double[nt]{0.0};
  double* myz = new double[nt]{0.0};
  double* myvx = new double[nt]{0.0};
  double* myvz = new double[nt]{0.0};
  double* mykx = new double[nt]{0.0};
  double* mykz = new double[nt]{0.0};
  double* nuei = new double[nt]{0.0};

  //Initializing Arrays
  for(int i = 0; i < nt; i++)
  {
    nuei[i] = 1.0;
  }
  uray[0] = urayinit;

  myx[0] = x_init;

  myz[0] = z_init;
  //determining the initial x values within the desired range to track the beam
  for(int i = 0;i < nx;i++)
  {
    if(myx[0] - x[i] <= ((0.5+1.0e-10)*dx + 1e-11) && myx[0] - x[i] >= -1*((0.5+1.0e-10)*dx  + 1e-11) )
  {
      thisx_0=i;
      thisx_00=i;
      break;
    }
  }
  //determining the initial z values within the desired range to track the beam

  for(int i = 0;i < nz;i++)
  {
    if(myz[0] - z[i] <= ((0.5+1.0e-10)*dz + 1e-11) && myz[0] - z[i] >= -1*((0.5+1.0e-10)*dz + 1e-11) )
    {

      thisz_0=i;
      thisz_00=i;
      break;
    }
  }
  double k = sqrt((pow(omega,2.0)-pow(wpe[thisx_0][thisz_0],2.0))/pow(c,2.0));
  double knorm = sqrt(pow(kx_init,2.0)+pow(kz_init,2.0));
  mykx[0]=(kx_init/knorm)*k;			// Normalized value for the ray's initial k_x
  mykz[0]=(kz_init/knorm)*k;			// Normalized value for the ray's initial k_z
  myvx[0] = pow(c,2.0)*mykx[0]/omega;                   // v_group, group velocity (dw/dk) from D(k,w).
  myvz[0] =  pow(c,2.0)*mykz[0]/omega;
  markingx[0] = thisx_0;
  markingz[0] = thisz_0;


  //__________Time Stepping__________
    int numcrossing = 0;
    //looping through time intervals
    for(int i = 1; i < nt;i++)
    {
      myvz[i] = myvz[i-1] - pow(c,2.0)/(2.0*ncrit)*dedendz[thisx_0][thisz_0]*dt;
      myvx[i] = myvx[i-1] - pow(c,2.0)/(2.0*ncrit)*dedendx[thisx_0][thisz_0]*dt;
      myx[i] = myx[i-1] + myvx[i]*dt;
      myz[i] = myz[i-1] + myvz[i]*dt;

      int search_index_x = 1;
      int search_index_z = 1;
      int thisx_m = fmax(0, thisx_0-search_index_x);
      int thisx_p = fmin(nx-1, thisx_0+search_index_x);
      int thisz_m = fmax(0, thisz_0-search_index_z);
      int thisz_p = fmin(nz-1, thisz_0+search_index_z);

      //assigning the current x value to be tracked

      for(int j = thisx_m; j <= thisx_p;j++)
      {

        if ( myx[i] - x[j] <= ((0.5+1.0e-10)*dx + 1e-12) && myx[i] - x[j] >= -1*((0.5+1.0e-10)*dx + 1e-12))
        {
          thisx = j;
          break;
        }
      }

      //assigning the current z value to be tracked
      for(int j = thisz_m; j <= thisz_p; j++)
      {
        if (myz[i] - z[j] <= ((0.5+1.0e-10)*dz + 1e-12) && myz[i] - z[j] >= -1*((0.5+1.0e-10)*dz + 1e-12))
        {

          thisz = j;
          break;
         }
      }

      double linez[2]={myz[i-1], myz[i]};
      double linex[2]={myx[i-1], myx[i]};
      int lastx = 10000;
      int lastz = 10000;
      //iterating through the selected portions of the x spatial tracking arrays
      //Boxes stores the spatial locations of each crossing of each ray
      //Marked = trajectory of a single ray, Boxes = coordinates of each ray intersection
      for(int j = thisx_m; j <= thisx_p;j++)
      {
        double currx = x[j]-dx/2;
        if((myx[i] > currx && myx[i-1] < (currx + 1e-12)) || (myx[i] < currx && myx[i-1] > (currx- 1e-12)))
        {
          double m = (myz[i] - myz[i-1])/(myx[i]-myx[i-1]);
          double b = myz[i] - myx[i]*m;
          double crossx = m*currx+b;
          if(abs(crossx-lastz)>1.0e-20)
          {
            ints[beam][raynum][numcrossing]=uray[i];
            crossesx[beam][raynum][numcrossing] = currx;
            crossesz[beam][raynum][numcrossing] = crossx;
            if(myx[i] < (xmax+dx/2 + 1e-11) && myx[i] > (xmin-dx/2 - 1e-12))
            {
              boxes[beam][raynum][numcrossing][0] = thisx+1;
              boxes[beam][raynum][numcrossing][1] = thisz+1;
            }
            lastx = currx;
            numcrossing += 1;
            break;
          }
        }
      }

      //iterating through the selected portions of the z spatial tracking arrays
        for(int j = thisz_m; j <= thisz_p;j++)
        {
          double currz = z[j]-dz/2;


          if((myz[i] > (currz) && myz[i-1] < (currz + 1e-11)) || (myz[i] < (currz+1e-11) && myz[i-1] > (currz - 1e-11)))
          {
            double m = (myx[i] - myx[i-1])/(myz[i]-myz[i-1]);
            double b = myx[i] - myz[i]*m;
            double crossz = m*currz+b;
            if(abs(crossz-lastx)>1.0e-20)
            {
              ints[beam][raynum][numcrossing]=uray[i];
              crossesz[beam][raynum][numcrossing] = currz;
              crossesx[beam][raynum][numcrossing] = crossz;
              if(myz[i] < (zmax+dz/2 +1e-11) && myz[i] > (zmin-dz/2-1e-11))
              {
                cout << numcrossing << endl;
                boxes[beam][raynum][numcrossing][0] = thisx + 1;
                boxes[beam][raynum][numcrossing][1] = thisz + 1;
              }
              lastz = currz;
              numcrossing += 1;
              break;
            }
          }
        }
        thisx_0 = thisx;
        thisz_0 = thisz;
        markingx[i] = thisx;
        markingz[i] = thisz;
        //ensuring that the beams are not identical
        double slope;
        double ztarg;
        double xtarg;
        if(markingx[i] != markingx[i-1] && markingz[i] != markingz[i-1])
        {
          ztarg = z[thisz] - (dz/2.0);
          if(myvz[i] < 0.0)
          {
            ztarg = z[thisz] + (dz/2.0);
          }
          slope = (myz[i] - myz[i-1])/(myx[i] - myx[i-1]+1.0e-10);

          xtarg = x[thisx]-(dx/2.0);
          if(myvx[i] >= 0.0)
          {
            xtarg = x[thisx]+(dx/2.0);
          }
          slope = (myx[i] - myx[i-1])/(myz[i] - myz[i-1]+1.0e-10);
          ztarg = myz[i-1]+(xtarg-myx[i-1])/slope;

          for(int j = 0; j < nrays;j++)
          {
            if(marked[thisx*nz + thisz][j*nbeams + beam] == 0 || marked[thisx*nz + thisz][j*nbeams + beam] > (raynum + 1))
            {
              #pragma omp atomic write
              marked[thisx*nz + thisz][j*nbeams + beam] = raynum+1;
              #pragma omp atomic update
              present[beam][thisx][thisz] += 1.0;
              break;
            }
          }

        }else if(markingz[i] != markingz[i-1] )
        {
          ztarg = z[thisz]-dz/2.0;
          if(myvz[i] < 0.0)
          {
            ztarg = z[thisz]+dz/2.0;
          }
          slope = (myz[i] - myz[i-1])/(myx[i] - myx[i-1]+1.0e-10);
          xtarg = myx[i]+(ztarg-myz[i-1])/slope;
          for(int j = 0; j < nrays;j++)
          {
            if(marked[thisx*nz + thisz][j*nbeams + beam] == 0 || marked[thisx*nz + thisz][j*nbeams + beam] > (raynum + 1))
            {
              #pragma omp atomic write
              marked[thisx*nz + thisz][j*nbeams + beam] = raynum+1;
              #pragma omp atomic update
              present[beam][thisx][thisz] += 1.0;
              break;
            }
          }
        }else if (markingx[i] != markingx[i-1])
        {
          xtarg = x[thisx]-dx/2.0;
          if(myvx[i] < 0.0)
          {
            xtarg = x[thisx]+dx/2.0;
          }
          slope = (myx[i] - myx[i-1])/(myz[i] - myz[i-1]+1.0e-10);
          ztarg = myz[i]+(xtarg-myx[i-1])/slope;

          for(int j = 0; j < nrays;j++)
          {
            if(marked[thisx*nz + thisz][j*nbeams + beam] == 0 || marked[thisx*nz + thisz][j*nbeams + beam] > (raynum + 1))
            {
              #pragma omp atomic write
              marked[thisx*nz + thisz][j*nbeams + beam] = raynum+1;
              #pragma omp atomic update
              present[beam][thisx][thisz] += 1.0;
              break;
            }
          }
        }



        //altering energy deposition values due to the launched rays
        uray[i] = uray[i-1];
  	    double increment = uray[i];
        double xp = (myx[i] - (x[thisx]+dx/2.0))/dx;
        double zp = (myz[i] - (z[thisz]+dz/2.0))/dz;
        if ( xp >= 0 && zp >= 0 ){
        double dl = zp;
        double dm = xp;
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x+1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z+1)
        double a4 = dl*dm;			// red 		: (x+1, z+1)
        #pragma omp atomic update
        edep[beam][thisx+1][thisz+1] += a1*increment;	// blue
        #pragma omp atomic update
        edep[beam][thisx+1+1][thisz+1] += a2*increment;	// green
        #pragma omp atomic update
        edep[beam][thisx+1][thisz+1+1] += a3*increment;	// yellow
        #pragma omp atomic update
        edep[beam][thisx+1+1][thisz+1+1] += a4*increment;	// red
      } else if ( xp < 0 && zp >= 0 ){
        double dl = zp;
        double dm = abs(xp);		// because xp < 0
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x-1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z+1)
        double a4 = dl*dm;			// red 		: (x-1, z+1)
        #pragma omp atomic update
        edep[beam][thisx+1][thisz+1] += a1*increment;	// blue
        #pragma omp atomic update
        edep[beam][thisx-1+1][thisz+1] += a2*increment;	// green
        #pragma omp atomic update
        edep[beam][thisx-1+1][thisz+1+1] += a4*increment;	// red
        #pragma omp atomic update
        edep[beam][thisx+1][thisz+1+1] += a3*increment;	// yellow

      } else if ( xp >= 0 && zp < 0 ){
        double dl = abs(zp);		// because zp < 0
        double dm = xp;
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x+1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z-1)
        double a4 = dl*dm;			// red 		: (x+1, z-1)
        #pragma omp atomic update
        edep[beam][thisx+1][thisz-1+1] += a3*increment;	// yellow
        #pragma omp atomic update
        edep[beam][thisx+1+1][thisz-1+1] += a4*increment;	// red
        #pragma omp atomic update
        edep[beam][thisx+1][thisz+1] += a1*increment;	// blue
        #pragma omp atomic update
        edep[beam][thisx+1+1][thisz+1] += a2*increment;	// green
      } else if ( xp < 0 && zp < 0 ){
        double dl = abs(zp);		// because zp < 0
        double dm = abs(xp);		// because xp < 0
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x-1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z-1)
        double a4 = dl*dm;			// red 		: (x-1, z-1)
        #pragma omp atomic update
        edep[beam][thisx+1][thisz+1] += a1*increment;	// blue
        #pragma omp atomic update
        edep[beam][thisx+1][thisz+1-1] += a3*increment;	// yellow
        #pragma omp atomic update
        edep[beam][thisx+1-1][thisz+1] += a2*increment;	// green
        #pragma omp atomic update
        edep[beam][thisx+1-1][thisz+1-1] += a4*increment;	// red
      } else {
        double store = edep[0][thisx][thisz];
        #pragma omp atomic write
        edep[0][thisx][thisz] = store + (nuei[i] * (*(eden[thisx]+thisz))/ncrit * uray[i-1]*dt);
        cout << "***** ERROR in interpolation of laser deposition to grid!! *****" << endl;
        break;
      }
      #pragma omp atomic write
      amplitude_norm[i] = (pow(omega,2.0)-pow(*(wpe[thisx_00]+thisz_00),2.0))/(pow(omega,2.0)-pow(pow(*(wpe[thisx]+thisz),2.0),(1./4.)));
      #pragma omp atomic write
      mytime[i] = dt*i;
      if ( (myx[i] < (xmin-(dx/2.0))) || (myx[i] > (xmax+(dx/2.0))))
      {
        break;                  // "breaks" out of the i loop once the if condition is satisfied
      } else if ( (myz[i] < (zmin-(dz/2.0))) || (myz[i] > (zmax+(dz/2.0)))){
           // the "|" means "or" (symbol above the return key)
        break;
    }
  }
  delete [] mytime;
  delete [] nuei;
  delete [] amplitude_norm;
  delete [] markingx;
  delete [] markingz;
  delete [] uray;
  delete [] myx;
  delete [] myz;
  delete [] mykx;
  delete [] mykz;
  delete [] myvx;
  delete [] myvz;

}






//use two threads here
void launch_ray_XZ(double x_init, double z_init, double kx_init, double kz_init,double urayinit, int raynum)
{
  rayLaunch(x_init,z_init, kx_init, kz_init, urayinit, raynum);

}
