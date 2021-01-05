#include "implSim.hpp"
#include "customMath.hpp"
using namespace std;

void track(int raynum, int xinit, int zinit, int kxinit, int kzinit, int** count, int beam, double** force, double urayinit)
{

  double xp = xinit;//initial x position
  double zp = zinit;//initial z position
  int currX = (int)((xp-xmin)*nx/(xmax-xmin));//initial x cell location
  int currZ = (int)((zp-zmin)*nz/(zmax-zmin));//initial z cell location
  //initialize ray variables
  double wpeInit = sqrt(eden[currX][currZ]*1e6*pow(ec,2.0)/(me*e0));
  double k = sqrt((pow(omega,2.0)-pow(wpeInit,2.0))/pow(c,2.0));
  double knorm = sqrt(pow(kxinit,2.0)+pow(kzinit,2.0));
  double kx = k*kxinit/knorm;
  double kz = k*kzinit/knorm;
  double rayV[] = {kx*pow(c,2)/omega, kz*pow(c,2)/omega};//ray velocity

  //storage variables to enable intelligent updates
  int prevX = 0;
  int prevZ = 0;
  double deltaX;
  double deltaZ;
  int numcrossing = 0;
  while(abs(xp) <= xmax && abs(zp) <= zmax)
  {
    //If the ray enters a new cell
    int cond = (prevX != currX || prevZ != currZ);
    if(prevX != currX || prevZ != currZ)
    {
      int delx = (prevX < currX);
      int delz = (prevZ < currZ);
      deltaX = xp-dx*currX;
      deltaZ = zp-dz*currZ;
      boxes[beam][raynum][numcrossing][0] = currX + 1;
      boxes[beam][raynum][numcrossing][1] = currZ + 1;
      marked[(currX*nz + currZ)*nbeams + beam].push(raynum+1);
      crossesx[beam][raynum][numcrossing] = xp;    //If the ray spends multiple iterations in the same cell
      crossesx[beam][raynum][numcrossing] = zp;
    }else
    {

    }
    rayV[0] -= force[0][currX*nz + currZ]*dt;
    rayV[1] -= force[1][currX*nz + currZ]*dt;
    deltaX = rayV[0]*dt;
    deltaZ = rayV[1]*dt;
    //Update position variables for comparison in the next iteration
    prevX = currX;
    prevZ = currZ;
    xp += rayV[0]*dt;
    zp += rayV[1]*dt;
    numcrossing++;
    currX = (int)((xp-xmin)*nx/(xmax-xmin));
    currZ = (int)((zp-zmin)*nz/(zmax-zmin));
    //if(raynum > nrays)
      //printf("XF <%e,%e>, %d\n", rayV[0], rayV[1], raynum);

  }
}


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

  //determining the initial x grid index within the desired range to track the beam
  for(int i = 0;i < nx;i++)
  {
    if(myx[0] - x[i] <= ((0.5+1.0e-10)*dx + 1e-11) && myx[0] - x[i] >= -1*((0.5+1.0e-10)*dx  + 1e-11) )
  {
      thisx_0=i;
      thisx_00=i;
      break;
    }
  }

  //determining the initial z grid index within the desired range to track the beam
  for(int i = 0;i < nz;i++)
  {
    if(myz[0] - z[i] <= ((0.5+1.0e-10)*dz + 1e-11) && myz[0] - z[i] >= -1*((0.5+1.0e-10)*dz + 1e-11) )
    {

      thisz_0=i;
      thisz_00=i;
      break;
    }
  }
  //determining the velocity characteristics of the ray based upon its initial position
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
      int thisx_m = fmax(0, thisx_0-search_index_x );
      int thisx_p = fmin(nx-1, thisx_0+search_index_x);
      int thisz_m = fmax(0, thisz_0-search_index_z);
      int thisz_p = fmin(nz-1, thisz_0+search_index_z);
      //determining the current x index of the ray
      for(int j = thisx_m; j <= thisx_p;j++)
      {
        if ( myx[i] - x[j] <= ((0.5+1.0e-10)*dx + 1e-12) && myx[i] - x[j] >= -1*((0.5+1.0e-10)*dx + 1e-12))
        {
          thisx = j;
          break;
        }
      }

      //determining the current z index of the ray
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
        //if the ray is currently between within the desired caustic zone for a crossing
        if((myx[i] > currx && myx[i-1] <= (currx + 1e-10)) || (myx[i] < currx && myx[i-1] >= (currx- 1e-10)))
        {
          double m = (myz[i] - myz[i-1])/(myx[i]-myx[i-1]);
          double b = myz[i] - myx[i]*m;
          double crossx = m*currx+b;
          //if the ray has moved since last update
          if(abs(crossx-lastz)>1.0e-20)
          {
            ints[beam][raynum][numcrossing]=uray[i];
            crossesx[beam][raynum][numcrossing] = currx;
            crossesz[beam][raynum][numcrossing] = crossx;
            //if ray is still within the grid
            if(myx[i] < (xmax+dx/2 + 1e-10) && myx[i] > (xmin-dx/2 - 1e-10))
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
      //Same idea as previous loop, but for the Z coordinate instead of x
        for(int j = thisz_m; j <= thisz_p;j++)
        {
          double currz = z[j]-dz/2;
          if((myz[i] > (currz) && myz[i-1] < (currz + 1e-10)) || (myz[i] < (currz) && myz[i-1] > (currz - 1e-10)))
          {

            double m = (myx[i] - myx[i-1])/(myz[i]-myz[i-1]);
            double b = myx[i] - myz[i]*m;
            double crossz = m*currz+b;
            if(abs(crossz-lastx) > 1.0e-20)
            {
              ints[beam][raynum][numcrossing]=uray[i];
              crossesz[beam][raynum][numcrossing] = currz;
              crossesx[beam][raynum][numcrossing] = crossz;
              if(myz[i] < (zmax+dz/2 +1e-10) && myz[i] > (zmin-dz/2-1e-10))
              {
                boxes[beam][raynum][numcrossing][0] = thisx + 1;
                boxes[beam][raynum][numcrossing][1] = thisz + 1;
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
        markingx[i] = thisx;
        markingz[i] = thisz;
        if(markingx[i] != markingx[i-1] && markingz[i] != markingz[i-1]) //ensure that the ray has left the previous grid zone during the last time step
        {
          //store that this ray has crossed this zone
          marked[(thisx*nz + thisz)*nbeams+beam].push(raynum+1);//this ray has crossed this zone from this beam
          #pragma omp atomic update
          present[beam][thisx][thisz] += 1.0;//another ray has crossed this zone from this beam
        }else if(markingz[i] != markingz[i-1])//ensure that the ray has left the previous grid zone during the last time step
        {
          //store that this ray has crossed this zone
          marked[(thisx*nz + thisz)*nbeams+beam].push(raynum+1);//this ray has crossed this zone from this beam
          #pragma omp atomic update
          present[beam][thisx][thisz] += 1.0;//another ray has crossed this zone from this beam
        }else if (markingx[i] != markingx[i-1])//ensure that the ray has left the previous grid zone during the last time step
        {
          //store that this ray has crossed this zone
          marked[(thisx*nz + thisz)*nbeams+beam].push(raynum+1);//this ray has crossed this zone from this beam
          #pragma omp atomic update
          present[beam][thisx][thisz] += 1.0;//another ray has crossed this zone from this beam
        }


        printf("%f %f\n",uray[i], uray[i-1]);
        //Deposit energy due to the incident ray
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
