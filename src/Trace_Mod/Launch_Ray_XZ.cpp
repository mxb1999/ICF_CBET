#include "Trace_interface.hpp"
using namespace std;

#define MAG(x,z) sqrt(pow(x,2) + pow(z,2))
void track(int raynum, double xinit, double zinit, double kx_init, double kz_init, double urayinit, int beam)
{

  double x_coor = xinit;//initial x position
  double z_coor = zinit;//initial z position

  int currX = (int)((x_coor-xmin)/(dx));//initial x cell location
  int currZ = (int)((z_coor-zmin)/(dz));//initial z cell location

  //initialize ray variables
  double wpeInit = sqrt(vec2D(eden,currX,currZ,nz)*1e6*pow(ec,2.0)/(me*e0));
  double k = sqrt((pow(omega,2.0)-pow(wpeInit,2.0))/pow(c,2.0));
  double knorm = sqrt(pow(kx_init,2.0)+pow(kz_init,2.0));
  double kx = k*kx_init/knorm;
  double kz = k*kz_init/knorm;
  double rayV[] = {kx*pow(c,2)/omega, kz*pow(c,2)/omega};//ray velocity
  //storage variables to enable intelligent updates
  double x_coor_prev = x_coor;
  double z_coor_prev = z_coor;
  int prevX = -1;
  int prevZ = -1;
  double deltaX = 0;
  double deltaZ = 0;
  int numcrossing = 0;
  double mag;
  double prevcross[2] = {currX * dx + xmin - dx/2, currZ*dz + zmin- dz/2};
  short dirX = (currX != 0);//1 if right, 0 if left
  short dirZ = (currZ != 0);//1 if top, 0 if bottom
  while((currX < nx && currZ < nz) && (currX >= 0 && currZ >= 0))
  {

    //If the ray enters a new cell
    short condx = (prevX != currX);//indicate which was changed
    short condz = (prevZ != currZ);

    if(condx || condz)
    {
      dirX = (prevX > currX);//indicate direction of propogation
      dirZ = (prevZ > currZ);
      double xnew;
      double znew;
      //set crossing coordinates depending on the side of the zone entered
      if(condx && condz)//both boundaries were crossed
      {
        if(raynum == 0 && beam == 1)
          printf("\t1\n");
        xnew = currX*dx+xmin+dx*dirX;//added dx*dirX term differentiates whether the ray entered in the top/bottom and left/right of the zone
        znew = currZ*dz+zmin+dz*dirZ;
      }else if(condx)//the x boundary was crossed, linearly interpolate x
      {
        if(raynum == 0 && beam == 1)
          printf("\t2\n");
        xnew = currX*dx+xmin + dx*dirX;
        double m = (z_coor-z_coor_prev)/(x_coor-x_coor_prev);
        double b = 0.5*(z_coor+z_coor_prev-m*(x_coor_prev + x_coor));
        znew = m*xnew+b;
      }else//the z boundary was crossed, linearly interpolate x
      {
        
        znew = currZ*dz+zmin+dz*dirZ;
        double m = (z_coor-z_coor_prev)/(x_coor-x_coor_prev);
        double b = 0.5*(z_coor+z_coor_prev-m*(x_coor_prev + x_coor));
        xnew = (znew-b)/m;
      }
      //printf("XCross %d, ZCross %d:(%f, %f), ([%f,%f], [%f, %f]) \n",(prevX != currX), prevZ != currZ, x_coor, z_coor, dx*currX+ xmin, dx*currX+dx + xmin, dz*currZ+ zmin, dz*currZ+dz+ zmin);
      //update boxes at each crossing into a new zone
      vec4DW(boxes, beam,raynum,numcrossing,0, nrays, ncrossings, 2, currX+1);//[beam,raynum,numcrossing,0]
      vec4DW(boxes, beam,raynum,numcrossing,1, nrays, ncrossings, 2,currZ+1);//[beam,raynum,numcrossing,1]
      //printf("%d\n",marked[(currX*nz + currZ)*nbeams + beam].size());
      vec2DW(crossesx,beam,raynum,numcrossing,xnew);    //If the ray spends multiple iterations in the same cell
      crossesz[beam,raynum,numcrossing] = znew;

      double xchange = xnew - prevcross[0];
      double zchange = znew - prevcross[1];
      dkx[beam,raynum,numcrossing] = xchange;
      dkz[beam,raynum,numcrossing] = zchange;
      dkmag[beam,raynum,numcrossing] = MAG(xchange, zchange);
      prevcross[0] = xnew;
      prevcross[1] = znew;
      numcrossing++;

      deltaX = 0;
      deltaZ = 0;
    }
    //printf("(%d, %d)\n",currX, currZ);
    rayV[0] -= pow(c,2.0)/(2.0*ncrit)*dedendx[currX,currZ]*dt;
    rayV[1] -= pow(c,2.0)/(2.0*ncrit)*dedendz[currX,currZ]*dt;
    //printf("%e\n",rayV[0]*dt);
    deltaX += rayV[0]*dt;
    deltaZ += rayV[1]*dt;
    //Update position variables for comparison in the next iteration
    prevX = currX;
    prevZ = currZ;
    x_coor_prev = x_coor;
    x_coor_prev = z_coor;
    x_coor += rayV[0]*dt;
    z_coor += rayV[1]*dt;
    //if(raynum > nrays)
      //printf("XF <%e,%e>, %d\n", rayV[0], rayV[1], raynum);
    double xp = (x_coor - (x[currX] + dx*dirX))/dx;
    double zp = (z_coor - (z[currZ] + dz*dirZ))/dz;
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
    double xabs = x_coor-xmin;
    double zabs = z_coor-zmin;
    currX = (int)((xabs)/dx);
    currZ = (int)((zabs)/dz);
    
  }
}

//initializing necessary arrays for the calculation
void rayLaunch(double x_init, double z_init, double kx_init, double kz_init, double urayinit, int raynum, int beam)
{
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
  double myx = x_init;
  double myz = z_init;
  double myvx;
  double myvz;
  double mykx;
  double mykz;
  double myvxprev;
  double myvzprev;
  //determining the initial x grid index within the desired range to track the beam


  //determining the initial z grid index within the desired range to track the beam
  for(int i = 0;i < nz;i++)
  {
    if(z_init - z[i] <= ((0.5+1.0e-10)*dz + 1e-11) && z_init - z[i] >= -1*((0.5+1.0e-10)*dz + 1e-11) )
    {

      thisz_0=i;
      break;
    }
  }

  for(int i = 0;i < nx;i++)
  {

    double initX = x_init;
    double currX = x[i];
    if(initX - currX <= ((0.5+1.0e-10)*dx + 1e-11) && initX - currX >= -1*((0.5+1.0e-10)*dx + 1e-11) )
    {
      thisx_0=i;
      break;
    }
  }
  
  //determining the velocity characteristics of the ray based upon its initial position
  double k = sqrt((pow(omega,2.0)-pow(vec2D(wpe,thisx_0, thisz_0, nz),2.0))/pow(c,2.0));
  double knorm = sqrt(pow(kx_init,2.0)+pow(kz_init,2.0));
  mykx=(kx_init/knorm)*k;			// Normalized value for the ray's initial k_x
  mykz=(kz_init/knorm)*k;			// Normalized value for the ray's initial k_z
  myvxprev = pow(c,2.0)*mykx/omega;                   // v_group, group velocity (dw/dk) from D(k,w).
  myvzprev =  pow(c,2.0)*mykz/omega;
  ////markingxprev = thisx_0;
  ////markingzprev = thisz_0;
  //__________Time Stepping__________
    int numcrossing = 0;

    //looping through time intervals
    for(int i = 1; i < nt;i++)
    {
      
      double forcex = pow(c,2.0)/(2.0*ncrit)*vec2D(dedendx,thisx_0,thisz_0, nz);
      double forcez = pow(c,2.0)/(2.0*ncrit)*vec2D(dedendz,thisx_0,thisz_0, nz);


      myvx = myvxprev - forcex*dt;
      myvz = myvzprev - forcez*dt;
      myx += myvx*dt;
      myz += myvz*dt;

      int search_index_x = 1;
      int search_index_z = 1;
      int thisx_m = fmax(0, thisx_0-search_index_x );
      int thisx_p = fmin(nx-1, thisx_0+search_index_x);
      int thisz_m = fmax(0, thisz_0-search_index_z);
      int thisz_p = fmin(nz-1, thisz_0+search_index_z);
      //determining the current x index of the ray
      for(int j = thisx_m; j <= thisx_p;j++)
      {
        if ( myx - x[j] <= ((0.5+1.0e-10)*dx + 1e-12) && myx - x[j] >= -1*((0.5+1.0e-10)*dx + 1e-12))
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
        if (myz - z[j] <= ((0.5+1.0e-10)*dz + 1e-12) && myz - z[j] >= -1*((0.5+1.0e-10)*dz + 1e-12))
        {

          thisz = j;
          break;
         }
      }
      //double linez[2]={z_init, myz};
      //double linex[2]={x_init, myx};
      int lastx = 10000;
      int lastz = 10000;
      //iterating through the selected portions of the x spatial tracking arrays
      //boxes stores the spatial locations of each crossing of each ray
      //Marked = trajectory of a single ray, boxes = coordinates of each ray intersection
      
      for(int j = thisx_m; j <= thisx_p;j++)
      {
        double currx = x[j]-dx/2;//crossing into 
        //if the ray is currently between within the desired caustic zone for a crossing
        if((myx > currx && x_init <= (currx + 1e-10)) || (myx < currx && x_init >= (currx- 1e-10)))
        {
          double m = (myz - z_init)/(myx-x_init);
          double b = myz - myx*m;
          double crossx = m*currx+b;
          //if the ray has moved since last update
          if(abs(crossx-lastz)>1.0e-20)
          {
            vec3DW(crossesx, beam, raynum, numcrossing, nrays, ncrossings,currx);
            vec3DW(crossesz, beam, raynum, numcrossing, nrays, ncrossings,crossx);
            //if ray is still within the grid
            if(myx < (xmax+dx/2 + 1e-10) && myx > (xmin-dx/2 - 1e-10))
            {
              //printf("(%f, %f), ([%f,%f], [%f, %f]) \n", myx, myz, dx*thisx+ xmin, dx*thisx+dx + xmin, dz*thisz+ zmin, dz*thisz+dz+ zmin);
              vec4DW(boxes, beam,raynum,numcrossing,0, nrays, ncrossings, 2, thisx+1);//[beam][raynum][numcrossing][0]
              vec4DW(boxes, beam,raynum,numcrossing,1, nrays, ncrossings, 2, thisz+1);//[beam][raynum][numcrossing][1]
              //vec4DW(marked, beam,raynum,thisx,thisz, nrays, nx, nz, 1);
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
          double currz = z[j]-dz/2;//center of the jth zone
          //printf("REE\n");
          if((myz > (currz) && z_init < (currz + 1e-10)) || (myz < (currz) && z_init > (currz - 1e-10)))//if myz is approximately equal to a zone crossing
          {

            double m = (myx - x_init)/(myz-z_init);
            double b = myx - myz*m;
            double crossz = m*currz+b;
            if(abs(crossz-lastx) > 1.0e-20)
            {
              vec3DW(crossesx, beam, raynum, numcrossing, nrays, ncrossings,crossz);
              vec3DW(crossesz, beam, raynum, numcrossing, nrays, ncrossings,currz);
              if(myz < (zmax+dz/2 +1e-10) && myz > (zmin-dz/2-1e-10))
              {
               vec4DW(boxes, beam,raynum,numcrossing,0, nrays, ncrossings, 2, thisx+1);//[beam][raynum][numcrossing][0]
               vec4DW(boxes, beam,raynum,numcrossing,1, nrays, ncrossings, 2, thisz+1);//[beam][raynum][numcrossing][1]
              // vec4DW(marked, beam,raynum,thisx,thisz, nrays, nx, nz, 1);
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
  	    /*double increment = thisInit.urayinit;
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

          break;
        }
        //printf("edep: %d %d %d %d %d\n",thisx+1, thisz+1, beam, raynum, RAYS);
        vec4DI(edep, beam, raynum, thisx+1, thisz+1, nrays, nx+2,nz+2, a1*increment);	// blue
        vec4DI(edep, beam, raynum, thisx+xadd+1, thisz+1, nrays, nx+2,nz+2, a2*increment);// green
        vec4DI(edep, beam, raynum, thisx+1, thisz+zadd+1, nrays, nx+2,nz+2, a3*increment);// yellow
        vec4DI(edep, beam, raynum, thisx+xadd+1, thisz+zadd+1, nrays, nx+2,nz+2, a4*increment);	// red

        */
        x_init = myx;
        z_init = myz;
        myvxprev = myvx;
        myvzprev = myvz;
        ////markingxprev = //markingx;
        ////markingzprev = //markingz;
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







//use two threads here
void launch_ray_XZ(rayinit initSettings, int raynum)
{
  double xinit = initSettings.xinit;
  double zinit = initSettings.zinit;
  double kx_init = initSettings.kxinit;
  double kz_init = initSettings.kzinit;
  int beam = initSettings.beam;
  double urayinit = initSettings.urayinit;    
  rayLaunch(xinit,zinit, kx_init, kz_init, urayinit, raynum,beam);
  
}
