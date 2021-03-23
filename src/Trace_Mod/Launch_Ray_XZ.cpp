#include "Trace_interface.hpp"

using namespace std;

#define MAG(x,z) sqrt(pow(x,2) + pow(z,2))
void track(int raynum, double xinit, double zinit, double kxinit, double kzinit, double urayinit, int beam)
{

  double x_coor = xinit;//initial x position
  double z_coor = zinit;//initial z position

  int currX = (int)((x_coor-xmin)/(dx));//initial x cell location
  int currZ = (int)((z_coor-zmin)/(dz));//initial z cell location

  //initialize ray variables
  double wpeInit = sqrt(vec2D(eden,currX,currZ,nz)*1e6*pow(ec,2.0)/(me*e0));
  double k = sqrt((pow(omega,2.0)-pow(wpeInit,2.0))/pow(c,2.0));
  double knorm = sqrt(pow(kxinit,2.0)+pow(kzinit,2.0));
  double kx = k*kxinit/knorm;
  double kz = k*kzinit/knorm;
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

    #pragma omp atomic update
    edep[beam,currX+1,currZ+1] += a1*increment;	// blue
    #pragma omp atomic update
    edep[beam,currX+xadd+1,currZ+1] += a2*increment;	// green
    #pragma omp atomic update
    edep[beam,currX+1,currZ+zadd+1] += a3*increment;	// yellow
    #pragma omp atomic update
    edep[beam,currX+xadd+1,currZ+zadd+1] += a4*increment;	// red

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
  int thisx_00 = 0;
  int thisz_00 = 0;
  double uray;
//  double* mytime = new double[nt]{0.0};
  //span(mytime, dt, nt*dt, nt);
  //double* amplitude_norm= new double[nt]{0.0};
  int markingxprev;// = new double[nt]{0.0};
  int markingzprev;// = new double[nt]{0.0};
  int markingx;// = new double[nt]{0.0};
  int markingz;
  double myx;
  double myz;
  double myxprev;
  double myzprev;
  double myvx;
  double myvz;
  double mykx;
  double mykz;
  double myvxprev;
  double myvzprev;
/*  double* nuei = new double[nt]{0.0};

  //Initializing Arrays
  for(int i = 0; i < nt; i++)
  {
    nuei[i] = 1.0;
  }*/
  myxprev = x_init;
  myzprev = z_init;

  //determining the initial x grid index within the desired range to track the beam
  for(int i = 0;i < nx;i++)
  {
    if(myxprev - x[i] <= ((0.5+1.0e-10)*dx + 1e-11) && myxprev - x[i] >= -1*((0.5+1.0e-10)*dx + 1e-11) )
    {
      thisx_0=i;
      thisx_00=i;
      break;
    }
  }

  //determining the initial z grid index within the desired range to track the beam
  for(int i = 0;i < nz;i++)
  {
    if(myzprev - z[i] <= ((0.5+1.0e-10)*dz + 1e-11) && myzprev - z[i] >= -1*((0.5+1.0e-10)*dz + 1e-11) )
    {

      thisz_0=i;
      thisz_00=i;
      break;
    }
  }
  //if(raynum == 0 && beam == 0)
  //{
 //   printf("Init: (%d, %d)\n", thisx_0, thisz_0);
 // }

  //determining the velocity characteristics of the ray based upon its initial position
  double k = sqrt((pow(omega,2.0)-pow(vec2D(wpe, thisx_0, thisz_0, nz),2.0))/pow(c,2.0));
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
      myvz = myvzprev - pow(c,2.0)/(2.0*ncrit)*vec2D(dedendz,thisx_0,thisz_0, nz)*dt;
      myvx = myvxprev - pow(c,2.0)/(2.0*ncrit)*vec2D(dedendx,thisx_0,thisz_0, nz)*dt;

      myx = myxprev + myvx*dt;
      myz = myzprev + myvz*dt;

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

      //determining the current z index of the ray
      for(int j = thisz_m; j <= thisz_p; j++)
      {
        if (myz - z[j] <= ((0.5+1.0e-10)*dz + 1e-12) && myz - z[j] >= -1*((0.5+1.0e-10)*dz + 1e-12))
        {

          thisz = j;
          break;
         }
      }

      if(raynum == 0 && beam == 1)
      {
        raypath[thisx*nz + thisz]++;
      }
      double linez[2]={myzprev, myz};
      double linex[2]={myxprev, myx};
      int lastx = 10000;//
      int lastz = 10000;
      //iterating through the selected portions of the x spatial tracking arrays
      //Boxes stores the spatial locations of each crossing of each ray
      //Marked = trajectory of a single ray, Boxes = coordinates of each ray intersection
      for(int j = thisx_m; j <= thisx_p;j++)
      {
        double currx = x[j]-dx/2;//-dx/2;//crossing into 
        //if the ray is currently between within the desired caustic zone for a crossing
        if((myx > currx && myxprev <= (currx + 1e-10)) || (myx < currx && myxprev >= (currx- 1e-10)))
        {
          double m = (myz - myzprev)/(myx-myxprev);
          double b = myz - myx*m;
          double crossx = m*currx+b;
          //if the ray has moved since last update
          if(abs(crossx-lastz)>1.0e-20)
          {
            ints[beam,raynum,numcrossing]=uray;
            crossesx[beam,raynum,numcrossing] = currx;
            crossesz[beam,raynum,numcrossing] = crossx;
            //if ray is still within the grid
            if(myx < (xmax+dx/2 + 1e-10) && myx > (xmin-dx/2 - 1e-10))
            {
              //printf("(%f, %f), ([%f,%f], [%f, %f]) \n", myx, myz, dx*thisx+ xmin, dx*thisx+dx + xmin, dz*thisz+ zmin, dz*thisz+dz+ zmin);
              vec4DW(boxes, beam,raynum,numcrossing,0, nrays, ncrossings, 2,thisx + 1);
              vec4DW(boxes, beam,raynum,numcrossing,1, nrays, ncrossings, 2,thisz + 1);
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
          double currz = z[j]-dz/2;//-dz/2;//center of the jth zone
          //printf("REE\n");
          if((myz > (currz) && myzprev < (currz + 1e-10)) || (myz < (currz) && myzprev > (currz - 1e-10)))//if myz is approximately equal to a zone crossing
          {

            double m = (myx - myxprev)/(myz-myzprev);
            double b = myx - myz*m;
            double crossz = m*currz+b;
            if(abs(crossz-lastx) > 1.0e-20)
            {
              ints[beam,raynum,numcrossing]=uray;
              crossesz[beam,raynum,numcrossing] = currz;
              crossesx[beam,raynum,numcrossing] = crossz;
              if(myz < (zmax+dz/2 +1e-10) && myz > (zmin-dz/2-1e-10))
              {
               // printf("Double ree\n");
                vec4DW(boxes, beam,raynum,numcrossing,0, nrays, ncrossings, 2,thisx + 1);//boxes[beam][raynum][numcrossing][0:1] = thisx+1:thisz+1; 
                vec4DW(boxes, beam,raynum,numcrossing,1, nrays, ncrossings, 2,thisz + 1);
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
        if(markingx != markingxprev && markingz != markingzprev) //ensure that the ray has left the previous grid zone during the last time step
        {
          //store that this ray has crossed this zone
          //marked[(thisx*nz + thisz)*nbeams+beam].push(raynum+1);//this ray has crossed this zone from this beam
        //  #pragma omp atomic update
         // present[beam,thisx,thisz] += 1.0;//another ray has crossed this zone from this beam
        }else if(markingz != markingzprev)//ensure that the ray has left the previous grid zone during the last time step
        {
          //store that this ray has crossed this zone
          //marked[(thisx*nz + thisz)*nbeams+beam].push(raynum+1);//this ray has crossed this zone from this beam
       //   #pragma omp atomic update
       //   present[beam,thisx,thisz] += 1.0;//another ray has crossed this zone from this beam
        }else if (markingx != markingxprev)//ensure that the ray has left the previous grid zone during the last time step
        {
          //store that this ray has crossed this zone
          //marked[(thisx*nz + thisz)*nbeams+beam].push(raynum+1);//this ray has crossed this zone from this beam
       //   #pragma omp atomic update
       //   present[beam,thisx,thisz] += 1.0;//another ray has crossed this zone from this beam
        }


        //Deposit energy due to the incident ray
        uray = urayinit;
  	    double increment = uray;
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
        if(a1*increment < 0 || a2*increment < 0 || a3*increment < 0 || a4*increment < 0  )
        {
          //printf("Kill me %e %e %d %d\n", vec2D(dedendx,thisx_0,thisz_0, nz),vec2D(dedendz,thisx_0,thisz_0, nz), thisx,thisz);
        }
        if((xadd != xadd) && (zadd != zadd))
        {
          printf("NaN Error: Error in Grid Interpolation \n");
          break;
        }
        
        vec3DIA(edep,beam,thisx+1,thisz+1,nx+2,nz+2,a1*increment);	// blue
        vec3DIA(edep,beam,thisx+1+xadd,thisz+1,nx+2,nz+2,a2*increment);	// blue
        vec3DIA(edep,beam,thisx+1,thisz+1+zadd,nx+2,nz+2,a3*increment);	// blue
        vec3DIA(edep,beam,thisx+1+xadd,thisz+1+zadd,nx+2,nz+2,a4*increment);	// blue

        
        
        myxprev = myx;
        myvxprev = myvx;
        myzprev = myz;
        myvzprev = myvz;
        markingxprev = markingx;
        markingzprev = markingz;
        urayinit = uray;
        /*if ( xp >= 0 && zp >= 0 ){

      } else if ( xp < 0 && zp >= 0 ){
        double dl = zp;
        double dm = abs(xp);		// because xp < 0
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x-1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z+1)
        double a4 = dl*dm;			// red 		: (x-1, z+1)
        #pragma omp atomic update
        edep[beam,thisx+1,thisz+1] += a1*increment;	// blue
        #pragma omp atomic update
        edep[beam,thisx-1+1,thisz+1] += a2*increment;	// green
        #pragma omp atomic update
        edep[beam,thisx-1+1,thisz+1+1] += a4*increment;	// red
        #pragma omp atomic update
        edep[beam,thisx+1,thisz+1+1] += a3*increment;	// yellow

      } else if ( xp >= 0 && zp < 0 ){
        double dl = abs(zp);		// because zp < 0
        double dm = xp;
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x+1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z-1)
        double a4 = dl*dm;			// red 		: (x+1, z-1)
        #pragma omp atomic update
        edep[beam,thisx+1,thisz+1] += a1*increment;	// blue
        #pragma omp atomic update
        edep[beam,thisx+1+1,thisz+1] += a2*increment;	// green
        #pragma omp atomic update
        edep[beam,thisx+1,thisz-1+1] += a3*increment;	// yellow
        #pragma omp atomic update
        edep[beam,thisx+1+1,thisz-1+1] += a4*increment;	// red
      } else if ( xp < 0 && zp < 0 ){
        double dl = abs(zp);		// because zp < 0
        double dm = abs(xp);		// because xp < 0
        double a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , z  )
        double a2 = (1.0-dl)*dm;		// green	: (x-1, z  )
        double a3 = dl*(1.0-dm);		// yellow	: (x  , z-1)
        double a4 = dl*dm;			// red 		: (x-1, z-1)
        #pragma omp atomic update
        edep[beam,thisx+1,thisz+1] += a1*increment;	// blue
        #pragma omp atomic update
        edep[beam,thisx+1-1,thisz+1] += a2*increment;	// green
        #pragma omp atomic update
        edep[beam,thisx+1,thisz+1-1] += a3*increment;	// yellow
        #pragma omp atomic update
        edep[beam,thisx+1-1,thisz+1-1] += a4*increment;	// red
      } else {
        double store = edep[0,thisx,thisz];
        #pragma omp atomic write
        edep[0,thisx,thisz] = store + (nuei[i] * (*(eden[thisx]+thisz))/ncrit * urayinit*dt);
        cout << "***** ERROR in interpolation of laser deposition to grid!! *****" << endl;
        break;
      }*/
    //  #pragma omp atomic write
    //  amplitude_norm[i] = (pow(omega,2.0)-pow(*(wpe[thisx_00]+thisz_00),2.0))/(pow(omega,2.0)-pow(pow(*(wpe[thisx]+thisz),2.0),(1./4.)));
    ///  #pragma omp atomic write
    //  mytime[i] = dt*i;
      printf("%f %f\n", myx, myz);
      if ( (myx < (xmin-(dx/2.0))) || (myx > (xmax+(dx/2.0))))
      {
        printf("Break X %f\n", myx);
        break;                  // "breaks" out of the i loop once the if condition is satisfied
      } else if ( (myz < (zmin-(dz/2.0))) || (myz > (zmax+(dz/2.0)))){
           // the "|" means "or" (symbol above the return key)
        printf("Break Z %f\n", myz);

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
  double kxinit = initSettings.kxinit;
  double kzinit = initSettings.kzinit;
  int beam = initSettings.beam;
  double urayinit = initSettings.urayinit;    
  if(switchvar)
  {
    track(raynum, xinit, zinit, kxinit, kzinit, urayinit, beam);
  }else
  {
    rayLaunch(xinit,zinit, kxinit, kzinit, urayinit, raynum,beam);
  }
  if(raynum == 0 && beam == 1)
  {
    for(int i = 0; i < ncrossings; i++)
    {
      int temp1 = vec4D(boxes,beam,raynum,i,0,nrays, ncrossings,2);
      int temp2 = vec4D(boxes,beam,raynum,i,1,nrays, ncrossings,2);
      if(temp1 || temp2)
      {

        //printf("Crossing %d at (%d,%d), (%e, %e), (%e,%e)\n",i ,temp1, temp2, crossesx[1,0,i], crossesz[1,0,i], dx*temp1+xmin, dz*temp2+zmin);
      }
    }
  }
}
