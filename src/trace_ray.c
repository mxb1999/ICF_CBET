#include "trace_rays.h"


TraceResult* new_TraceResult(path_t* trajectories, laser_t* vectors, laser_t* areas)
{

};
void trace_ray(Grid* grid, TraceInput input, int nt)
{
    int thisx = 0;
    int thisy = 0;
    int thisx_0 = 0;
    int thisy_0 = 0;
    spatial_t x_init = input.init_location[0], y_init = input.init_location[1];
    spatial_t kx_init = input.init_vector[0], ky_init = input.init_vector[1];
    plasma_t ncrit = NCRIT;
    spatial_t xmin = XMIN, xmax = XMAX, ymin = YMIN, ymax = YMAX;
    int ny = NY, nx = NX;
    spatial_t dx = DX, dy = DY;
    laser_t omega = OMEGA;
    int raynum = input.ray;
    path_t* trajectory = input.trajectories + raynum*CROSSINGS;
    spatial_t* area = input.areas + raynum*CROSSINGS;
    spatial_t* wave_vector = input.wave_vectors + raynum*CROSSINGS;
    //spatial_t uray;
    //int //markingxprev;// = new spatial_t[nt]{0.0};
    //int //markingyprev;// = new spatial_t[nt]{0.0};
    //int //markingx;// = new spatial_t[nt]{0.0};
    //int //markingy;
    spatial_t my_x = x_init;
    spatial_t my_y = y_init;
    spatial_t myvx;
    spatial_t myvy;
    spatial_t mykx;
    spatial_t myky;
    spatial_t myvxprev;
    spatial_t myvyprev;
    spatial_t lastcrossx = x_init;
    spatial_t lastcrossy = y_init;
    //determining the initial x grid index within the desired range to track the beam
    double dt = COURANT*dx/c;

    //determining the initial y grid index within the desired range to track the beam
    for(int i = 0;i < ny;i++)
    {
      if(y_init - (dy*i - ymin) <= ((0.5+1.0e-10)*dy + 1e-11) && y_init - (dy*i - ymin) >= -1*((0.5+1.0e-10)*dy + 1e-11) )
      {

        thisy_0=i;
        break;
      }
    }

    for(int i = 0;i < nx;i++)
    {

      spatial_t initX = x_init;
      spatial_t currX = (dx*i - xmin);
      if(initX - currX <= ((0.5+1.0e-10)*dx) && initX - currX >= -1*((0.5+1.0e-10)*dx) )
      {
        thisx_0=i;
        break;
      }
    }
    //determining the velocity characteristics of the ray based upon its initial position
    spatial_t k = sqrt((pow(omega,2.0)-pow(WPE_2D(thisx_0, thisy_0),2.0))/pow(c,2.0));
    spatial_t knorm = sqrt(pow(kx_init,2.0)+pow(ky_init,2.0));
    mykx=(kx_init/knorm)*k;			// Normaliyed value for the ray's initial k_x
    myky=(ky_init/knorm)*k;			// Normaliyed value for the ray's initial k_y
    myvxprev = pow(c,2.0)*mykx/omega;                   // v_group, group velocity (dw/dk) from D(k,w).
    myvyprev =  pow(c,2.0)*myky/omega;
    ////markingxprev = thisx_0;
    ////markingyprev = thisy_0;
    //__________Time Stepping__________
      int numcrossing = 0;

      //looping through time intervals
      for(int i = 1; i < nt;i++)
      {
        int offsetx = (thisx_0 == nx-1);
        int offsety = (thisy_0 == ny-1);
        plasma_t dedendx = (EDEN_2D(thisx_0+1-offsetx, thisy_0) - EDEN_2D(thisx_0-offsetx, thisy_0))/dx;
        plasma_t dedendy = (EDEN_2D(thisx_0, thisy_0+1-offsety) - EDEN_2D(thisx_0, thisy_0-offsety))/dy;
        spatial_t forcex = pow(c,2.0)/(2.0*ncrit)*dedendx;
        spatial_t forcey = pow(c,2.0)/(2.0*ncrit)*dedendy;

        myvx = myvxprev - forcex*dt;
        myvy = myvyprev - forcey*dt;
        my_x += myvx*dt;
        my_y += myvy*dt;

        int search_index_x = 1;
        int search_index_y = 1;
        int thisx_m = fmax(0, thisx_0-search_index_x );
        int thisx_p = fmin(nx-1, thisx_0+search_index_x);
        int thisy_m = fmax(0, thisy_0-search_index_y);
        int thisy_p = fmin(ny-1, thisy_0+search_index_y);
      /* //determining the current x index of the ray
        for(int j = thisx_m; j <= thisx_p;j++)
        {
          if ( my_x - (dx*j - xmin) <= ((0.5+1.0e-10)*dx + 1e-12) && my_x - (dx*j - xmin) >= -1*((0.5+1.0e-10)*dx + 1e-12))
          {
            thisx = j;
            break;
          }
        }
        if(raynum == 0 && beam == 0)
        {

        // printf("check 2 %d\n",i);
        }
        //determining the current y index of the ray
        for(int j = thisy_m; j <= thisy_p; j++)
        {
          if (my_y - (dy*j - ymin) <= ((0.5+1.0e-10)*dy + 1e-12) && my_y - (dy*j - ymin) >= -1*((0.5+1.0e-10)*dy + 1e-12))
          {

            thisy = j;
            break;
          }
        }*/
        int discreteX = ((my_x-xmin)/dx);
        int discretey = ((my_y-ymin)/dy);

        thisx = discreteX + (my_x-(discreteX+xmin)*dx > dx/2);
        thisy = discretey + (my_y-(discretey+ymin)*dy > dy/2);
        //spatial_t liney[2]={y_init, my_y};
        //spatial_t linex[2]={x_init, my_x};
        spatial_t lastx = 10000;
        spatial_t lasty = 10000;
        //iterating through the selected portions of the x spatial tracking arrays
        //boxes stores the spatial locations of each crossing of each ray
        //Marked = trajectory of a single ray, boxes = coordinates of each ray intersection

        for(int j = thisx_m; j <= thisx_p;j++)
        {
          spatial_t currx = (dx*j - xmin)-dx/2;//crossing into
          //if the ray is currently between within the desired caustic yone for a crossing
          if((my_x > currx && x_init <= (currx)) || (my_x < currx && x_init >= (currx)))
          {
            spatial_t m = (my_y - y_init)/(my_x-x_init);
            spatial_t b = my_y - my_x*m;
            spatial_t crossx = m*currx+b;
            //if the ray has moved since last update
            if(abs(crossx-lasty)>1.0e-20)
            {
              vec3DW(crossesx, beam, raynum, numcrossing, nrays, ncrossings,currx);
              vec3DW(crossesy, beam, raynum, numcrossing, nrays, ncrossings,crossx);
              //if ray is still within the grid
              //if(my_x < (xmax+dx/2) && my_x > (xmin-dx/2))
             // {
                //printf("(%f, %f), ([%f,%f], [%f, %f]) \n", my_x, my_y, dx*thisx+ xmin, dx*thisx+dx + xmin, dy*thisy+ ymin, dy*thisy+dy+ ymin);
              ACCESS3D(trajectory)
                //vec4DW(marked, beam,raynum,thisx,thisy, nrays, nx, ny, 1);
          //    }

              lastx = currx;
              numcrossing += 1;
              break;
            }
          }
        }

        //iterating through the selected portions of the y spatial tracking arrays
        //Same idea as previous loop, but for the y coordinate instead of x
          for(int j = thisy_m; j <= thisy_p;j++)//for [thisy_m, thisy_p] previous y locations, iterate through spatial locations centered on thisy
          {
            spatial_t curry = (dy*j - ymin)-dy/2;//center of the jth yone
            //printf("REE\n");
            if((my_y > (curry) && y_init < (curry + 1e-10)) || (my_y < (curry) && y_init > (curry - 1e-10)))//if my_y is approximately equal to a yone crossing
            {

              spatial_t m = (my_x - x_init)/(my_y-y_init);
              spatial_t b = my_x - my_y*m;
              spatial_t crossy = m*curry+b;
              if(abs(crossy-lastx) > 1.0e-20)
              {
                vec3DW(crossesx, beam, raynum, numcrossing, nrays, ncrossings,crossy);
                vec3DW(crossesy, beam, raynum, numcrossing, nrays, ncrossings,curry);
                if(my_y < (ymax+dy/2 +1e-10) && my_y > (ymin-dy/2-1e-10))
                {
                vec4DW(boxes, beam,raynum,numcrossing,0, nrays, ncrossings, 2, thisx+1);//[beam][raynum][numcrossing][0]
                vec4DW(boxes, beam,raynum,numcrossing,1, nrays, ncrossings, 2, thisy+1);//[beam][raynum][numcrossing][1]
                // vec4DW(marked, beam,raynum,thisx,thisy, nrays, nx, ny, 1);
                }
                lasty = curry;

                numcrossing += 1;
                break;
              }
            }
          }
          //Sets the "initial conditions" for the next iteration
          thisx_0 = thisx;
          thisy_0 = thisy;
          //markingx = thisx;
          //markingy = thisy;
          //Deposit energy due to the incident ray
          /*spatial_t increment = thisInit.urayinit;
          spatial_t xp = (my_x - (x[thisx]+dx/2.0))/dx;
          spatial_t yp = (my_y - (y[thisy]+dy/2.0))/dy;
          int xadd = (xp >= 0) ? 1 : -1;
          int yadd = (yp >= 0) ? 1 : -1;
          spatial_t dl = yp * yadd;
          spatial_t dm = xp * xadd;
          spatial_t a1 = (1.0-dl)*(1.0-dm);		// blue		: (x  , y  )
          spatial_t a2 = (1.0-dl)*dm;		// green	: (x+1, y  )
          spatial_t a3 = dl*(1.0-dm);		// yellow	: (x  , y+1)
          spatial_t a4 = dl*dm;			// red 		: (x+1, y+1)
          if((xadd != xadd) && (yadd != yadd))
          {

            break;
          }
          //printf("edep: %d %d %d %d %d\n",thisx+1, thisy+1, beam, raynum, RAYS);
          vec4DI(edep, beam, raynum, thisx+1, thisy+1, nrays, nx+2,ny+2, a1*increment);	// blue
          vec4DI(edep, beam, raynum, thisx+xadd+1, thisy+1, nrays, nx+2,ny+2, a2*increment);// green
          vec4DI(edep, beam, raynum, thisx+1, thisy+yadd+1, nrays, nx+2,ny+2, a3*increment);// yellow
          vec4DI(edep, beam, raynum, thisx+xadd+1, thisy+yadd+1, nrays, nx+2,ny+2, a4*increment);	// red

          */
          x_init = my_x;
          y_init = my_y;
          myvxprev = myvx;
          myvyprev = myvy;
          ////markingxprev = //markingx;
          ////markingyprev = //markingy;
        if ( (my_x < (xmin-(dx/2.0))) || (my_x > (xmax+(dx/2.0))))
        {

          break;                  // "breaks" out of the i loop once the if condition is satisfied
        } else if ( (my_y < (ymin-(dy/2.0))) || (my_y > (ymax+(dy/2.0)))){
            // the "|" means "or" (symbol above the return key)

          break;
        }
      }
}