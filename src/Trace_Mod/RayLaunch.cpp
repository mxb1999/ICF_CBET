#include "Trace_interface.hpp"

using namespace std;
//uses x0 and z0 arrays to initialize the ray positions, then launches via Launch_Ray_XZ()
void trackRays()
{
  cout<<"Check\n";

  if(printUpdates)
  {
    cout << "Tracking Rays" << endl;
  }
  cout<<"Check 2\n";
  int beam;
  //tracking arrays
  double x0[nrays*nbeams];//initial x position of ray
  double z0[nrays*nbeams];//initial z position of ray
  double kx0[nrays*nbeams];//initial x velocity of ray
  double kz0[nrays*nbeams];//initial z velocity of ray
  double phase_x[nrays];//phase of ray
  double pow_x[nrays];//power delivery of ray
  //initializing arrays for beam 1
  span(phase_x,beam_min_z, beam_max_z, nrays);
  span(z0, beam_min_z, beam_max_z, nrays);
  rayinit* raycoor = new rayinit[nrays*nbeams];
  #pragma omp parallel for num_threads(threads)
  for(int i = 0; i < nrays; i++)
  {
    raycoor[i].xinit = xmin-(dt/courant_mult*c*0.5);;
    raycoor[i].zinit = z0[i] + offset-(dz/2)-(dt/courant_mult*c*0.5);//initial z position of ray
    raycoor[i].kxinit = 1.0;
    raycoor[i].kzinit = 0;
    raycoor[i].beam = 0;
    double absX = raycoor[i].xinit-xmin;
    double absZ = raycoor[i].zinit-zmin;
    int currX = (absX)/dx;
    currX += ((raycoor[i].xinit-currX*dx) > 0.5);
    int currZ = (absX)/dx;
    currZ += ((raycoor[i].xinit-currX*dx) > 0.5);
    raycoor[i].wpeinit = 0;
    kz0[i] = 0;
    pow_x[i] = exp(-1*pow(pow(phase_x[i]/sigma,2.0),4.0/2.0));
    phase_x[i] += offset;
    //Beam one lies along the z axis, x axis is constant
  }
  int finalts[nrays][nbeams];
  cout<<"Check 3 \n" << printUpdates;

  if(printUpdates)
  {
    cout << "Launching Beam 1" << endl;
  }
  cout <<  scientific;
  beam = 0;
  //Loop to launch the rays for beam 1, parallelized using OpenMP
  //#pragma omp parallel for num_threads(threads)
  double interpTerm[nrays*nbeams];
  span(x0 + nrays, beam_min_z, beam_max_z, nrays);
  #pragma omp parallel for num_threads(threads)
  for(int i = 0; i < nrays; i++)
  {
    interpTerm[i] = interp(phase_x, pow_x, z0[i], nrays);
    raycoor[i].kxinit = 0;
    raycoor[i].kzinit = 1.0;
    raycoor[i].beam = 1;
    raycoor[i].xinit = x0[i] - (dx/2)+(dt/courant_mult*c*0.5);
    raycoor[i].zinit = zmin-(dt/courant_mult*c*0.5);
    interpTerm[i+nrays] = interp(phase_x, pow_x, z0[i], nrays);
    raycoor[i+nrays].xinit = x0[i] - (dx/2)+(dt/courant_mult*c*0.5);
    raycoor[i+nrays].zinit = zmin-(dt/courant_mult*c*0.5);
    interpTerm[i] = interp(phase_x, pow_x, z0[i], nrays);
    interpTerm[i+nrays] = interp(phase_x, pow_x, z0[i], nrays);
    //phase_x[i] +=offset;
  }
 // for(int i = nrays; i < nrays*2; i++)
  //{
  //  kx0[i] = 0.0;
 //   kz0[i] = 1.0;
  //  x0[i] -= x0[i] (dx/2)+(dt/courant_mult*c*0.5);
 //   z0[i] = zmin-(dt/courant_mult*c*0.5);
 //   //phase_x[i] -= offset;
 // }
  rayLaunch<<<1024,512>>>(rayinit* init, double* edepcu, double* force, double* crossx_cu, double* crossz_cu, int nrays);
  for(int i = 0; i < nrays*nbeams;i++)
  {
    beam = i/nrays;
    int rnum = i%nrays;
    printf("%d %d %d\n", beam, rnum, nrays);
    launch_ray_XZ(x0[i],z0[i],kx0[i],kz0[i],uray_mult*interpTerm[i],i%nrays,beam);
  }
  //Reset the intial conditions for beam 2
  
  

  beam = 1;

  if(printUpdates)
  {
    cout << "Launching Beam 2" << endl;
  }
  cout<<"Check 4\n";
/*
  //Loop to launch beam 2 rays
  //#pragma omp parallel for num_threads(threads)
  for(int i = 0; i < nrays;i++)
  {
    z0[i] = zmin-(dt/courant_mult*c*0.5);
    x0[i] += (dx/2)-(dt/courant_mult*c*0.5);
    double interpNum = interp(phase_x, pow_x, x0[i], nrays);
    #pragma omp atomic update
    injected += uray_mult*interpNum;
    launch_ray_XZ(x0[i],z0[i],kx0[i],kz0[i], uray_mult*interpNum, i,beam);
  }
*/
  if(printUpdates)
  {
    cout << "Finished Launching Rays" << endl;
  }
}



//updating the intersections array to account for array intersections found via marked
void updateIntersections()
{

  if(printUpdates)
  {
    cout << "Updating Rays Intersections" << endl;
  }
  for(int i = 0; i < nx;i++)
  {
    for(int j = 0; j < nz; j++)
    {
      int s1 = 0;
      int s2 = 0;
      for(int q = 0; q < ncrossings;q++)
      {
          s1 += (vec4D(marked, 0, i, j, q, nx, nz, ncrossings) != 0);
          s2 += (vec4D(marked, 1, i, j, q, nx, nz, ncrossings) != 0);
      }
      vec2DW(intersections,i,j,nz,s1*s2);
    }
  }
}
void fillMarked()
{

  for(int i = 0; i < nbeams; i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      for(int q = 0; q < ncrossings;q++)
      {
        int tempx = vec4D(boxes, i,j,q,0, nrays, ncrossings, 2);//[i][j][q][0]
        int tempz = vec4D(boxes, i,j,q,1, nrays, ncrossings, 2);//[i][j][q][1]
        if(tempx && tempz)
        {
          tempx--;
          tempz--;
          vec4DW(marked, i, tempx,tempz,q, nx, nz, ncrossings, j+1);
          vec3DI(present, i, tempx, tempz, nx,nz,1);
          if(i == 1 && j == 0)
          {
            vec2DW(raypath, tempx, tempz, nz,1);
          }
          
          //marked[(tempx*nz + tempz)*nbeams + i].push(j+1);
        }else
        {
          break;
        }
        //if(tempx ||)
      }
    }
  }
}

void launchRays()
{
  trackRays();
  //fillMarked();

 // updateIntersections();
}
