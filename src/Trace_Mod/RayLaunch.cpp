#include "Trace_interface.hpp"
#include "parallelConfig.hpp"
using namespace std;
//uses x0 and z0 arrays to initialize the ray positions, then launches via Launch_Ray_XZ()
void trackRays()
{
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
  rayinit* raycoor;
  cudaMallocManaged(&raycoor, sizeof(rayinit)*RAYS);
    double interpTerm[nrays*nbeams];

  //#pragma omp parallel for num_threads(threads)
  for(int i = 0; i < nrays; i++)
  {
    rayinit* curr = raycoor + i;
    curr->xinit = xmin-(dt/courant_mult*c*0.5);
    curr->zinit = z0[i] + offset-((dz/2)+(dt/courant_mult*c*0.5));//initial z position of ray
    curr->kxinit = 1.0;
    curr->kzinit = 0;
    curr->beam = 0;
    double absX = curr->xinit-xmin;
    double absZ = curr->zinit-zmin;
    int currX = (absX)/dx;
    currX += ((curr->xinit-currX*dx) > 0.5);
    int currZ = (absX)/dx;
    currZ += ((curr->xinit-currX*dx) > 0.5);
    curr->wpeinit = 0;
    kz0[i] = 0;
    pow_x[i] = exp(-1*pow(pow(phase_x[i]/sigma,2.0),4.0/2.0));
    phase_x[i] += offset;
    //Beam one lies along the z axis, x axis is constant
  }
  int finalts[nrays][nbeams];

  if(printUpdates)
  {
    cout << "Launching Beam 1" << endl;
  }
  cout <<  scientific;
  beam = 0;
  //Loop to launch the rays for beam 1, parallelized using OpenMP
  //#pragma omp parallel for num_threads(threads)
  span(x0 + nrays, beam_min_z, beam_max_z, nrays);
  //#pragma omp parallel for num_threads(threads)
  for(int i = nrays; i < nrays*2; i++)
  {
    rayinit* curr = raycoor + i;
    curr->kxinit = 0;
    curr->kzinit = 1.0;
    curr->beam = 1;
    curr->xinit = x0[i] - ((dx/2)+(dt/courant_mult*c*0.5));
    curr->zinit = zmin-(dt/courant_mult*c*0.5);
    curr->urayinit =  interp(phase_x, pow_x, curr->xinit, nrays)*uray_mult;
    raycoor[i-nrays].urayinit =  interp(phase_x, pow_x, raycoor[i-nrays].zinit, nrays)*uray_mult;
    interpTerm[i] = interp(phase_x, pow_x, z0[i], nrays);
    //phase_x[i] +=offset;
  }
  //rayLaunch<<<1024,512>>>(rayinit* init, double* edepcu, double* force, double* crossx_cu, double* crossz_cu, int nrays);

  //Reset the intial conditions for beam 2
  printf("Check 1\n");
  fflush(stdout);
  LaunchCUDARays(raycoor);
  /*
  switch(deviceConfiguration->getType())
  {
    case CUDA:
      cout << "Launching Beam 1 via CUDA" << endl;
      LaunchCUDARays(raycoor);
      return;
      break;
    case OPENCL:
      break;
  }
  */
 /*
  for(int i = 0; i < nrays*nbeams;i++)
  {
    int rnum = i%nrays;
    double xinit = curr->xinit;
    double zinit = curr->zinit;
    double kxinit = curr->kxinit;
    double kzinit = curr->kzinit;
    beam = curr->beam;
    double urayinit = curr->urayinit;

    launch_ray_XZ(raycoor[i],rnum);
  }
  beam = 1;

  if(printUpdates)
  {
    cout << "Launching Beam 2" << endl;
  }

  //Loop to launch beam 2 rays
  //#pragma omp parallel for num_threads(threads)
*/
  if(printUpdates)
  {
    cout << "Finished Launching Rays" << endl;
  }
}



//updating the intersections array to account for array intersections found via marked
void updateIntersections()
{
  //intersections = new int
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
          //printf("Ree: %d %d %p\n", tempx,tempz, marked);
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

  updateIntersections();
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {
      if(vec3D(edep_flat,0,i,j,nx+2,nz+2) + vec3D(edep_flat,1,i,j,nx+2,nz+2) > 0)
      {
       //printf("(%d %d)\n", i,j);

      }
    }
  }
}
