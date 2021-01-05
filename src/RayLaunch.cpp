#include "implSim.hpp"
#include "customMath.hpp"
using namespace std;
//uses x0 and z0 arrays to initialize the ray positions, then launches via Launch_Ray_XZ()
void trackRays()
{
  if(printUpdates)
  {
    cout << "Tracking Rays" << endl;
  }

  //tracking arrays
  double x0[nrays];//initial x position of ray
  double z0[nrays];//initial z position of ray
  double kx0[nrays];//initial x velocity of ray
  double kz0[nrays];//initial z velocity of ray
  double phase_x[nrays];//phase of ray
  double pow_x[nrays];//power delivery of ray
  //initializing arrays for beam 1
  span(phase_x,beam_min_z, beam_max_z, nrays);
  span(z0, beam_min_z, beam_max_z, nrays);
  for(int i = 0; i < nrays; i++)
  {
    kx0[i] = 1.0;
    kz0[i] = 0;
    pow_x[i] = exp(-1*pow(pow(phase_x[i]/sigma,2.0),4.0/2.0));
    phase_x[i] += offset;
    //Beam one lies along the z axis, x axis is constant
    x0[i] = xmin-(dt/courant_mult*c*0.5);
    z0[i] += offset-(dz/2)-(dt/courant_mult*c*0.5);//initial z position of ray
  }
  int finalts[nrays][nbeams];

  if(printUpdates)
  {
    cout << "Launching Beam 1" << endl;
  }
  cout <<  scientific;
  beam = 0;
  //Loop to launch the rays for beam 1, parallelized using OpenMP
  #pragma omp parallel for num_threads(threads)
  for(int i = 0; i < nrays;i++)
  {
    double interpNum = interp(phase_x, pow_x, z0[i], nrays);
    #pragma omp atomic update
    injected += uray_mult*interpNum;
    launch_ray_XZ(x0[i],z0[i],kx0[i],kz0[i],uray_mult*interpNum,i);
  }
  //Reset the intial conditions for beam 2
  span(x0, beam_min_z, beam_max_z, nrays);
  for(int i = 0; i < nrays; i++)
  {
    kx0[i] = 0.0;
    kz0[i] = 1.0;
    x0[i] -= (dx/2)+(dt/courant_mult*c*0.5);
    z0[i] = zmin-(dt/courant_mult*c*0.5);
    phase_x[i] -= offset;
  }

  beam = 1;

  if(printUpdates)
  {
    cout << "Launching Beam 2" << endl;
  }

  //Loop to launch beam 2 rays
  #pragma omp parallel for num_threads(threads)
  for(int i = 0; i < nrays;i++)
  {
    z0[i] = zmin-(dt/courant_mult*c*0.5);
    x0[i] += (dx/2)-(dt/courant_mult*c*0.5);
    double interpNum = interp(phase_x, pow_x, x0[i], nrays);
    #pragma omp atomic update
    injected += uray_mult*interpNum;
    launch_ray_XZ(x0[i],z0[i],kx0[i],kz0[i], uray_mult*interpNum, i);
  }

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
  for(int i = 0; i < nx; i++)
  {
    for(int j = 0; j < nz; j++)
    {
      intersections[i][j] = 0;
    }
  }
  for(int i = 1; i < nx;i++)
  {
    for(int j = 1; j < nz; j++)
    {
      int temp1;
      int size1 = marked[(i*nz+j)*nbeams].size();
      //if at least two beams are in the same x and z coordinates, update intersections
      if((temp1 = marked[(i*nz+j)*nbeams].front()) == 0)
      {
        marked[(i*nz+j)*nbeams].pop();
        marked[(i*nz+j)*nbeams].push(temp1);
        break;
      }else
      {
        marked[(i*nz+j)*nbeams].pop();
        marked[(i*nz+j)*nbeams].push(temp1);
        int temp2;
        int size2 = marked[(i*nz+j)*nbeams+1].size();
        for(int l = 0; l < size2; l++)
        {
          //break if no ray from the other beam is found
          if((temp2 = marked[(i*nz+j)*nbeams+1].front()) == 0)
          {
            marked[(i*nz+j)*nbeams+1].pop();
            marked[(i*nz+j)*nbeams+1].push(temp2);
            break;
          }else
          {
            marked[(i*nz+j)*nbeams+1].pop();
            marked[(i*nz+j)*nbeams+1].push(temp2);
            intersections[i][j] += 1.0;
          }
        }
      }
    }
  }
}

void launchRays()
{
  trackRays();
  updateIntersections();
}
