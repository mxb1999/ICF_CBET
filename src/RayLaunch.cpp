#include "include/implSim.h"
#include "include/customMath.h"
using namespace std;
//uses x0 and z0 arrays to initialize the ray positions, then launches via Launch_Ray_XZ()
void trackRays()
{
  cout << "Tracking Rays" << endl;
  //tracking arrays
  double x0[nrays];
  double z0[nrays];
  double kx0[nrays];
  double kz0[nrays];
  double phase_x[nrays];
  double pow_x[nrays];
  span(phase_x,beam_min_z, beam_max_z, nrays);
  for(int i = 0; i < nrays; i++)
  {
    kx0[i] = 1.0;
    kz0[i] = -0.1;
    pow_x[i] = exp(-1*pow(pow(phase_x[i]/sigma,2.0),4.0/2.0));
    phase_x[i] += offset;
  }
  int finalts[nrays][nbeams];

  cout << "BEAMNUM is" << beam << endl;
  span(z0, beam_min_z, beam_max_z, nrays);
  cout <<  scientific;

  for(int i = 0; i < nrays; i++)
  {
    x0[i] = xmin-(dt/courant_mult*c*0.5);
    //cout << z0[i] << endl;
    z0[i] += offset-(dz/2)-(dt/courant_mult*c*0.5);

  }
  beam = 0;
  //#pragma omp parallel for num_threads(12)
  for(int i = 0; i < nrays;i++)
  {
    double interpNum = interp(phase_x, pow_x, z0[i], nrays);


    #pragma omp atomic update
    injected += uray_mult*interpNum;

    launch_ray_XZ(x0[i],z0[i],kx0[i],kz0[i],uray_mult*interpNum,i);
  }
  for(int i = 0; i < nrays; i++)
  {
    kx0[i] = 0.0;
    kz0[i] = 1.0;
    phase_x[i] -= offset;
  }

  beam = 1;

  cout << "BEAMNUM is" << beam << endl;

  span(x0, beam_min_z, beam_max_z, nrays);
  for(int i = 0; i < nrays;i++)
  {
    x0[i] -= (dx/2)+(dt/courant_mult*c*0.5);
    z0[i] = zmin-(dt/courant_mult*c*0.5);
  }
  //#pragma omp parallel for num_threads(12)
  for(int i = 0; i < nrays;i++)
  {
    z0[i] = zmin-(dt/courant_mult*c*0.5);
    x0[i] += (dx/2)-(dt/courant_mult*c*0.5);
    double interpNum = interp(phase_x, pow_x, x0[i], nrays);
    #pragma omp atomic update
    injected += uray_mult*interpNum;
    launch_ray_XZ(x0[i],z0[i],kx0[i],kz0[i], uray_mult*interpNum, i);
  }

  printf("%s\n", "Finished Launching Rays");
}



//updating the intersections array to account for array intersections found via marked
void updateIntersections()
{
  cout << "Updating Rays Intersections" << endl;
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
      for(int m = 0; m < nrays; m++)
      {
        //if at least two beams are in the same x and z coordinates, update intersections
        if(marked[i*nz+j][m*nbeams+0] == 0)
        {
          break;
        }else
        {
          for(int l = 0; l < nrays; l++)
          {
            if(marked[i*nz+j][l*nbeams+1] == 0)
            {
              break;
            }else
            {
              intersections[i][j] += 1.0;
            }
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
