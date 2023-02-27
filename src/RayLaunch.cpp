#include "Trace_interface.hpp"
#include "parallelConfig.hpp"
using namespace std;

void fillMarked()//marked temp indexed by rays
{
  for(int beam = 0; beam < nbeams; beam++)
  {
    //#pragma omp parallel for num_threads(threads)
    for(int raynum = 0; raynum < nrays; raynum++)
    {
      for(int i = 0; i < ncrossings;i++)
      {
        int id = vec3D(boxes, beam,raynum,i, nrays, ncrossings);
        int cx, cz;
        if(!id)
        {
          break;
        }
        cx = (id - 1) % nx;
        cz = (id - cx - 1)/nx;
        int val = vec3D(present, beam, cx,cz, nx,nz);
        vec4DW(marked, beam,cx,cz, val, nx,nz, numstored, raynum+1);
        vec3DI(present, beam, cx,cz, nx,nz,1);
      }
    }
  }
}


double* calculate_multiplier(double* areas)
{
  double* resultarr = (double*)malloc(sizeof(double)*CROSS);
  if(resultarr == NULL)
  {
    return NULL;
  }
  for(int i = 0; i < nbeams*nrays; i++)
  {
    for(int j = 0; j < ncrossings; j++)
    {
      int offset = (j == ncrossings-1);
      double averagearea = areas[i*ncrossings + j + offset] - areas[i*ncrossings + j] + (!offset)*areas[i*ncrossings + j];
      int ix = vec4D(boxes, i/nrays, i % nrays, j, 0, nrays, ncrossings, 2);
      int iz = vec4D(boxes, i/nrays, i % nrays, j, 1, nrays, ncrossings, 2);
      if(!ix || !iz)
      {
        break;
      }
      ix--;
      iz--;
      double neOverNc = eden[ix*nz + iz]/ncrit;
      double nediffinv = 1/sqrt(1-neOverNc);
      double epseff = 1/(nediffinv*nediffinv);
      resultarr[i*ncrossings + j] = nediffinv*averagearea*1/sqrt(epseff);
    }
  }
  return resultarr;
}
//uses x0 and z0 arrays to initialize the ray positions, then launches via Launch_Ray_XZ()
void trackRays()
{
  fillTraceArrays();
  if(printUpdates)
  {
    //cout << "Tracking Rays" << endl;
  }
  int beam;
  //tracking arrays
  double* x0 = new double[nrays*nbeams];//initial x position of ray
  double* z0 = new double[nrays*nbeams];//initial z position of ray
  double* phase_x = new double[nrays];//phase of ray
  double* pow_x = new double[nrays];//power delivery of ray
  //initializing arrays for beam 1
  span(phase_x,beam_min_z, beam_max_z, nrays);
  span(z0, beam_min_z, beam_max_z, nrays);
  rayinit* raycoor;
  if(cudaCalc)
  {
    cudaMallocManaged(&raycoor, sizeof(rayinit)*RAYS);
  }else
  {
    raycoor = (rayinit*)malloc(sizeof(rayinit)*RAYS);
  }

  //#pragma omp parallel for num_threads(threads)
  for(int i = 0; i < nrays; i++)
  {
    rayinit* curr = raycoor + i;
    curr->xinit = xmin;
    curr->zinit = z0[i] + offset;//initial z position of ray
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
    pow_x[i] = exp(-1*pow(pow(phase_x[i]/sigma,2.0),4.0/2.0));
    phase_x[i] += offset;

    //Beam one lies along the z axis, x axis is constant
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
    curr->xinit = x0[i];// - ((dx/2)+(dt/courant_mult*c*0.5));
    curr->zinit = zmin;//-(dt/courant_mult*c*0.5);
    curr->urayinit =  interp(phase_x, pow_x, curr->xinit, nrays)*uray_mult;
    raycoor[i-nrays].urayinit =  interp(phase_x, pow_x, raycoor[i-nrays].zinit, nrays)*uray_mult;
    //phase_x[i] +=offset;
  }
  //rayLaunch<<<1024,512>>>(rayinit* init, double* edepcu, double* force, double* crossx_cu, double* crossz_cu, int nrays);

  //Reset the intial conditions for beam 2
  delete [] x0;
  delete [] z0;
  delete [] phase_x;
  delete [] pow_x;

  if(cudaCalc)
  {
    LaunchCUDARays(raycoor);
    return;
  }
  auto startL = std::chrono::high_resolution_clock::now();
  areas = new double[nrays*ncrossings*nbeams];

  //#pragma omp parallel for num_threads(threads)
  for(int i = 0; i < nrays*nbeams;i++)
  {
    int rnum = i%nrays;
    rayinit* curr = raycoor + i;
    double xinit = curr->xinit;
    double zinit = curr->zinit;
    double kxinit = curr->kxinit;
    double kzinit = curr->kzinit;
    beam = curr->beam;
    double urayinit = curr->urayinit;

    launch_ray_XZ(raycoor[i],rnum);
  }
  auto stopL = std::chrono::high_resolution_clock::now();
  //Loop to launch beam 2 rays
  //#pragma omp parallel for num_threads(threads)

  if(printUpdates)
  {
    //cout << "Finished Launching Rays" << endl;
  }
  auto startI = std::chrono::high_resolution_clock::now();
  fillMarked();
  auto stopI = std::chrono::high_resolution_clock::now();

  //*output << "CPUTrace " << threads << " " << nrays << " " << chrono::duration_cast<chrono::milliseconds>(stopL-startL).count() << " " << chrono::duration_cast<chrono::milliseconds>(stopI-startI).count() << std::endl;
  cudaFree(raycoor);
}





void launchRays()
{
  //trackRays();
  //fillMarked();
}
