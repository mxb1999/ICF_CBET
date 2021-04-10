#include "Trace_interface.hpp"
#include "parallelConfig.hpp"
using namespace std;

void fillTempMarked(int* markedTemp)//marked temp indexed by rays
{
  for(int beam = 0; beam < nbeams; beam++)
  {
    #pragma omp parallel for num_threads(threads)
    for(int raynum = 0; raynum < nrays; raynum++)
    {
      for(int i = 0; i < ncrossings;i++)
      {
        int cx = vec4D(boxes, beam,raynum,i,0, nrays, ncrossings, 2);
        int cz = vec4D(boxes, beam,raynum,i,1, nrays, ncrossings, 2);
        if(!cx || !cz)
        {
          break;
        }
        cx--;
        cz--;
        vec4DW(markedTemp, beam, cx,cz,raynum, nx, nz, nrays,1);
      }
    }
  }
  
}
void fillMarked(int* markedTemp)//marked temp indexed by rays
{
  #pragma omp parallel for num_threads(threads)
  for(int ix = 0; ix < nx; ix++)
  {
    for(int iz = 0; iz < nz; iz++)
    {
      for(int i = 0; i < nbeams;i++)
      {
        vec3DW(present, i, ix,iz, nx,nz,0);
        int cnt = 0;
        for(int j = 0; j < nrays;j++)
        {
          if(vec4D(markedTemp, i,ix,iz, j, nx, nz, nrays))
          {
            vec4DW(marked, i,ix,iz, cnt, nx,nz, numstored, j+1);
            vec3DI(present, i, ix,iz, nx,nz,1);
            cnt++;
          }
        }
    
      }
    }
  }
  //printf("(%d, %d): %p %p\n", ix, iz,present, present+nx*nz*nbeams);
  
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
  if(cudaCalc)
  {
    LaunchCUDARays(raycoor);
    return;
  }
  auto startL = std::chrono::high_resolution_clock::now();
  

  #pragma omp parallel for num_threads(threads)
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
  int* markedTemp = new int[GRID*RAYS];
  fillTempMarked(markedTemp);
  fillMarked(markedTemp);
  auto stopI = std::chrono::high_resolution_clock::now();

  //std::cout << threads << " " << nrays << " " << chrono::duration_cast<chrono::milliseconds>(stopL-startL).count() << " " << chrono::duration_cast<chrono::milliseconds>(stopI-startI).count() << std::endl;
  delete [] markedTemp;
}





void launchRays()
{
  trackRays();
}
