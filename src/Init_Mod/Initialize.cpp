#include "init_interface.hpp"



int* getVarI(string target)
{
  
  std::map<std::string, int*> stringmap;
  stringmap["printUpdates"] = &printUpdates;
  stringmap["printTimings"] = &printTimings;
  stringmap["printCBETDiagnostics"] = &printCBETDiagnostics;
  stringmap["printRayTrackDiagnostics"] = &printRayTrackDiagnostics;
  stringmap["printHDF5Diagnostics"] = &printHDF5Diagnostics;
  stringmap["printSpecificTimings"] = &printSpecificTimings;
  stringmap["pyPlot"] = &pyPlot;
  stringmap["iterate"] = &iterate;
  stringmap["calcCBET"] = &calcCBET;
  stringmap["nx"] = &nx;
  stringmap["nz"] = &nz;
  stringmap["maxIter"] = &maxIter;
  stringmap["threads"] = &threads;
  stringmap["nbeams"] = &nbeams;
  stringmap["rays_per_zone"] = &rays_per_zone;
  stringmap["switchvar"] = &switchvar;

  int s = stringmap.size();
  int* result = stringmap[target];
  if(s != stringmap.size())
  {
    stringmap.erase(target);
    return NULL;
  }
  return result;
}
double* getVarD(string target)
{
  std::map<std::string, double*> stringmap;
  stringmap["lambda"] = &lambda;
  stringmap["estat"] = &estat;
  stringmap["mach"] = &mach;
  stringmap["Z"] = &Z;
  stringmap["Te"] = &Te;
  stringmap["Te_eV"] = &Te_eV;
  stringmap["Ti"] = &Ti;
  stringmap["Ti_eV"] = &Ti_eV;
  stringmap["iaw"] = &iaw;
  stringmap["xmin"] = &xmin;
  stringmap["xmax"] = &xmax;
  stringmap["zmin"] = &zmin;
  stringmap["zmax"] = &zmax;
  stringmap["maxIncr"] = &maxIncr;
  stringmap["converge"] = &converge;
  stringmap["courant_mult"] = &courant_mult;
  stringmap["intensity"] = &intensity;
  stringmap["offset"] = &offset;
  stringmap["beam_max_z"] = &beam_max_z;
  stringmap["beam_min_z"] = &beam_min_z;

  int s = stringmap.size();
  double* result = stringmap[target];
  if(s != stringmap.size())
  {
    stringmap.erase(target);
    return NULL;
  }
  return result;

}

void setInitParams()//import configuration file
{
  std::ifstream file;
  file.open("config/defaultconfig.txt", std::ifstream::in);
  char temp[256];
  if(file.is_open())
  {
    while(!file.eof())
    {
      file.getline(temp, sizeof(temp));
      std::string line(temp);
      int segment = line.find_first_of('=');
      std::string varname  = line.substr(0,segment);

      std::string val  = line.substr(segment+1);

      int* intName = getVarI(varname);
      double* dName = getVarD(varname);
      if(intName)
      {
        int flag = 1;
        if(val[0] == '-')
        {
          val = val.substr(1);
          flag = -1;
        }
        *intName = flag * std::stoi(val, NULL);
      }else if(dName)
      {
        int flag = 1;
        *dName = flag * std::stod(val, NULL);
      }else
      {
        std::cout << "Variable " << varname << " not found\n";
      }
    }

  }else
  {
    printf("Unable to open file\n");
  }
  fflush(stdout);
    file.close();
  //fill in dependent variables
  mi_kg = 10230.0*me;	   // Mass of ion in kg
  mi = 10230*(1.0e3*me);          // Mass of ion in g
  uray_mult = intensity*(courant_mult)*pow(double(rays_per_zone),-1.0); //multiplier which determines intensity deposited in a given zone
  dz = (zmax-zmin)/(nz-1);//dimensions of a single cell
  dx = (xmax-xmin)/(nx-1);
  nrays= int(rays_per_zone*(beam_max_z-beam_min_z)/dz)+0;//number of rays per beam
  dt=courant_mult*fmin(dx,dz)/c;//time stepping
  nt=int(pow(courant_mult,-1.0)*fmax(nx,nz)*2.0)+1;//number of time steps to track for a given ray
  numstored = nx*6;//number of rays stored per grid zone
  ncrossings = nx * 3;//max number of ray crossings that can be stored per ray
  freq = c/lambda;		// frequency of light, in Hz
  omega = 2*pi*freq;	// frequency of light, in rad/s
  ncrit = 1e-6*(pow(omega,2.0)*me*e0/pow(ec,2.0));
}
//dynamically allocate and initialize the arrays
void initialize()
{
  setInitParams();
    printf("(%d,%d)\n", nx, nz);

  cudaMallocManaged(&intersections, sizeof(int)*GRID);
  //intersections = new int[GRID]{0}; //nx nz
  cudaMallocManaged(&marked, sizeof(int)*GRID*RAYS);
  //marked = new int[GRID*nbeams*numstored]{0}; //nx nz nrays nbeams
  cudaMallocManaged(&dedendx, sizeof(int)*GRID);
  //dedendx = new double[GRID]; //nx nz
  cudaMallocManaged(&dedendz, sizeof(int)*GRID);
  //dedendz = new double[GRID]; //nx nz
  cudaMallocManaged(&x, sizeof(int)*nx);
  //x = new double[nx]{0.0}; //nx nz
  cudaMallocManaged(&marked, sizeof(int)*GRID*numstored*nbeams);
  
  cudaMallocManaged(&present, sizeof(int)*GRID*nbeams);

  cudaMallocManaged(&z, sizeof(int)*nz);
  //z = new double[nz]{0.0}; //nx nz
  cudaMallocManaged(&eden, sizeof(int)*GRID);
  //eden = new double[GRID]; //nx nz
  cudaMallocManaged(&edep, sizeof(int)*RAYS* (nx+2) * (nz+2));
  //edep = new double[nbeams * nx+2 * nz+2]{0.0}; //nx+2 nz+2 nrays
  cudaMallocManaged(&present, sizeof(int)*GRID);
  //present = new int[nbeams*GRID]; //nx nz nbeams
  cudaMallocManaged(&machnum, sizeof(int)*GRID);
  //machnum = new double[GRID]; //nx nz
  cudaMallocManaged(&boxes, sizeof(int)*RAYS*ncrossings);
  //boxes = new int[nbeams*nrays*ncrossings*2]{0}; //nbeams nrays ncrossings 2
  cudaMallocManaged(&u_flow, sizeof(int)*GRID);
  //u_flow = new double[GRID]; //nx nz
  cudaMallocManaged(&dkx, sizeof(int)*RAYS*ncrossings);
  //dkx = new double[nbeams*nrays*ncrossings]; //nbeams nrays 2
  cudaMallocManaged(&dkz, sizeof(int)*RAYS*ncrossings);
  //dkz = new double[nbeams*nrays*ncrossings]; //nbeams nrays 2
  cudaMallocManaged(&dkmag, sizeof(int)*RAYS*ncrossings);
  //dkmag = new double[nbeams*nrays*ncrossings]; //nbeams nrays 2
  cudaMallocManaged(&W, sizeof(int)*GRID*nbeams);
  //W = new double[nbeams*GRID];//nx nz
  cudaMallocManaged(&W_new, sizeof(int)*GRID*nbeams);
  //W_new = new double[nbeams*GRID];//nx nz
  cudaMallocManaged(&wpe, sizeof(int)*GRID);
  //wpe = new double[GRID]; //nx nz
  cudaMallocManaged(&crossesz, sizeof(int)*RAYS*ncrossings);
  //crossesz = new double[nbeams*nrays*ncrossings]; //nbeams nrays ncrossings
  cudaMallocManaged(&crossesx, sizeof(int)*RAYS*ncrossings);
  //crossesx = new double[nbeams*nrays*ncrossings]; //nbeams nrays ncrossings
  cudaMallocManaged(&ints, sizeof(int)*RAYS*ncrossings*numstored);
  //ints = new int[nbeams*nrays*ncrossings*numstored]; //nbeams nrays ncrossings
  cudaMallocManaged(&raypath, sizeof(int)*GRID);
  //raypath = new int[GRID];
  cudaMallocManaged(&mult, sizeof(int)*RAYS*ncrossings);
  //mult = new double[nbeams*nrays*ncrossings]{1.0};
  for(int i = 0; i < RAYS*ncrossings;i++)
  {
    mult[i] = 1;
  }
  gain2arr = new double[GRID];
  gain1arr = new double[GRID];
  mag = new double[GRID];

  auto check1 = chrono::high_resolution_clock::now();
  //Calculating the initial energy density, wpe, and machnum values
  span(x, xmin, xmax, nx);
  span(z, zmin, zmax, nz);
  
  for(int i = 0; i < nx;i++)
  {
    for(int j = 0; j < nz;j++)
    {
      //eden[i][j] = temp[i];
      double temp = fmax(0.0,((0.3*ncrit-0.1*ncrit)/(xmax-xmin))*(x[i]-xmin)+(0.1*ncrit));
      vec2DW(eden,i,j,nz,temp);

      temp = sqrt(vec2D(eden,i,j,nz)*1e6*pow(ec,2.0)/(me*e0));
      vec2DW(wpe,i,j,nz, temp);

      temp = fmax(0.0,(((-0.4)-(-2.4))/(xmax-xmin))*(x[i]-xmin))+(-2.4);
      vec2DW(machnum,i,j,nz, temp);
    }
  }
  for(int i = 0; i < nx-1; i++)
  {
    for(int j = 0; j < nz-1; j++)
    {
      double tempx = (vec2D(eden,i+1,j,nz)-vec2D(eden,i,j,nz))/(x[i+1]-x[i]);
      double tempz = (vec2D(eden,i,j+1,nz)-vec2D(eden,i,j,nz))/(z[j+1]-z[j]);
      vec2DW(dedendx,i,j,nz, tempx);
      vec2DW(dedendz,i,j,nz, tempz);
    }
  }
  for(int i = 0; i < max(nx,nz);i++)
  {
    if(i < nx)
    {
      double temp = vec2D(dedendz,i,nz-2,nz);
      vec2DW(dedendz,i,nz-1,nz,temp);

    }
    if(i < nx-1)
    {
      double tempx = (vec2D(eden,i+1,nz-1,nz)-vec2D(eden,i,nz-1,nz))/(x[i+1]-x[i]);
      vec2DW(dedendx,i,nz-1,nz, tempx);
    }
    if(i < nz-1)
    {
      double tempz = (vec2D(eden,nx-1,i+1,nz)-vec2D(eden,nx-1,i,nz))/(z[i+1]-z[i]);
      vec2DW(dedendz,nx-1,i,nz, tempz);
    }
    if(i < nz)
    {

      double temp = vec2D(dedendx,nx-2,i,nz);
      vec2DW(dedendx,nx-1,i,nz,temp);
    }
  }
}
