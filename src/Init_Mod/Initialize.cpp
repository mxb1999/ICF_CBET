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
}
//dynamically allocate and initialize the arrays
void initialize()
{
  setInitParams();
  auto start = chrono::high_resolution_clock::now();
  mi_kg = 10230.0*me;	   // Mass of ion in kg
  mi = 10230*(1.0e3*me);          // Mass of ion in g
  uray_mult = intensity*(courant_mult)*pow(double(rays_per_zone),-1.0); //multiplier which determines intensity deposited in a given zone
  dz = (zmax-zmin)/(nz-1);
  dx = (xmax-xmin)/(nx-1);
  nrays= int(rays_per_zone*(beam_max_z-beam_min_z)/dz)+0;//number of rays per beam
  dt=courant_mult*fmin(dx,dz)/c;//time stepping
  nt=int(pow(courant_mult,-1.0)*fmax(nx,nz)*2.0)+1;//number of time steps to track for a given ray
  numstored = nx*6;//number of rays stored per grid zone
  ncrossings = nx * 3;//max number of ray crossings that can be stored per ray

  freq = c/lambda;		// frequency of light, in Hz
  omega = 2*pi*freq;	// frequency of light, in rad/s
  ncrit = 1e-6*(pow(omega,2.0)*me*e0/pow(ec,2.0));
  //dynamic allocation (using pointers to access arrays so the stack is not filled)
  intersections = new double*[nx]; //nx nz
  marked = new int[nx*nz*nbeams*numstored]{0}; //nx nz nrays nbeams
  dedendx = new double*[nx]; //nx nz
  dedendz = new double*[nx]; //nx nz
  x = new double[nx]{0.0}; //nx nz
  z = new double[nz]{0.0}; //nx nz
  eden = new double*[nx]; //nx nz
  edep = new double**[nbeams]; //nx+2 nz+2 nbeams
  present = new int**[nbeams]; //nx nz nbeams
  machnum = new double*[nx]; //nx nz
  boxes = new int[nbeams*nrays*ncrossings*2]{0}; //nbeams nrays ncrossings 2
  u_flow = new double*[nx]; //nx nz
  dkx = new double**[nbeams]; //nbeams nrays 2
  dkz = new double**[nbeams]; //nbeams nrays 2
  dkmag = new double**[nbeams]; //nbeams nrays 2
  W = new double**[nbeams];//nx nz
  W_init = new double**[nbeams];//nx nz
  W_new = new double**[nbeams];//nx nz
  wpe = new double*[nx]; //nx nz
  crossesz = new double**[nbeams]; //nbeams nrays ncrossings
  crossesx = new double**[nbeams]; //nbeams nrays ncrossings
  ints = new int**[nbeams]; //nbeams nrays ncrossings
  raypath = new int*[nx];
  mult = new double[nbeams*nrays*ncrossings]{1.0};
  for(int i = 0; i < nbeams*nrays*ncrossings;i++)
  {
    mult[i] = 1;
  }
  gain2arr = new double*[nx];
  gain1arr = new double*[nx];
  mag = new double*[nx];

  auto check1 = chrono::high_resolution_clock::now();
  for(int i = 0; i < nx*nz*nbeams; i++)
  {
    //marked[i] = queue<int>();
  }
  auto check1_2 = chrono::high_resolution_clock::now();
  for(int i = 0; i < nx; i++)
  {
    gain2arr[i] = new double[nz];
    gain1arr[i] = new double[nz];
    mag[i] = new double[nz];
    raypath[i] = new int[nz]{0};
    intersections[i] = new double[nz]{0.0};
    eden[i] = new double[nz]{0.0};
    machnum[i] = new double[nz]{0.0};
    dedendx[i] = new double[nz]{0.0};
    dedendz[i] = new double[nz]{0.0};
    u_flow[i] = new double[nz]{0.0};
    wpe[i] = new double[nz]{0.0};
  }

  auto check2 = chrono::high_resolution_clock::now();

  for(int i = 0; i < nbeams; i++)
  {
    edep[i] = new double*[nx+2];
    //boxes[i] = new int**[nrays];
    present[i] = new int*[nx];
    dkx[i] = new double*[nrays];
    dkz[i] = new double*[nrays];
    dkmag[i] = new double*[nrays];
    crossesz[i] = new double*[nrays];
    crossesx[i] = new double*[nrays];
    ints[i] = new int*[nrays];
    for(int j = 0;j < nx+2;j++)
    {
      if(j < nx)
      {
        present[i][j] = new int[nz]{0};
      }
      edep[i][j] = new double[nz+2]{0.0};
    }
    for(int j = 0; j < nrays; j++)
    {
      dkx[i][j] = new double[ncrossings-1]{0.0};
      dkz[i][j] = new double[ncrossings-1]{0.0};
      dkmag[i][j] = new double[ncrossings-1]{0.0};
      crossesz[i][j] = new double[ncrossings]{0.0};
      crossesx[i][j] = new double[ncrossings]{0.0};
    //  boxes[i][j] = new int*[ncrossings];
      ints[i][j] = new int[ncrossings]{0};
      for(int m = 0; m < ncrossings;m++)
      {
    //    boxes[i][j][m] = new int[2]{0};
      }
    }
  }
  auto check3 = chrono::high_resolution_clock::now();
  if(printUpdates)
  {
    cout << "Setting initial conditions for ray tracker..." <<endl;
    cout << "nrays per beam is"<< nrays <<endl;
  }
  //Calculating the initial energy density, wpe, and machnum values
  span(x, xmin, xmax, nx);
  span(z, zmin, zmax, nz);
  double* temp1 = new double[nx];
  double* temp2 = new double[nx];
  span(temp1, 0.1, 0.4, nx);
  span(temp2, 1, 1.8, nx);
  for(int i = 0; i < nx;i++)
  {
    for(int j = 0; j < nz;j++)
    {
      //eden[i][j] = temp[i];
      eden[i][j] = fmax(0.0,((0.3*ncrit-0.1*ncrit)/(xmax-xmin))*(x[i]-xmin)+(0.1*ncrit));
      wpe[i][j] = sqrt(eden[i][j]*1e6*pow(ec,2.0)/(me*e0));
      machnum[i][j] = fmax(0.0,(((-0.4)-(-2.4))/(xmax-xmin))*(x[i]-xmin))+(-2.4);
    }
  }
  for(int i = 0; i < nx-1; i++)
  {
    for(int j = 0; j < nz-1; j++)
    {
      dedendx[i][j] = (eden[i+1][j]-eden[i][j])/(x[i+1]-x[i]);
      dedendz[i][j] = (eden[i][j+1]-eden[i][j])/(z[j+1]-z[j]);
    }
  }
  auto check4 = chrono::high_resolution_clock::now();
  for(int i = 0; i < fmax(nx,nz);i++)
  {
    if(i < nx)
    {
      dedendz[i][nz-1] = dedendz[i][nz-2];

    }
    if(i < nz)
    {
      dedendx[nx-1][i] = dedendx[nx-2][i];
    }
  }
  auto check5 = chrono::high_resolution_clock::now();
  if(printSpecificTimings)
  {
    cout << "Initialize CPU Time 1: " << chrono::duration_cast<chrono::milliseconds>(check1-start).count() << " seconds" << endl;
    cout << "Initialize CPU Time 2: " << chrono::duration_cast<chrono::milliseconds>(check2-check1).count() << " seconds" << endl;
    cout << "Initialize CPU Time 3: " << chrono::duration_cast<chrono::milliseconds>(check3-check2).count() << " seconds" << endl;
    cout << "Initialize CPU Time 4: " << chrono::duration_cast<chrono::milliseconds>(check4-check3).count() << " seconds" << endl;
    cout << "Initialize CPU Time 5: " << chrono::duration_cast<chrono::milliseconds>(check5-check4).count() << " seconds" << endl;
  }


}