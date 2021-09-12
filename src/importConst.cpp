#include "io_interface.hpp"
#include "IOVars.hpp"



int* getVarI(std::string target)
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
  stringmap["maxIterations"] = &maxIter;
  stringmap["threads"] = &threads;
  stringmap["nbeams"] = &nbeams;
  stringmap["rays_per_zone"] = &rays_per_zone;
  stringmap["switchvar"] = &switchvar;
  stringmap["cudaCalc"] = &cudaCalc;
  int s = stringmap.size();
  int* result = stringmap[target];
  if(s != stringmap.size())
  {
    stringmap.erase(target);
    return NULL;
  }
  return result;
}
double* getVarD(std::string target)
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
        //std::cout << "Variable " << varname << " not found\n";
      }
    }

  }else
  {
    printf("Unable to open file\n");
  }
  fflush(stdout);
    file.close();
  //fill in dependent variables
  double evToKelvin = 11604.5052;
  Ti = evToKelvin*Ti_eV;
  Te = evToKelvin*Te_eV;
  mi_kg = 10319.0*me;	   // Mass of ion in kg
  mi = 10319.0*(1.0e3*me);          // Mass of ion in g
  dz = (zmax-zmin)/(nz-1);//dimensions of a single cell
  dx = (xmax-xmin)/(nx-1);
  uray_mult = intensity*(courant_mult)*pow(double(rays_per_zone),-1.0); //multiplier which determines intensity deposited in a given zone
  nrays= int(rays_per_zone*(beam_max_z-beam_min_z)/dz)+0;//number of rays per beam
  dt=courant_mult*fmin(dx,dz)/c;//time stepping
  nt=int(pow(courant_mult,-1.0)*fmax(nx,nz)*2.0)+1;//number of time steps to track for a given ray
  numstored = nx*6;//number of rays stored per grid zone
  ncrossings = 222;//nx * 3;//max number of ray crossings that can be stored per ray
  freq = c/lambda;		// frequency of light, in Hz
  omega = 2*pi*freq;	// frequency of light, in rad/s
  ncrit = 1e-6*(pow(omega,2.0)*me*e0/pow(ec,2.0));
  double A = mi_kg/me/1836.0;
  double mp = 938273;
  printf("rays %d\n", nrays);
  cs = c*sqrt((Z*Te_eV/1e3+3.0*Ti_eV/1e3)/(A * mp));
}
//dynamically allocate and initialize the arrays
void initialize()
{
  setInitParams();

}
