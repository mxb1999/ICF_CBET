 
 GConfig* deviceConfiguration;

  //execution option variables
  int printUpdates,printTimings,printCBETDiagnostics,printRayTrackDiagnostics,printHDF5Diagnostics,printSpecificTimings, pyPlot,iterate,calcCBET, cudaCalc, optimize;
  //CBET constant variables
  double lambda,estat,mach,Z,mi,mi_kg,Te,Te_eV,Ti,Ti_eV,iaw,ncrit,freq,omega;
  int switchvar;
  //spatial information
  int nx, nz;
  double xmin, xmax, zmin, zmax, dx, dz;

  double maxIncr, converge;
  int threads,maxIter;
std::ofstream* output;
  int* raypath;
  int nbeams, nrays, nt, numstored, rays_per_zone, ncrossings;
  //arrays used only for plotting
  double* i_bplot;//nx nz
  double* convergeplot; //maxWrites
  double* orderplot1; //nx nz
  double* orderplot2; //nx nz
  double* i_b1Error;
  double* i_b2Error;
  double* mult;//nbeams nbeams-1 nrays ncrossings
  double* i_b_newplot;//nx nz
  double* edenplot; //the array is eden/ncrit,  nx nz
  double* edepplot; //nx nz
  double* ib_orig;//nx nz
  int* anyInt;//nx nz
  double* perturbation;//nx nz
  double intensity, offset, uray_mult, beam_min_z, beam_max_z, dt, courant_mult;
  //debugging arrays/values
  double* gain2arr;
  double* gain1arr;
  double* mag;
  short* rayAdjList;
  double cs;
  int ray1num;
  double maxInc;
  int iter;
  int count = 0;
  double injected = 0;
  int gridcount = 0;
  //Fundamental Constants
  double sigma = 1.7e-4;
  double e0 =8.85418782e-12;
  double me =9.10938291e-31;
  double pi =3.14159265359;
  double kb= 1.3806485279e-16;   //Boltzmann constant in erg/K
  double kb2= 1.3806485279e-23;   //Boltzmann constant in J/K
  double ec= 1.60217662e-19;
  double c= 29979245800.0;              // Speed of light in cm/s