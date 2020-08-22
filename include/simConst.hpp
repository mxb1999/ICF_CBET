#include <vector>
#include <string>
#include <iostream>
#include <cmath>


//Used to determine print statements, 1 = on, 0 = off
const int printUpdates = 1;
const int printTimings = 1;
const int printCBETDiagnostics = 0;
const int printRayTrackDiagnostics = 0;
const int printHDF5Diagnostics = 0;
const int printSpecificTimings = 1;
const int pyPlot = 1; //Determines whether the python script will automatically
const int iterate = 1;//determine whether the simulation will iterate




//Constants
const double sigma = 1.7e-4;
const double e0 =8.85418782e-12;
const double me =9.10938356e-31;
const double pi =3.14159265359;
const double ec = 1.60217662e-19;
const double lambda = 3.51e-5;	// wavelength of light, in cm. This is frequency-tripled "3w" or "blue" (UV) light
const double estat=4.80320427e-10; 	       // electron charge in statC
const double mach = -1.0*sqrt(2);                 // Mach number for max resonance
const double Z = 3.1;                        // ionization state
const double mi = 10230*(1.0e3*me);          // Mass of ion in g
const double mi_kg = 10230.0*me;	   // Mass of ion in kg
const double Te = 2.0e3*11604.5052;          // Temperature of electron in K
const double Te_eV = 2.0e3;
const double Ti = 1.0e3*11604.5052;          // Temperature of ion in K
const double c = 29979245800.0;              // Speed of light in cm/s

const double Ti_eV = 1.0e3;
const double iaw = 0.2;                      // ion-acoustic wave energy-damping rate (nu_ia/omega_s)!!
const double kb = 1.3806485279e-16;   //Boltzmann constant in erg/K
const double kb2 = 1.3806485279e-23;   //Boltzmann constant in J/K
const double freq = c/lambda;		// frequency of light, in Hz
const double omega = 2*pi*freq;	// frequency of light, in rad/s
const double ncrit = 1e-6*(pow(omega,2.0)*me*e0/pow(ec,2.0));

//Grid Constants for X by Z grid
const int nx=100; const float xmin = -5.0e-4; const float xmax=5.0e-4; const float dx = (xmax-xmin)/(nx-1);
const int nz=100; const float zmin = -5.0e-4; const float zmax=5.0e-4; const float dz = (zmax-zmin)/(nz-1);

//Constants pertaining to iteration and parallelization
const double maxIncrement = 0.2;
const int maxIterations = 100*iterate;
//Number of parallel threads
const int threads = 12;
//Fractional convergence cutoff (when the field stops changing by (converge*100) %)
const double converge = 1e-6;

//Beam/Ray Tracking parameters
const int nbeams = 2; //number of interacting beams
const int rays_per_zone = 5; //Rays launched per grid zone
const double intensity = 1.0e17; // W/cm^2
const float courant_mult = 0.2; // 0.37 // 0.25 // 0.36 // 0.22; Multiplier used to determine time stepping
const double uray_mult = intensity*(courant_mult)*pow(double(rays_per_zone),-1.0); //multiplier which determines intensity deposited in a given zone
const double offset = 0;//Determines offset of beam along axis
const double beam_max_z = 3.0e-4; const double beam_min_z = -3.0e-4;//determines the width of the beam
const int nrays= int(rays_per_zone*(beam_max_z-beam_min_z)/dz)+0;//number of rays per beam
const double dt=courant_mult*fmin(dx,dz)/c;//time stepping
const int nt=int(pow(courant_mult,-1.0)*fmax(nx,nz)*2.0)+1;//number of time steps to track for a given ray
const int numstored = nx*6;//number of rays stored per grid zone
const int ncrossings = nx * 3;//max number of ray crossings that can be stored per ray
