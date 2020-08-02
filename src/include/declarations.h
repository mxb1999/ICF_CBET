#include <vector>
double maxDev;
int beam;

//Launch Ray Values


double cs;
int ray1num;
double maxInc;
int iter;
int count = 0;
double injected = 0;
int gridcount = 0;
int** counter;
//Pointers for necessary arrays

double** intersections; //nx nz
int** marked; //nx nz numstored nbeams
double** dedendx; //nx nz
double** dedendz; //nx nz
double* x; //nx
double* z; //nz
double** eden; //nx nz

double*** edep; //nx+2 nz+2 nbeams
int*** present; //nx nz nbeams
double** machnum; //nx nz
int**** boxes; //nbeams nrays nx*3 2
bool**** boxTrack;
double**** W_storage; //nbeams nx nz nrays
double** u_flow; //nx nz
double*** dkx; //nbeams nrays 2
double*** dkz; //nbeams nrays 2
double*** dkmag; //nbeams nrays 2

//Launch_Ray_XZ specific arrays (all have a length of nt)

//CBET specific arrays
double*** W;//nx nz
double*** W_init;//nx nz
double*** W_new;//nx nz
double*** i_b;//nx nz
double*** i_b_prev; //nbeams nx nz
double*** i_b_new;//nx nz
double** wpe; //nx nz
double*** crossesz; //nbeams nrays ncrossings
double*** crossesx; //nbeams nrays ncrossings
int*** ints; //nbeams nrays ncrossings
int iteration = 0;
//arrays used only for plotting
double** i_bplot;//nx nz
double* convergeplot; //maxWrites
double** orderplot1; //nx nz
double** orderplot2; //nx nz
double** i_b1Error;
double** i_b2Error;
double** i_b_newplot;//nx nz
double** edenplot; //the array is eden/ncrit,  nx nz
double** edepplot; //nx nz
