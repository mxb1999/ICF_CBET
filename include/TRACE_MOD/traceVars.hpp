  //Declare variables here to avoid multiple instantiations
  int* intersections; //nx nz
  int* marked; //nx nz numstored nbeams
  double* dedendx; //nx nz
  double* dedendz; //nx nz
  double* x; //nx
  double* z; //nz
  double* eden; //nx nz
  double* wpe; //nx nz
  double* crossesz; //nbeams nrays ncrossings
  double* crossesx; //nbeams nrays ncrossings
  int* ints; //nbeams nrays ncrossings
  double* edep; //nx+2 nz+2 nbeams
  int* present; //nx nz nbeams
  int* boxes; //nbeams nrays nx*3 2
  double* edep_flat;
