  //CBET specific arrays
  double* machnum; //nx nz
  double* u_flow; //nx nz
  double* dk; //nbeams nrays 2
  double* dkmag; //nbeams nrays 2
  double* W;//nx nz
  double* W_init;//nx nz
  double* W_new;//nx nz
  double* i_b;//nx nz
  double* i_b_prev; //nbeams nx nz
  double* i_b_new;//nx nz
  double* wMult;
  double* spatialLog;
  double* neovernc;
  int* interactions_ML;