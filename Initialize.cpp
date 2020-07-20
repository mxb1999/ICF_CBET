#include "implSim.h"
#include "declarations.h"

using namespace std;
//dynamically allocate and initialize the arrays
void initialize()
{

  auto start = chrono::high_resolution_clock::now();

  //dynamic allocation (using pointers to access arrays so the stack is not filled)
  intersections = new double*[nx]; //nx nz
  marked = new int*[nx * nz]; //nx nz nrays nbeams
  dedendx = new double*[nx]; //nx nz
  dedendz = new double*[nx]; //nx nz
  x = new double[nx]{0.0}; //nx nz
  z = new double[nz]{0.0}; //nx nz
  eden = new double*[nx]; //nx nz
  edep = new double**[nbeams]; //nx+2 nz+2 nbeams
  present = new int**[nbeams]; //nx nz nbeams
  machnum = new double*[nx]; //nx nz
  boxes = new int***[nbeams]; //nbeams nrays ncrossings 2
  W_storage = new double***[nbeams]; //nbeams nrays nx nz
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
  auto check1 = chrono::high_resolution_clock::now();
  for(int i = 0; i < nx*nz; i++)
  {
    marked[i] = new int[nrays*nbeams]{0};
  }

  auto check1_2 = chrono::high_resolution_clock::now();
    for(int i = 0; i < nx; i++)
    {
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
    W_storage[i] = new double**[nrays];
    W[i] = new double*[nx];
    W_new[i] = new double*[nx];
    W_init[i] = new double*[nx];
    edep[i] = new double*[nx+2];
    boxes[i] = new int**[nrays];
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
        W[i][j] = new double[nz]{0};
        W_new[i][j] = new double[nz]{0};
        W_init[i][j] = new double[nz]{0};
        present[i][j] = new int[nz]{0};
      }
      edep[i][j] = new double[nz+2]{0.0};
    }
    for(int j = 0; j < nrays; j++)
    {
      W_storage[i][j] = new double*[nx];
      dkx[i][j] = new double[ncrossings-1]{0.0};
      dkz[i][j] = new double[ncrossings-1]{0.0};
      dkmag[i][j] = new double[ncrossings-1]{0.0};
      crossesz[i][j] = new double[ncrossings]{0.0};
      crossesx[i][j] = new double[ncrossings]{0.0};
      boxes[i][j] = new int*[ncrossings];
      ints[i][j] = new int[ncrossings]{0};
      for(int m = 0; m < ncrossings;m++)
      {
        if(m < nx)
        {
          W_storage[i][j][m] = new double[nz]{0.0};
        }
        boxes[i][j][m] = new int[2]{0};
      }
    }
  }
  auto check3 = chrono::high_resolution_clock::now();

  cout << "Setting initial conditions for ray tracker..." <<endl;
  cout << "nrays per beam is"<< nrays <<endl;

  //Calculating the initial energy density, wpe, and machnum values
  span(x, xmin, xmax, nx);
  span(z, zmin, zmax, nz);
  for(int i = 0; i < nx;i++)
  {
    for(int j = 0; j < nz;j++)
    {
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
  cout << "Initialize CPU Time 1: " << chrono::duration_cast<chrono::milliseconds>(check1-start).count() << " seconds" << endl;
  cout << "Initialize CPU Time 1.5: " << chrono::duration_cast<chrono::milliseconds>(check1_2-check1).count() << " seconds" << endl;
  cout << "Initialize CPU Time 2: " << chrono::duration_cast<chrono::milliseconds>(check2-check1).count() << " seconds" << endl;
  cout << "Initialize CPU Time 3: " << chrono::duration_cast<chrono::milliseconds>(check3-check2).count() << " seconds" << endl;
  cout << "Initialize CPU Time 4: " << chrono::duration_cast<chrono::milliseconds>(check4-check3).count() << " seconds" << endl;
  cout << "Initialize CPU Time 5: " << chrono::duration_cast<chrono::milliseconds>(check4-check3).count() << " seconds" << endl;

}
