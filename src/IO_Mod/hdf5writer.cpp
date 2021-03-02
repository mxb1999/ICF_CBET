#include "io_interface.hpp"


using namespace H5;


//determines the type based off of a user string
int determineType(string type)
{
  if(type.compare("double*") == 0)
  {
    return 0;
  }else if(type.compare("double**") == 0)
  {
    return 1;
  }else if(type.compare("double***") == 0)
  {
    return 2;
  }else if(type.compare("double****") == 0)
  {
    return 3;
  }else if(type.compare("int*") == 0)
  {
    return 4;
  }else if(type.compare("int**") == 0)
  {
    return 5;
  }else if(type.compare("int***") == 0)
  {
    return 6;
  }else if(type.compare("int****") == 0)
  {
    return 7;
  }else
  {
    return -1;
  }
};
//function for writing a dataset to HDF5: arr is the file to be written, type denotes whether it is a double (0) or int (1) (can be expanded), store is the file being written to,
//name is the name of the dataset to be created for reference, dimnum denotes the number of dimensions in the set, dimlist is the size of those dimensions
void writeArr(void* arr, int type, H5File* store, string name, int dimnum, int* dimlist)
{
  if(!arr)
  {
    return;
  }
  //Finds the overall size of the dataset
  int prod = 1;
  for(int i = 0; i < dimnum;i++)
  {
    prod *= dimlist[i];
  }
  hsize_t dims[dimnum];
  for(int i = 0; i < dimnum; i++)
  {
    dims[i] = dimlist[i];
  }
  //If the dataset contains doubles
  if(type == 0)
  {
    DataSpace* arrspace = new DataSpace(dimnum, dims);    //Allocate the memory in an array space in the .hdf file
    DataSet* arrData = new DataSet(store->createDataSet(name, PredType::NATIVE_DOUBLE, *arrspace));    //write data in the case that arr is simply an int
    //double* copy = new double[prod];//declare a new 1D array to be written to temporarily
    //Copies arr to copy using a 1D indexing scheme based on the number of dimensions
    double* arrcast = (double*)arr;

    /*switch(dimnum)
    {
          double* arrp;
    double* copyp;
      //1D array
      case 1:
      {
        arrData->write(arr, PredType::NATIVE_DOUBLE);
        break;
      }
      //2D array
      case 2:
      {
        double** arrcast = (double**)arr;
        for(int i = 0; i < dimlist[0];i++)
        {
          for(int j = 0; j < dimlist[1];j++)
          {
            copy[i*dimlist[1]+j] = arrcast[i][j];
          }
        }
        break;
      }
      //3D array
      case 3:
      {
        for(int i = 0; i < dimlist[0];i++)
        {
          for(int j = 0; j < dimlist[1];j++)
          {
            for(int m = 0; m < dimlist[2];m++)
            {
              copy[(i*dimlist[1]+j)*dimlist[2]+m] = arrcast[i][j][m];
            }
          }
        }
        arrData->write(copy, PredType::NATIVE_DOUBLE);
        break;
      }
      //4D array
      case 4:
      {
        double**** arrcast = (double****)arr;
        for(int i = 0; i < dimlist[0];i++)
        {
          for(int j = 0; j < dimlist[1];j++)
          {
            for(int m = 0; m < dimlist[2];m++)
            {
              for(int n = 0; n < dimlist[3];n++)
              {
                copy[((i*dimlist[1]+j)*dimlist[2]+m)*dimlist[3]+n] = (double)arrcast[i][j][m][n];
              }
            }
          }
        }
        arrData->write(copy, PredType::NATIVE_DOUBLE);
        break;
      }
    }*/
    arrData->write(arrcast, PredType::NATIVE_DOUBLE);
    arrData->close();
    arrspace->close();
    //cout << name << endl;
  }else
  {
    //Same as above, just with the int datatype instead of double
    int fill = 0;
    DSetCreatPropList plist;
    plist.setFillValue(PredType::NATIVE_INT, &fill);
    DataSpace arrspace(dimnum, dims);
    DataSet* arrData = new DataSet(store->createDataSet(name, PredType::NATIVE_INT, arrspace,plist));
    int* arrcast = (int*)arr;
    //int* copy = new int[prod];
    /*switch(dimnum)
    {
      case 1:
      {
        arrData->write(arr, PredType::NATIVE_INT);
        break;
      }
      case 2:
      {
        int** arrcast = (int**)arr;
        for(int i = 0; i < dimlist[0];i++)
        {
          for(int j = 0; j < dimlist[1];j++)
          {
            copy[i*dimlist[1]+j] = arrcast[i][j];
          }
        }
        arrData->write(copy, PredType::NATIVE_INT);
        break;
      }
      case 3:
      {
        int*** arrcast = (int***)arr;
        for(int i = 0; i < dimlist[0];i++)
        {
          for(int j = 0; j < dimlist[1];j++)
          {
            for(int m = 0; m < dimlist[2];m++)
            {
              copy[(i*dimlist[1]+j)*dimlist[2]+m] = arrcast[i][j][m];
            }
          }
        }
        arrData->write(copy, PredType::NATIVE_INT);
        break;
      }
      case 4:
      {
        int**** arrcast = (int****)arr;
        for(int i = 0; i < dimlist[0];i++)
        {
          for(int j = 0; j < dimlist[1];j++)
          {
            for(int m = 0; m < dimlist[2];m++)
            {
              for(int n = 0; n < dimlist[3];n++)
              {
                copy[((i*dimlist[1]+j)*dimlist[2]+m)*dimlist[3]+n] = (int)arrcast[i][j][m][n];
              }
            }
          }
        }
        break;
      }
    }
    */
    arrData->write(arrcast, PredType::NATIVE_INT);
    arrData->close();
    arrspace.close();
  }
}
void writePlotArrays()
{
    //Arrays to be plotted via Python
    i_b_newplot = new double[nx*nz]{0.0};
    edepplot = new double[nx*nz]{0.0};
    edenplot = new double[nx*nz]{0.0};
    raytrace = new int[nx*nz]{0};
    ib_orig = new double[nx*nz]{0.0};
    anyInt = new int[nx*nz];
    perturbation = new double[nx*nz];
    for(int i = 0; i < nx;i++)
    {
      for(int j = 0; j < nz;j++)
      {
        vec2DW(edenplot,i,j,nz, vec2D(eden,i,j,nz)/ncrit);

        if(calcCBET)
        {
        vec2DW(i_b_newplot,i,j,nz, 8.53e-10*sqrt(fmax(1.0e-10,vec3D(i_b_new,0,i,j,nx,nz))+fmax(1.0e-10,vec3D(i_b_new,1,i,j,nx,nz)))*(1.053/3.0));
        vec2DW(perturbation,i,j,nz, fmax(vec3D(W_new,0,i,j,nx,nz), vec3D(W_new,1,i,j,nx,nz)) - sqrt(1.0-vec2D(eden,i,j,nz)/ncrit)/double(rays_per_zone));

        }
        vec2DW(raytrace,i,j,nz, vec3D(present,1,i,j,nx,nz)+vec3D(present,0,i,j,nx,nz));//+vec3D(present,0,i,j,nx,nz);
        vec2DW(ib_orig,i,j,nz, 8.53e-10*sqrt(vec3D(edep,0,i,j,nx,nz)+vec3D(edep,1,i,j,nx,nz)+1.0e-10)*(1.053/3.0));
        if(vec2D(intersections,i,j,nz)> 0)
        {
          vec2DW(anyInt,i,j,nz, 1);
        }
        for(int m = 0;m<nbeams;m++)
        {
          vec2DI(edepplot,i,j,nz,vec3D(edep,m,i,j,nx,nz));
        }
      }
    }
}
void updateH5()
{
  std::string name = "output/implSim.hdf";
  if(printUpdates)
  {
    cout << "Initializing plot arrays..." << endl;
  }
  writePlotArrays();
  static H5File* store = new H5File(name, H5F_ACC_TRUNC);
  if(printUpdates)
  {
    cout << "File Opened" << endl;
    cout << "Starting Write..." << endl;
  }
  double** multArr = new double*[nx];
  for(int i = 0; i < nx;i++)
  {
    multArr[i] = new double[nz];
  }
  for(int i = 0; i < nrays;i++)
  {
    for(int j = 0; j < ncrossings; j++)
    {
      int boxx = vec4D(boxes, 0,i,j,0,nrays,ncrossings,2);
      int boxz = vec4D(boxes, 0,i,j,1,nrays,ncrossings,2);
      if(!boxx || !boxz)
      {
        break;
      }
      boxx--;
      boxz--;
      multArr[boxx][boxz] = mult[i*ncrossings+j];
    }
  }
  //Output arrays to be plotted in Python using included script'
  printf("Check 1\n");
  fflush(stdout);
  //Core output arrays
  if(calcCBET)
  {
  writeArr(ib_orig, 0, store, "/Original_Field", 2, new int[2]{nx,nz});//original electric field
  writeArr(i_b_newplot, 0, store, "/New_Field", 2, new int[2]{nx,nz});//Post-CBET electric field
  writeArr(perturbation, 0, store, "/density_perturbation", 2, new int[2]{nx,nz});//electron density pertubation (or indicative of it)
  writeArr(W_new, 0, store, "/new_energy1", 2, new int[2]{nx,nz});//energy deposited, CBET multiplier for beam 1

  writeArr(W_new+(nx*nz), 0, store, "/new_energy2", 2, new int[2]{nx,nz});//energy deposited, CBET multiplier for beam 2
  writeArr(W, 0, store, "/original_energy1", 2, new int[2]{nx,nz});//original energy deposited
  writeArr(W+(nx*nz), 0, store, "/original_energy2", 2, new int[2]{nx,nz});//original energy deposited
   
  writeArr(anyInt, 1, store, "/nonZero", 2, new int[2]{nx,nz});//stores any location where a ray has been as one
  writeArr(i_b, 0, store, "/beam1_intensity", 2, new int[2]{nx,nz});//beam 1 intensity post-CBET
  writeArr(i_b+(nx*nz), 0, store, "/beam2_intensity", 2, new int[2]{nx,nz});//beam 2 intensity post-CBET
  writeArr(multArr, 0, store, "/gain1", 2, new int[2]{nx,nz});//beam 2 intensity post-CBET
  writeArr(gain2arr, 0, store, "/gain2", 2, new int[2]{nx,nz});//beam 2 intensity post-CBET
  writeArr(mag, 0, store, "/mag", 2, new int[2]{nx,nz});//beam 2 intensity post-CBET
  writeArr(u_flow, 0, store, "/u_flow", 2, new int[2]{nx,nz});//beam 2 intensity post-CBET
  }
  writeArr(edepplot, 0, store, "/edep",2, new int[2]{nx,nz});
  writeArr(x, 0, store, "/x", 1, new int[1]{nx});//x coordinates
  writeArr(z, 0, store, "/z", 1, new int[1]{nz});//z coordinates
  writeArr(eden, 0, store, "/eden", 2, new int[2]{nx,nz});//electron density gradient
  writeArr(edenplot, 0, store, "/eden_ncrit", 2, new int[2]{nx,nz});//normalized electron density gradient
  writeArr(machnum, 0, store, "/machnum", 2, new int[2]{nx,nz});//implosion velocity relative to Mach 1
  writeArr(edepplot, 0, store, "/total_intensity", 2, new int[2]{nx,nz});//total intensity deposited

  //Arrays useful for debugging
  //writeArr(orderplot1, 0, store, "/updated", 2, new int[2]{nx,nz});//stores all updated ray locations
  writeArr(raypath, 1, store, "/raypath", 2, new int[2]{nx,nz});//total intensity deposited

  writeArr(raytrace, 1, store, "/raydist", 2, new int[2]{nx,nz});//all ray paths found

  writeArr(intersections, 1, store, "/intersections", 2, new int[2]{nx,nz});//all ray paths found
      printf("Check 2\n");
  fflush(stdout);
 
  store->close();//close hdf file
  if(printUpdates)
  {
    cout << "Write Finished" << endl;
    cout << "File Closed" << endl;
  }

}
