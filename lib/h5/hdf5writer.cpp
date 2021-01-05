#include "implSim.hpp"
#include <H5Cpp.h>


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
    double* copy = new double[prod];//declare a new 1D array to be written to temporarily
    //Copies arr to copy using a 1D indexing scheme based on the number of dimensions
    switch(dimnum)
    {
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
        arrData->write(copy, PredType::NATIVE_DOUBLE);
        break;
      }
      //3D array
      case 3:
      {
        double*** arrcast = (double***)arr;
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
    }
    arrData->close();
    arrspace->close();
    //cout << name << endl;
    delete [] copy;
  }else
  {
    //Same as above, just with the int datatype instead of double
    int fill = 0;
    DSetCreatPropList plist;
    plist.setFillValue(PredType::NATIVE_INT, &fill);
    DataSpace arrspace(dimnum, dims);
    DataSet* arrData = new DataSet(store->createDataSet(name, PredType::NATIVE_INT, arrspace,plist));
    int* copy = new int[prod];
    switch(dimnum)
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
        arrData->write(copy, PredType::NATIVE_INT);
        break;
      }
    }
    arrData->close();
    arrspace.close();

  }
}
void writePlotArrays()
{
    //Arrays to be plotted via Python
    i_b_newplot = new double*[nx];
    edepplot = new double*[nz];
    edenplot = new double*[nz];
    raytrace = new double*[nx];
    ib_orig = new double*[nx];
    anyInt = new int*[nx];
    perturbation = new double*[nx];
    for(int i = 0; i < nz;i++)
    {
      i_b_newplot[i] = new double[nz];
      edepplot[i] = new double[nx]{0.0};
      edenplot[i] = new double[nx]{0.0};
      perturbation[i] = new double[nz];
      raytrace[i] = new double[nz]{0.0};
      ib_orig[i] = new double[nz]{0.0};
      anyInt[i] = new int[nz];
    }
    for(int i = 0; i < nx;i++)
    {
      for(int j = 0; j < nz;j++)
      {
        edenplot[i][j] = eden[i][j]/ncrit;
        i_b_newplot[i][j] = 8.53e-10*sqrt(fmax(1.0e-10,i_b_new[0][i][j])+fmax(1.0e-10,i_b_new[1][i][j]))*(1.053/3.0);
        perturbation[i][j] = fmax(W_new[0][i][j], W_new[1][i][j]) - sqrt(1.0-eden[i][j]/ncrit)/double(rays_per_zone);
        raytrace[i][j] = ((present[0][i][j]+present[1][i][j]) != 0);
        ib_orig[i][j] = 8.53e-10*sqrt(edep[0][i][j]+edep[1][i][j]+1.0e-10)*(1.053/3.0);
        if(intersections[i][j] > 0)
        {
          anyInt[i][j] = 1;
        }
        for(int m = 0;m<nbeams;m++)
        {
          edepplot[i][j]+=edep[m][i][j];
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

  //Output arrays to be plotted in Python using included script'

  //Core output arrays
  writeArr(ib_orig, 0, store, "/Original_Field", 2, new int[2]{nz,nx});//original electric field
  writeArr(i_b_newplot, 0, store, "/New_Field", 2, new int[2]{nz,nx});//Post-CBET electric field
  writeArr(x, 0, store, "/x", 1, new int[1]{nx});//x coordinates
  writeArr(z, 0, store, "/z", 1, new int[1]{nz});//z coordinates
  writeArr(eden, 0, store, "/eden", 2, new int[2]{nz,nx});//electron density gradient
  writeArr(edenplot, 0, store, "/eden_ncrit", 2, new int[2]{nz,nx});//normalized electron density gradient
  writeArr(machnum, 0, store, "/machnum", 2, new int[2]{nx,nz});//implosion velocity relative to Mach 1
  writeArr(perturbation, 0, store, "/density_perturbation", 2, new int[2]{nx,nz});//electron density pertubation (or indicative of it)
  writeArr(edepplot, 0, store, "/total_intensity", 2, new int[2]{nz,nx});//total intensity deposited

  //Arrays useful for debugging
  //writeArr(orderplot1, 0, store, "/updated", 2, new int[2]{nz,nx});//stores all updated ray locations
  writeArr(raytrace, 0, store, "/raydist", 2, new int[2]{nz,nx});//all ray paths found
  writeArr(intersections, 0, store, "/intersections", 2, new int[2]{nz,nx});//all ray paths found
  writeArr(W_new[0], 0, store, "/new_energy1", 2, new int[2]{nz,nx});//energy deposited, CBET multiplier for beam 1
  writeArr(W_new[1], 0, store, "/new_energy2", 2, new int[2]{nz,nx});//energy deposited, CBET multiplier for beam 2
  writeArr(W[0], 0, store, "/original_energy1", 2, new int[2]{nz,nx});//original energy deposited
  writeArr(W[1], 0, store, "/original_energy2", 2, new int[2]{nz,nx});//original energy deposited
  writeArr(anyInt, 1, store, "/nonZero", 2, new int[2]{nx,nz});//stores any location where a ray has been as one
  writeArr(i_b_new[0], 0, store, "/beam1_intensity", 2, new int[2]{nz,nx});//beam 1 intensity post-CBET
  writeArr(i_b_new[1], 0, store, "/beam2_intensity", 2, new int[2]{nz,nx});//beam 2 intensity post-CBET


  store->close();//close hdf file
  if(printUpdates)
  {
    cout << "Write Finished" << endl;
    cout << "File Closed" << endl;
  }

}
