#include "include/implSim.h"
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
void writeArr(void* arr, int type, H5File* store, string name, int dimnum, int* dimlist)
{

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
  if(type == 0)
  {
    //double fill = 0.0;
    //DSetCreatPropList plist;
    //plist.setFillValue(PredType::NATIVE_DOUBLE, &fill);
    DataSpace* arrspace = new DataSpace(dimnum, dims);
    DataSet* arrData = new DataSet(store->createDataSet(name, PredType::NATIVE_DOUBLE, *arrspace));
    arrData->write(arr, PredType::NATIVE_DOUBLE);

    double* copy = new double[prod];
    switch(dimnum)
    {
      case 1:
      {
        arrData->write(arr, PredType::NATIVE_DOUBLE);
        break;
      }
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
void updateH5()
{
  cout << "Started HDF5 Write" << endl;
  std::string name = "Output/implSim.hdf";
  static H5File* store = new H5File(name, H5F_ACC_TRUNC);
  cout << "File Opened" << endl;
  cout << "Starting Write..." << endl;
  writeArr(x, 0, store, "x (cm)", 1, new int[1]{nx});
  writeArr(z, 0, store, "z (cm)", 1, new int[1]{nz});
  writeArr(eden, 0, store, "eden", 2, new int[2]{nx,nz});
  writeArr(orderplot1, 0, store, "order1", 2, new int[2]{nx,nz});
  writeArr(orderplot2, 0, store, "order2", 2, new int[2]{nx,nz});
  writeArr(edenplot, 0, store, "eden_ncrit", 2, new int[2]{nx,nz});

  writeArr(edepplot, 0, store, "total_intensity", 2, new int[2]{nx,nz});
  writeArr(i_b[0], 0, store, "i_b1", 2, new int[2]{nx,nz});
  writeArr(i_b[1], 0, store, "i_b2", 2, new int[2]{nx,nz});
  writeArr(i_b1Error, 0, store, "i_b1Error", 2, new int[2]{nx,nz});
  writeArr(i_b2Error, 0, store, "i_b2Error", 2, new int[2]{nx,nz});
  writeArr(i_b_new[0], 0, store, "i_b2_new", 2, new int[2]{nx,nz});
  writeArr(i_b_new[1], 0, store, "i_b1_new", 2, new int[2]{nx,nz});
  writeArr(i_b_newplot, 0, store, "i_b_newplot", 2, new int[2]{nx,nz});
  writeArr(i_bplot, 0, store, "i_bplot", 2, new int[2]{nx,nz});
  cout << "Write Finished" << endl;
  store->close();
  cout << "File Closed" << endl;

}
