#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
//#include "customMath.hpp"
#include <limits>
//#include "implSim.hpp"
#include <algorithm>
#include <functional>
#include <iostream>
#include <cstring>
//implementations of custom math functions


//standard binary search algorithm, searches given array in log(n) time
int binSearch(double tgt, double* in, int lIndex, int rIndex)
{
  double left = in[lIndex];
  int midIndex = lIndex+(rIndex-lIndex)/2;
  double mid = in[midIndex];
  if((midIndex - lIndex == 1 && tgt <= mid) || rIndex - lIndex == 1)
  {
    return lIndex;
  }
  if(tgt <= mid)
  {
    return binSearch(tgt, in, lIndex,midIndex);
  }
  if(tgt > mid)
  {
    return binSearch(tgt, in, midIndex,rIndex);
  }

  return -1;
};

//returns the interpolated value from the given sorted arrays
double interp(double* xArr, double* yArr, double target, int xSize)
{
  if(target <= xArr[0])
  {
      return yArr[0];
  }
  if(target >= xArr[xSize-1])
  {
      return yArr[xSize-1];
  }
  int index = binSearch(target, xArr, 0, xSize - 1);

  double m = (yArr[index+1]-yArr[index])/(xArr[index+1]-xArr[index]);
  double b = yArr[index]-m*xArr[index];
  return m*target+b;
};
//used for interpolating all values of a given array, relies on previous interp function
//xsize is the length of "xArr" and "yArr", size is the length of the target array
double* interpArr(double* xArr, double* yArr, double* target, int xsize, int size)
{
  double* result = new double[size];
  for(int i = 0; i < size; i++)
  {
    result[i] = interp(xArr, yArr, target[i],xsize);

  }
  return result;
}
//fills an array of size num with evenly spaced double values between start and stop
void span(double* target, double start, double stop, int num)
{
  float increment = (stop-start)/(num - 1);
  for(int i = 0; i < num; i++)
  {
    target[i] = start + (increment * i);
  }
}

int partition(double* arr, int min, int max)
{
  double pivot = arr[max-1];
  int i = min - 1;
  for(int j = min; j < max; j++)
  {
    if(arr[j] < pivot)
    {
      i++;
      double temp = arr[i];
      arr[i] = arr[j];
      arr[j] = temp;
    }
  }
  double temp = arr[i+1];
  arr[i+1] = arr[max-1];
  arr[max-1] = temp;
  return i+1;
}
void sort(double* arr, int min, int max)
{
  if(min < max)
  {
    int partInd = partition(arr, min, max);

    sort(arr, min, partInd-1);
    sort(arr, partInd+1, max);
  }
}

//find the median index of an unsorted array
double median(double* arr, int N)
{
  int cnt = 0;
  for(int i = 0; i < N; i++)
  {
    if(arr[i] == arr[i])
    {
      cnt++;
    }
  }
  double copy[cnt];
  cnt = 0;
  for(int i = 0; i < N; i++)
  {
    if(arr[i] == arr[i])
    {
      copy[cnt] = arr[i];
      cnt++;
    }
  }
  sort(copy, 0, cnt);
  return copy[cnt/2];
}
/*
int main(int argc, char** argv)
{
  double tester[] = {2.32,3.3243, 34234.34, 324, 9645.0, 0.0032};
  sort(tester, 0, 6);
  for(int i = 0; i < 6; i++)
  {
    printf("%f\n", tester[i]);
  }
  return 0;
}*/