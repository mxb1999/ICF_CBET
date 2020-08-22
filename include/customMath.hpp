#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;
//custom math functions used in simulation code
int binSearch(double tgt, double* in, int index, int size);//Binary search for interp function, returns the index of the element such that in[index] <= tgt <= in[index+1]
double interp(double* xArr, double* yArr, double target, int xSize);//Returns an interpolated value based on target using xArr and yArr as a piecewise function
double* interpArr(double* xArr, double* yArr, double* target, int xsize, int size);//Returns an array filled with interpolated values using xArr and yArr as a piecewise function
bool areEqual(double a, double b);//determines if the two doubles are equal, not used but can be more convenient at times
int compareDub(double a, double b);//Compares doubles
void span(double* target, double start, double stop, int num);//Fills target with evenly spaced doubles running from start to stop
