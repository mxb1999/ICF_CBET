#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;
//custom math functions used in simulation code
int binSearch(double tgt, double* in, int index, int size);
double interp(double* xArr, double* yArr, double target, int xSize);
double* interpArr(double* xArr, double* yArr, double* target, int xsize, int size);
bool areEqual(double a, double b);
int compareDub(double a, double b);
void span(double* target, double start, double stop, int num);
