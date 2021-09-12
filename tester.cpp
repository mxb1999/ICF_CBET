#include <cmath>
#include <stdio.h>
const int nbeams = 1;
const int nrays = 1;
const int ncrossings = 11;
#define vec3D(p, i, j, k, s2, s3) p[((i)*(s2) + j)*(s3) + k]
#define vec3DW(p, i, j, k, s2, s3, val) p[((i)*(s2) + j)*(s3) + k] = val

double limitEnergy(double multiplier_acc, double i0, double currMax, int beam, int ray, int crossing, double* maxChange, double* i_b)
{
  double i_prev = vec3D(i_b, beam, ray, crossing, nrays, ncrossings);
  double i_curr = i0*multiplier_acc;
  double fractionalChange = fabs(i_curr-i_prev)/i_prev;
  printf("(%d %e) ",crossing, fractionalChange);
  *maxChange = fmax(fractionalChange, *maxChange);
  double correction = 0.0;
  if(fractionalChange > currMax)
  {
    int sign = (i_curr - i_prev > 0) ? 1 : -1;
    correction = 1 + currMax*sign;
    i_curr = i_prev*correction;
  }
  //printf("{%d %e}, ", crossing, (i_curr-i0*multiplier_acc)/1e14);
  return i_curr;
}
void updateIntensities(int iteration, double currMax, double* i_b, double* wMult, double* maxChange)
{
  for(int i = 0; i < nbeams; i++)
  {
    for(int j = 0; j < nrays; j++)
    {
      //printf("%d: ", j + 1 + 16*i);
      double i0 = vec3D(i_b, i, j, 0, nrays, ncrossings);
      double multAcc = 1.0;
      for(int k = 1; k < ncrossings; k++)
      {
        double mult = vec3D(wMult, i, j, k, nrays, ncrossings);
        multAcc *= mult;
        double new_intensity = limitEnergy(multAcc, i0, currMax, i, j, k, maxChange, i_b);
        //printf("{%d %e}, ", k+1, new_intensity);
        vec3DW(i_b, i, j, k, nrays, ncrossings, i0*multAcc);
      }
      printf("\n");
    }
  }
}
#define ABSLIM 1000
#define maxIncr 0.2
double limitMultiplier(double mult, int beam, int ray, int crossing, int iteration, double* stepSize)
{
  double stepRatio = vec3D(stepSize, beam, ray, crossing, nrays, ncrossings);
  double multPerStep = (1-mult)/stepRatio;
  double multiplierLim = maxIncr*iteration;
  multPerStep = (multPerStep > multiplierLim) ? multiplierLim : multPerStep;
  multPerStep = (multPerStep < -1*multiplierLim) ? -1*multiplierLim : multPerStep;
  mult = 1 - multPerStep*stepRatio;
  mult = (mult > ABSLIM) ? ABSLIM : mult;
  mult = (mult < 1/ABSLIM) ? 1/ABSLIM : mult;

  return mult;
};
int main(int argc, char** argv)
{
    double intensity[] ={0.0400653,	0.0400653,	0.0400653,	0.0400653,	0.0400653,	0.0400653,	0.0400653,	0.0400653,	0.0400653,	0.0400653,	0.0400653,	0.0400653};
    double mult[] = {1,	1,	1.0148592,	1.2551054,	35.561184,	44.299797,	2.6886756,	1.0001408,	1,	1,	1};
    double maxChange = 0;
    double stepSize[] = {0.410007412552480,	1.00029611586112,	1.00121053451488,	1.00277308212409,	1.00497968173124,	1.00782453976514,	1.01130304282651,	1.01541347430165,	0.918236402595692,	0.101905443875530,	0.555452308394799};
    //updateIntensities(1, 0.2, intensity, mult, &maxChange);
    for(int i = 0; i < ncrossings; i++)
    {
      printf("%f, ", limitMultiplier(mult[i], 0, 0, i, 1, stepSize));
    }
    printf("%f\n", maxChange);
    //for(int i = 0; i < ncrossings; i++)
    //{
        //printf("%f, ", intensity[i]);
    //}
}