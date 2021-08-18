#include "implSim.hpp"

/*Algorithm for mapping rays to field
 For each ray of each beam -> sample intensity and set to grid location (Beam Grid)
 Initialize as 0
 Update Field data w/ += beamGrid**2
 Finalize field data as sqrt(field data)
*/
double* direct_map(double* cbetAmp)
{
    double* fieldData = (double*)calloc(sizeof(double), GRID);
    double* beamGrid = (double*)calloc(sizeof(double), GRID);
    for(int i = 0; i < nbeams; i++)
    {
        for(int j = 0; j < nrays; j++)
        {
            for(int k = 0; k < ncrossings; k++)
            {
                int ix = vec4D(boxes, i, j, k, 0, nrays, ncrossings, 2);
                int iz = vec4D(boxes, i, j, k, 1, nrays, ncrossings, 2);
                if(!ix || !iz)
                {
                    break;
                }
                ix--;
                iz--;
                double amp = vec3D(cbetAmp, i, j, k, nrays, ncrossings)/1e14;
                if(ix == 0)
                {
                    //printf("Hello %e\n", amp);
                }
                vec2DW(beamGrid, ix, iz, nz, amp);
            }
        }
        for(int j = 0; j < nx; j++)
        {
            for(int k = 0; k < nz; k++)
            {
                double gridentry = vec2D(beamGrid, j, k, nz);

                vec2DI(fieldData, j, k, nz, pow(gridentry, 2));
                vec2DW(beamGrid, j, k, nz, 0.0);
            }
        }
    }
    for(int j = 0; j < nx; j++)
    {
        for(int k = 0; k < nz; k++)
        {
            double fieldentry = vec2D(fieldData, j, k, nz);
            vec2DW(fieldData, j, k, nz, sqrt(fieldentry));
        }
    }
    free(beamGrid);
    return fieldData;
}