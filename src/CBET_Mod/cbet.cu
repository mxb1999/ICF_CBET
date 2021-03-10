#include "CBET_Interface.hpp"
#include "parallelConfig.hpp"



//define CBET Gain Function
__global__ void cbetGain(CBETVars& constants, CBETArrs& arrays)
{
    int nbeams_cu = constants.nbeams_cu;
    int nrays_cu = constants.nrays_cu;
    int ncrossings_cu = constants.ncrossings_cu;
    int nx_cu = constants.nx_cu;
    int nz_cu = constants.nz_cu;
    double* i_b_cu = arrays.i_b_cu;
    double* W_cu = arrays.W_cu;
    double* W_new_cu = arrays.W_new_cu;
    double* dkx_cu = arrays.dkx_cu;
    double* dkz_cu = arrays.dkz_cu;
    double* dkmag_cu = arrays.dkmag_cu;
    int* ints_cu = arrays.ints_cu;

    //iterate over each ray beam (excepting the last one)
    //each beam will be treated as a pump beam for those preceeding, as a seed beam for those following
    for(int i = 0; i < nbeams_cu-1;i++)
    {
        //iterate over all rays in the beam
        for(int j = 0; j < nrays_cu; j++)
        {
            for(int m = 0; m < ncrossings_cu; m++)
            {
                int ix = vec4D(boxes, i,j,m,0, nrays, ncrossings, 2);
                int iz = vec4D(boxes, i,j,m,1, nrays, ncrossings, 2);
                if(!ix || !iz || !vec4D(ints, i, j, m, 0, nrays_cu, ncrossings_cu, numstored_cu,0))
                {
                    break;
                }
                ix--;
                iz--;
                for(int i = 0; i < )
            }
        }
    }
}
//define CBET Update function
__global__ void cbetUpdate(CBETVars& constants, CBETArrs& arrays)
{

}