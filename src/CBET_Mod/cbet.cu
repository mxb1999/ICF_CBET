#include "CBET_Interface.hpp"
#include "parallelConfig.hpp"


//define CBET Gain Function
__global__ void 
cbetGain(CBETVars* constants, CBETArrs* arrays)
{

    int nbeams_cu = constants->nbeams_cu;
    int nrays_cu = constants->nrays_cu;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int beam = index / nrays_cu;
    if(beam >= (nbeams_cu))
    {
        return;
    }

    int raynum = index % nrays_cu;
    //imported constants
    int ncrossings_cu = constants->ncrossings_cu;
    int nx_cu = constants->nx_cu;
    int nz_cu = constants->nz_cu;
    int numstored_cu = constants->numstored_cu;
    double dx_cu = constants->dx_cu;
    double dz_cu = constants->dz_cu;
    double ncrit_cu = constants->ncrit_cu;
    double c_cu = constants->c_cu;
    double pi_cu = constants->pi_cu;
    double iaw_cu = constants->iaw_cu;
    double cs_cu = constants->cs_cu;
    double estat_cu = constants->estat_cu;
    double Ti_cu = constants->Ti_cu;
    double Te_cu = constants->Te_cu;
    double Z_cu = constants->Z_cu;
    double omega_cu = constants->omega_cu;
    double kb_cu = constants->kb_cu;
    double me_cu = constants->me_cu;

    //imported arrays
    double* i_b_cu = arrays->i_b_cu;
    double* i_b_new_cu = arrays->i_b_new_cu;

    double* W_cu = arrays->W_cu;

    double* x_cu = arrays->x_cu;

    double* z_cu = arrays->z_cu;
    double* W_new_cu = arrays->W_new_cu;
    double* dkx_cu = arrays->dkx_cu;
    double* dkz_cu = arrays->dkz_cu;
    double* dkmag_cu = arrays->dkmag_cu;
    double* uflow_cu = arrays->uflow_cu;

    int* ints_cu = arrays->ints_cu;
    double* eden_cu = arrays->eden_cu;
    int* boxes_cu = arrays->boxes_cu;
    int* numrays_cu = arrays->numrays_cu;
    int* present_cu = arrays->present_cu;    
    //iterate over each ray beam (excepting the last one)
    //each beam will be treated as a pump beam for those preceeding, as a seed beam for those following
    double constant1 = (pow(estat_cu,2.0))/(4*(1.0e3*me_cu)*c_cu*omega_cu*kb_cu*Te_cu*(1+3*Ti_cu/(Z_cu*Te_cu)));
    for(int m = 0; m < ncrossings_cu; m++)
    {
        int ix = vec4D_cu(boxes_cu, beam,raynum,m,0, nrays_cu, ncrossings_cu, 2);
        int iz = vec4D_cu(boxes_cu, beam,raynum,m,1, nrays_cu, ncrossings_cu, 2);
        if(!ix || !iz)
        {
            break;
        }
        if(!vec4D_cu(ints_cu, beam, raynum, m, 0, nrays_cu, ncrossings_cu, numstored_cu))
        {
            continue;
        }

        ix--;
        iz--;
        int index = 0;
        //find all rays that interact with raynum in the other beam(s) of higher order
        for(int q = 0; q < nbeams_cu; q++)
        {
            if(q == beam)
            {
                continue;
            }
            //find the number of rays interacted with
            int cnt = vec4D_cu(numrays_cu, beam,raynum,m,q, nrays_cu, ncrossings_cu, nbeams_cu);
            int* crossInd = new int[cnt]{0};

            for(int l = index; l < index+cnt; l++)//get the crossing index for all rays of beam q in (ix,iz)
            {
                for(int p = 0; p < ncrossings_cu;p++)
                {
                    
                    int ox = vec4D_cu(boxes_cu, q,l,p,0, nrays_cu, ncrossings_cu, 2)-1;
                    int oz = vec4D_cu(boxes_cu, q,l,p,1, nrays_cu, ncrossings_cu, 2)-1;
                    if(!ox || !oz)
                    {
                        break;
                    }
                    
                    if(ox == ix && oz == iz)
                    {
                        crossInd[l] = p;
                        break;
                    }
                    
                }
            }
            //double xprev = x_cu[ix];
            //double zprev = z_cu[ix];
            int nlim = (vec3D_cu(present_cu,beam,ix,iz, nx_cu, nz_cu) > cnt) ? cnt : vec3D_cu(present_cu,beam,ix,iz, nx_cu, nz_cu);
            for(int n = 0; n < nlim; n++)
            {
                
                int rayOther = vec4D_cu(ints_cu, beam, raynum, m, index + n, nrays_cu, ncrossings_cu, numstored_cu);
                int rayCross = crossInd[n];
                double mag1 = vec3D_cu(dkmag_cu, beam, raynum, m, nrays_cu, ncrossings_cu);
                double mag2 = vec3D_cu(dkmag_cu, q, rayOther, rayCross, nrays_cu, ncrossings_cu);
                if(mag2 < 1.0*dx_cu)
                {
                    continue;
                }
                double ne = vec2D_cu(eden_cu, ix,iz,nz_cu);
                double epsilon = 1.0-ne/ncrit_cu;
                double kmag = (omega_cu/c_cu)*sqrt(epsilon);

                double kx1 = kmag * vec3D_cu(dkx_cu, beam, raynum, m, nrays_cu, ncrossings_cu) / mag1;
                double kx2 = kmag * vec3D_cu(dkx_cu, q, rayOther, rayCross, nrays_cu, ncrossings_cu) / mag2;

                double kz1 = kmag * vec3D_cu(dkz_cu, beam, raynum, m, nrays_cu, ncrossings_cu) / mag1;
                double kz2 = kmag * vec3D_cu(dkz_cu, q, rayOther, rayCross, nrays_cu, ncrossings_cu) / mag2;
                double kiaw = sqrt(pow(kx2-kx1,2.0)+pow(kz2-kz1,2.0));

                double ws = kiaw*cs_cu;
                double omega1 = omega_cu;
                double omega2 = omega_cu;
                double eta = ((omega2-omega1)-(kx2-kx1)*vec2D_cu(uflow_cu,ix,iz,nz_cu))/(ws+1.0e-10);
                double efield2 = sqrt(8.*pi_cu*1.0e7*vec3D_cu(i_b_new_cu, q, rayOther, rayCross, nrays_cu, ncrossings_cu)/c_cu);   
                double P = (pow(iaw_cu,2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw_cu),2)*pow(eta,2));  
                double gain1 = constant1*pow(efield2,2)*(ne/ncrit_cu)*(1/iaw_cu)*P;               //L^-1 from Russ's paper
                double oldEnergy1 = vec3D_cu(W_new_cu, beam,raynum,m,nrays_cu, ncrossings_cu);
                double oldEnergy2 = ((beam < q) ? 1 : -1)*vec3D_cu(W_new_cu, q,rayOther,rayCross,nrays_cu, ncrossings_cu);
                double newEnergy1Mult = exp(oldEnergy2*mag1*gain1/sqrt(epsilon));
                vec3DW_cu(W_cu, beam, raynum, m, nrays_cu, ncrossings_cu,oldEnergy1);
                vec3DM_cu(W_new_cu, beam, raynum, m, nrays_cu, ncrossings_cu,newEnergy1Mult);
            }
        }
    }

    
}
__device__ void
elimCrossEffect(int currcross, int beamOther, int rayOther, int rayCross, CBETVars* constants, CBETArrs* arrays)
{
    int nbeams_cu = constants->nbeams_cu;
    int nrays_cu = constants->nrays_cu;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    int beam = index / nrays_cu;
    if(beam >= (nbeams_cu - 1))
    {
        return;
    }
    printf("Rectifying\n");
    int raynum = index % nrays_cu;
    //imported constants
    int ncrossings_cu = constants->ncrossings_cu;
    int nx_cu = constants->nx_cu;
    int nz_cu = constants->nz_cu;
    int numstored_cu = constants->numstored_cu;
    double dx_cu = constants->dx_cu;
    double dz_cu = constants->dz_cu;
    double ncrit_cu = constants->ncrit_cu;
    double c_cu = constants->c_cu;
    double pi_cu = constants->pi_cu;
    double iaw_cu = constants->iaw_cu;
    double cs_cu = constants->cs_cu;
    double estat_cu = constants->estat_cu;
    double Ti_cu = constants->Ti_cu;
    double Te_cu = constants->Te_cu;
    double Z_cu = constants->Z_cu;
    double omega_cu = constants->omega_cu;
    double kb_cu = constants->kb_cu;
    double me_cu = constants->me_cu;
    double constant1 = (pow(estat_cu,2.0))/(4*(1.0e3*me_cu)*c_cu*omega_cu*kb_cu*Te_cu*(1+3*Ti_cu/(Z_cu*Te_cu)));

    //imported arrays
    double* i_b_cu = arrays->i_b_cu;
    double* i_b_new_cu = arrays->i_b_new_cu;

    double* W_cu = arrays->W_cu;

    double* x_cu = arrays->x_cu;

    double* z_cu = arrays->z_cu;
    double* W_new_cu = arrays->W_new_cu;
    double* dkx_cu = arrays->dkx_cu;
    double* dkz_cu = arrays->dkz_cu;
    double* dkmag_cu = arrays->dkmag_cu;
    double* uflow_cu = arrays->uflow_cu;
    int m = currcross;
    int q = beamOther;
    int* ints_cu = arrays->ints_cu;
    double* eden_cu = arrays->eden_cu;
    int* boxes_cu = arrays->boxes_cu;
    int* numrays_cu = arrays->numrays_cu;
    int* present_cu = arrays->present_cu;    
    int ix = vec4D_cu(boxes_cu, beam,raynum,m,0, nrays_cu, ncrossings_cu, 2)-1;
    int iz = vec4D_cu(boxes_cu, beam,raynum,m,1, nrays_cu, ncrossings_cu, 2)-1;
    double mag1 = vec3D_cu(dkmag_cu, beam, raynum, m, nrays_cu, ncrossings_cu);
    double mag2 = vec3D_cu(dkmag_cu, q, rayOther, rayCross, nrays_cu, ncrossings_cu);
    double ne = vec2D_cu(eden_cu, ix,iz,nz_cu);
    double epsilon = 1.0-ne/ncrit_cu;
    double kmag = (omega_cu/c_cu)*sqrt(epsilon);

    double kx1 = kmag * vec3D_cu(dkx_cu, beam, raynum, m, nrays_cu, ncrossings_cu) / mag1;
    double kx2 = kmag * vec3D_cu(dkx_cu, q, rayOther, rayCross, nrays_cu, ncrossings_cu) / mag2;

    double kz1 = kmag * vec3D_cu(dkz_cu, beam, raynum, m, nrays_cu, ncrossings_cu) / mag1;
    double kz2 = kmag * vec3D_cu(dkz_cu, q, rayOther, rayCross, nrays_cu, ncrossings_cu) / mag2;
    double kiaw = sqrt(pow(kx2-kx1,2.0)+pow(kz2-kz1,2.0));

    double ws = kiaw*cs_cu;
    double omega1 = omega_cu;
    double omega2 = omega_cu;
    double eta = ((omega2-omega1)-(kx2-kx1)*vec2D_cu(uflow_cu,ix,iz,nz_cu))/(ws+1.0e-10);
    double efield2 = sqrt(8.*pi_cu*1.0e7*vec3D_cu(i_b_new_cu, q, rayOther, rayCross, nrays_cu, ncrossings_cu)/c_cu);   
    double P = (pow(iaw_cu,2)*eta)/(pow((pow(eta,2)-1.0),2)+pow((iaw_cu),2)*pow(eta,2));  
    double gain1 = constant1*pow(efield2,2)*(ne/ncrit_cu)*(1/iaw_cu)*P;               //L^-1 from Russ's paper
    double oldEnergy2 = ((beam < q) ? 1 : -1)*vec3D_cu(W_new_cu, q,rayOther,rayCross,nrays_cu, ncrossings_cu);
    double energyCorrection = exp(-1*oldEnergy2*mag1*gain1/sqrt(epsilon));
    vec3DM_cu(W_new_cu, beam, raynum, m, nrays_cu, ncrossings_cu,energyCorrection);
}
__global__ void 
cbetUpdate(CBETVars* constants, CBETArrs* arrays, int* cbetKill, int beam)
{
    
    int nbeams_cu = constants->nbeams_cu;
    int nrays_cu = constants->nrays_cu;
    int index = blockIdx.x*blockDim.x+threadIdx.x;
    if(index >= (nrays_cu))
    {
        return;
    }

    int raynum = index % nrays_cu;
    //imported constants
    int ncrossings_cu = constants->ncrossings_cu;
    int nx_cu = constants->nx_cu;
    int nz_cu = constants->nz_cu;
    int numstored_cu = constants->numstored_cu;
    double dx_cu = constants->dx_cu;
    double dz_cu = constants->dz_cu;
    double ncrit_cu = constants->ncrit_cu;
    double c_cu = constants->c_cu;
    double pi_cu = constants->pi_cu;
    double iaw_cu = constants->iaw_cu;
    double cs_cu = constants->cs_cu;
    double estat_cu = constants->estat_cu;
    double Ti_cu = constants->Ti_cu;
    double Te_cu = constants->Te_cu;
    double Z_cu = constants->Z_cu;
    double omega_cu = constants->omega_cu;
    double kb_cu = constants->kb_cu;
    double me_cu = constants->me_cu;

    //imported arrays
    double* i_b_cu = arrays->i_b_cu;
    double* i_b_new_cu = arrays->i_b_new_cu;

    double* W_cu = arrays->W_cu;

    double* x_cu = arrays->x_cu;

    double* z_cu = arrays->z_cu;
    double* W_new_cu = arrays->W_new_cu;
    double* dkx_cu = arrays->dkx_cu;
    double* dkz_cu = arrays->dkz_cu;
    double* dkmag_cu = arrays->dkmag_cu;
    double* uflow_cu = arrays->uflow_cu;
    printf("Rectifying\n");

    int* ints_cu = arrays->ints_cu;
    double* eden_cu = arrays->eden_cu;
    int* boxes_cu = arrays->boxes_cu;
    int* numrays_cu = arrays->numrays_cu;
    int* present_cu = arrays->present_cu;    
    //iterate over each ray beam (excepting the last one)
    //each beam will be treated as a pump beam for those preceeding, as a seed beam for those following
    double constant1 = (pow(estat_cu,2.0))/(4*(1.0e3*me_cu)*c_cu*omega_cu*kb_cu*Te_cu*(1+3*Ti_cu/(Z_cu*Te_cu)));
    int contact = 0;
    for(int m = 0; m < ncrossings_cu; m++)
    {
        int ix = vec4D_cu(boxes_cu, beam,raynum,m,0, nrays_cu, ncrossings_cu, 2);
        int iz = vec4D_cu(boxes_cu, beam,raynum,m,1, nrays_cu, ncrossings_cu, 2);
        if(!ix || !iz)
        {
            break;
        }

        if(!vec4D_cu(ints_cu, beam, raynum, m, 0, nrays_cu, ncrossings_cu, numstored_cu))
        {
            continue;
        }
        if(!contact)
        {
            contact = m;
        }
        ix--;
        iz--;
        int index = 0;
        //find all rays that interact with raynum in the other beam(s) of higher order
        for(int q = 0; q < nbeams_cu; q++)
        {
            if(q == beam)
            {
                continue;
            }
            //find the number of rays interacted with
            int cnt = vec4D_cu(numrays_cu, beam,raynum,m,q, nrays_cu, ncrossings_cu, nbeams_cu);
            int* crossInd = new int[cnt]{0};
            for(int l = index; l < cnt; l++)//get the crossing index for all rays of beam q in (ix,iz)
            {
                for(int p = 0; p < ncrossings_cu;p++)
                {
                    int ox = vec4D_cu(boxes_cu, q,l,p,0, nrays_cu, ncrossings_cu, 2)-1;
                    int oz = vec4D_cu(boxes_cu, q,l,p,1, nrays_cu, ncrossings_cu, 2)-1;
                    if(ox == ix && oz == iz)
                    {
                        crossInd[l] = p;
                        break;
                    }
                }
                int killed = vec2D_cu(cbetKill, q,l,nrays_cu);
                if(killed && ((killed-1) < crossInd[l]))
                {
                    elimCrossEffect(m, q,l,crossInd[l], constants, arrays);
                }
            }
            double prevIntensitySeed = vec3D_cu(i_b_new_cu, beam,raynum,m,nrays_cu, ncrossings_cu);
            double intensityMultSeed = (1.0 - (vec3D_cu(W_new_cu, beam,raynum,m,nrays_cu, ncrossings_cu)/vec3D_cu(W_cu, beam,raynum,m,nrays_cu, ncrossings_cu)));
            double dI1 = (-1.0*intensityMultSeed*prevIntensitySeed);
            int kill = 0;
            int transfer = 0;
            if(abs(dI1) > prevIntensitySeed && dI1 < 0)
            {
                kill = 1;
                printf("killed at: %d\n", m-contact);
            }
            if(prevIntensitySeed <= 0)
            {
                kill = 1;
                prevIntensitySeed = 0;
            }
            //printf("%e\n", vec3D_cu(W_new_cu, beam,raynum,m,nrays_cu, ncrossings_cu)/vec3D_cu(W_cu, beam,raynum,m,nrays_cu, ncrossings_cu));

            for(int l = m; l < ncrossings_cu; l++)
            {
                double newVal;
                if(!kill)
                {
                    vec3DI_cu(i_b_new_cu, beam,raynum,l,nrays_cu, ncrossings_cu,dI1);
                    newVal = prevIntensitySeed+dI1;
                }else
                {
                    vec3DW_cu(i_b_new_cu, beam,raynum,l,nrays_cu, ncrossings_cu,0.0);
                    newVal = 0.0;
                }
            }
            if(kill)
            {
                vec2DW_cu(cbetKill, beam, raynum, nrays_cu, m+1);
            }
            
            
        }
    }

    
}

void launchCBETKernel()
{
    printf("CBET\n");
    initArrays();
    CBETVars* vars = new_cbetVars();
    CBETArrs* arrays = new_cbetArrs();
    int* cbetKill;
    cudaMallocManaged(&cbetKill, sizeof(int)*RAYS);
    for(int i = 0; i < RAYS; i++)
    {
        cbetKill[i] = 0;
    }
    int B1 = nrays*nbeams/256+1;
    int B2 = nrays/256+1;

    for(int i = 0; i < 1; i++)
    {
        cbetGain<<<B1, 256>>>(vars, arrays);
        cudaDeviceSynchronize();
        for(int i = 0; i < nbeams;i++)
        {
            cbetUpdate<<<B2, 256>>>(vars, arrays,cbetKill, i);
            cudaDeviceSynchronize();
        }
        
    }

    printf("%e\n", cs);
}

