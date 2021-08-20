#include "Trace_interface.hpp"
#include "traceVars.hpp"
void fillTraceArrays()
{

    if(cudaCalc)
    {
        //marked = new int[GRID*nbeams*numstored]{0}; //nx nz nrays nbeams
        cudaMallocManaged(&dedendx, sizeof(double)*GRID);
        //dedendx = new double[GRID]; //nx nz
        cudaMallocManaged(&dedendz, sizeof(double)*GRID);
        //dedendz = new double[GRID]; //nx nz
        cudaMallocManaged(&x, sizeof(double)*nx);
        //x = new double[nx]{0.0}; //nx nz
        cudaMallocManaged(&z, sizeof(double)*nz);
        //z = new double[nz]{0.0}; //nx nz
        cudaMallocManaged(&eden, sizeof(double)*GRID);
        //eden = new double[GRID]; //nx nz
        cudaMallocManaged(&marked, sizeof(int)*GRID*numstored*nbeams);
        cudaMallocManaged(&present, sizeof(int)*GRID*nbeams);
        cudaMallocManaged(&boxes, sizeof(int)*RAYS*ncrossings*2);
        cudaMemset(boxes, 0,sizeof(int)*RAYS*ncrossings*2);
        cudaMemset(marked, 0,sizeof(int)*GRID*numstored*nbeams);
        cudaMemset(present, 0,sizeof(int)*GRID*nbeams);
        //cudaMallocManaged(&edep, sizeof(double)*RAYS* (nx+2) * (nz+2));
        //cudaMallocManaged(&ints, sizeof(int)*CROSS*numstored);
        cudaMallocManaged(&wpe, sizeof(double)*GRID);
        cudaMallocManaged(&crossesx, sizeof(double)*RAYS*ncrossings);
        cudaMallocManaged(&crossesz, sizeof(double)*RAYS*ncrossings);
    }else
    {
        dedendx = new double[GRID];
        dedendz = new double[GRID];
        x = new double[nx];
        z = new double[nz];
        eden = new double[GRID];
        marked = new int[GRID*numstored*nbeams]{0};
        present = new int[GRID*nbeams]{0};
       // boxes = new int[CROSS*2]{0};
        wpe = new double[GRID];
        crossesx = new double[CROSS]{0.0};
        crossesz = new double[CROSS]{0.0};
    }
    printf("%d\n", nrays);
    double* edenTemp = new double[nx];
    span(edenTemp, 0.1,0.4, nx);
    span(x, xmin, xmax, nx);
    span(z, zmin, zmax, nz);
    for (int i = 0; i < nx; i++)
    {
        double eVal = edenTemp[i]*ncrit;
        for (int j = 0; j < nz; j++)
        {
            vec2DW(present, i, j, nz, 0);
            //eden[i][j] = temp[i];
            vec2DW(eden, i, j, nz, eVal);
            double temp = sqrt(vec2D(eden,i,j,nz)*1e6*pow(ec,2.0)/(me*e0));
            vec2DW(wpe,i,j,nz, temp);
        }
    }
    for(int i = 0; i < nx; i++)
    {
        for(int j = 0; j < nz; j++)
        {
            printf("(%e %d) ", eden[i*nz + j], i*nz+j);
        }
        printf("\n");
    }
    for (int i = 0; i < nx - 1; i++)
    {
        for (int j = 0; j < nz - 1; j++)
        {
            double tempx = (vec2D(eden, i + 1, j, nz) - vec2D(eden, i, j, nz)) / (x[i + 1] - x[i]);
            double tempz = (vec2D(eden, i, j + 1, nz) - vec2D(eden, i, j, nz)) / (z[j + 1] - z[j]);
            vec2DW(dedendx, i, j, nz, tempx);
            vec2DW(dedendz, i, j, nz, tempz);
        }
    }
    for (int i = 0; i < max(nx, nz); i++)
    {
        if (i < nx)
        {
            double temp = vec2D(dedendz, i, nz - 2, nz);
            vec2DW(dedendz, i, nz - 1, nz, temp);
        }
        if (i < nx - 1)
        {
            double tempx = (vec2D(eden, i + 1, nz - 1, nz) - vec2D(eden, i, nz - 1, nz)) / (x[i + 1] - x[i]);
            vec2DW(dedendx, i, nz - 1, nz, tempx);
        }
        if (i < nz - 1)
        {
            double tempz = (vec2D(eden, nx - 1, i + 1, nz) - vec2D(eden, nx - 1, i, nz)) / (z[i + 1] - z[i]);
            vec2DW(dedendz, nx - 1, i, nz, tempz);
        }
        if (i < nz)
        {

            double temp = vec2D(dedendx, nx - 2, i, nz);
            vec2DW(dedendx, nx - 1, i, nz, temp);
        }
    }
}

