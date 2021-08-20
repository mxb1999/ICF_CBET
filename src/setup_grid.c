#include "cbet_sim.h"
#define GRIDRES 50
Grid* setupGrid2D()
{
    Grid* grid = (Grid*)malloc(sizeof(Grid));
    IFNOMEM(grid);
    NX = GRIDRES;
    NY = GRIDRES;
    XMIN = -5e-4;
    XMAX = 5e-4;
    YMIN = -5e-4;
    YMAX = 5e-4;
    NDIMS = 2;
    NBEAMS = 2;
    RAYS_PER_ZONE = 3;
    grid->eden = (spatial_t*)malloc(sizeof(spatial_t)*NX*NY);
    grid->machnum = (spatial_t*)malloc(sizeof(spatial_t)*NX*NY);
    grid->wpe = (spatial_t*)malloc(sizeof(spatial_t)*NX*NY);
    grid->beam_centers = (spatial_t*)malloc(sizeof(spatial_t)*NBEAMS*NDIMS);
};