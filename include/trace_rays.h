#include "cbet_sim.h"
#ifndef TRACE

    #define TRACE
    typedef struct
    {
        path_t* trajectories;
        spatial_t* wave_vectors;
        spatial_t* areas;
    }TraceResult;

    typedef struct
    {
        spatial_t init_location[3];
        spatial_t init_vector[3];
        path_t* trajectories;
        spatial_t* wave_vectors;
        spatial_t* areas;
        int ray, nt;
        double dt;
    }TraceInput;


    TraceResult* new_TraceResult(path_t* trajectories, laser_t* vectors, laser_t* areas);
    TraceResult* launch_rays(Grid* grid);
    void trace_ray(Grid* grid, TraceInput input);
#endif
