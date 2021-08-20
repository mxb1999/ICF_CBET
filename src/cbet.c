#include "cbet.h"
#include <math.h>


void find_common_crossing_2d()
{

}


void apply_recurrance(Grid* grid, int seed_low, int seed_high, int pump_low, int pump_high, spatial_t* vectors, spatial_t* areas, path_t* trajectories)
{
    int diff = seed_high-seed_low;
    if(diff != 0)
    {
        int offset = ceil((double)diff/2);
        apply_recurrance(grid, seed_low, seed_high-offset, pump_low, pump_high-offset, vectors, areas, trajectories);
        apply_recurrance(grid, seed_low+offset, seed_high, pump_low, pump_high-offset, vectors, areas, trajectories);
        apply_recurrance(grid, seed_low, seed_high-offset, pump_low+offset, pump_high, vectors, areas, trajectories);
        apply_recurrance(grid, seed_low+offset, seed_high, pump_low+offset, pump_high, vectors, areas, trajectories);
    }
    if(seed_low == pump_low)
    {
        return;
    }
    

};


laser_t* calculate_multipliers(Grid* grid, spatial_t* vectors, spatial_t* areas, path_t* trajectories)
{
    int nraycross = NRAYS*CROSSINGS;
    laser_t* multipliers = (laser_t*)malloc(sizeof(laser_t)*nraycross);
    laser_t* old_multipliers = (laser_t*)malloc(sizeof(laser_t)*nraycross);
    int i;
    for(i = 0; i < nraycross; i++)
    {
        multipliers[i] = 1.0;
        old_multipliers[i] = 1.0;
    }


}