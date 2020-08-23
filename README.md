# Inertial Confinement Fusion Cross Beam Energy Tranfer Simulation (ICF_CBET)
ICF_CBET is a simulation primarily written in C++ with a Python plotting post-processor. The simulation models the phenomenon of "cross beam energy transfer" (CBET) and is based upon the work of
Russ Follett (DOI: 10.1103/PhysRevE.98.043202), and the more basic algorithm used here is derived from a previous simulation written in Yorick by Professor Adam Sefkow (University of Rochester). CBET occurs as a consequence of the ion-acoustic waves created by the beat frequency of overlapping laser beams. These waves mediate an energy transfer between the waves, which has significant consequences for research into inertial confinement fusion, as uniform energy deposition across the target is essential. This simulation is an in-production attempt to model this phenomenon in its simplest case, a 2 beam crossing modeled in 2 dimensions. Note that this was developed in a Linux environment, and therefore the library dependencies are configured for a standard Linux file organization.
## Configuration and Use
The following packages are required to run the simulation and use the post-processor:
* C++11 along with g++
* Python 3 (along with the following packages)
  * Tkinter
  * Matplotlib
  * Numpy
* HDF5
  * h5c++ (found in the hdf5-helpers package)
* OpenMP\n
To install these dependencies automatically (on an Ubuntu based distribution), run:
```bash
sudo make install
```
Alternatively, these packages can be installed separately on any other Linux based distribution.

## Usage
To build the simulation, run ```make``` in the top directory. The object files and dependency files are stored in the Bin folder at compile time. Run the simulation with ```./implSim```. The results of the simulation is stored in ```output/implSim.hdf``` in the HDF5 format. The file I/O is handled by ```lib/h5/hdf5writer.cpp```. All other source files are found in the ```src``` directory, with headers in the ```include``` directory. Various settings for the simulation execution can be found labeled in ```include/simConst.hpp```. These allow to configure for print output, set a maximum number of iterations for the CBET calculation, set the number of threads for OMP vectorization, and determine if the python plotting script should be automatically run upon completion. If the Python script is used, the primary output plots will be shown upon simulation completion. After that window is closed, a separate window will appear that will allow for plotting specific quantities saved in HDF5.
If any of the header files are changed, the command ```make reset``` is necessary to recompile the dependency files. This will delete all .o and .d files in the ```Bin``` directory and rerun ```make```. To delete all .o and .d files without the recompilation, run ```make clean```.

##Limitations:
Note that this software is incomplete, and thus should not be used for any results-dependent work. The field amplitudes at intensities > 1e16 W/cm<sup>2</sup> in the pump beam are far larger than reasonable. Similarly, when the beam radius <= 2e-5 cm, the ray tracer contained in the files ```src/RayLaunch.cpp``` and ```src/Launch_Ray_XZ.cpp``` does not properly record all ray intersections, leaving large portions of beams unchanged even at high intensities.
