Authors:
	Prof. A. Sefkow: Algorithm
	Matthew Burns: C++ Implementation
	University of Rochester

The included simulation models the energy deposition on a fusion target in an intertial confinement fusion (ICF) direct drive scenario. This release is serialized and developed on/for UNIX systems. To compile the program, run "make" in the
folder directory, then ./implSim for any further calculation runs. The python script is used for plotting the values output to HDF5 (implSim.hdf by default) and can be run using "python3 matplotting.py." Note: the -I paths within the makefile
may need to be altered to match the system location of "Python.h" and "H5cpp.h," otherwise a compiletime error will occur. When plotting with the GUI, the first two drop down boxes are for the dimension arrays, and the third is for a 2D array
to be spatially plotted. For instance, to get a plot of intensity distribution over the target before the CBET calculation was updated, plot x, z, and i_b_plot. 
