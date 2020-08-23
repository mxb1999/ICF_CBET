CPP=g++#compilers being used
H5=h5c++

#directories
IDIR=include
ODIR=Bin
SRC_DIR=src
HDIR=lib/h5
OPDIR=output

H5FLAGS = -g -Wall -Werror -fopenmp -Iinclude#Compiler flags for h5c++
CPPFLAGS= -g  -Wall -std=c++11 -Werror -fopenmp #compiler flags for g++
LIBS = 	-Iinclude -MMD -MP -lm -I/src/include -I/usr/include/python3.8 -lpython3.8#Library Dependecies
HLIBS =  -I/usr/include/hdf5/serial -L/usr/include/hdf5/serial#hdf5 libraries

_COBJ = Initialize.o cbet.o customMath.o hdf5writer.o implSim.o Launch_Ray_XZ.o RayLaunch.o #Core C++ files being used
COBJ = $(patsubst %,$(ODIR)/%,$(_COBJ))

_HOBJ = hdf5writer.o #HDF5 related files
HOBJ = $(patsubst %,$(ODIR)/%,$(_HOBJ))

$(ODIR)/%.o:  $(SRC_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c -fopenmp -g -o $@ $^ $(LIBS)

$(ODIR)/%.o:  $(HDIR)/%.cpp#Compile instructions for HDF5 I/O files
	$(H5) -c -fopenmp -g -o $@ $^ $(LIBS)



implSim: $(COBJ) $(HOBJ) #Program compile
	$(H5) $(CPPFLAGS) -o $@ $^  $(LIBS)

.phony: clean

#removes all object and dependency files, must be run when a change is made to .h files
clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
	rm -f $(ODIR)/*.d *~ core $(INCDIR)/*~

.phony: reset

reset:
	make clean
	make


.phony: install

install:
	apt-get install libhdf5-serial-dev
	apt-get install hdf5-helpers
	apt-get install python3
	apt-get install python3-matplotlib
	apt-get install python3-numpy
	apt-get install python3-tk
	apt-get install libomp-dev
