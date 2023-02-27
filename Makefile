CPP=g++-8#compilers being usedi
NV=nvcc

#directories
IDIR=include
ODIR=Bin
SRCDIR=src

HDIR=lib/h5
LDIR=lib
OPDIR=output
CU_ODIR=Bin/cuda

LIBS = -lpython3.8   -L/usr/local/cuda/lib64/ -lcudadevrt -lcudart -L/usr/local/cuda/lib64/ -L/usr/lib64/hdf5 -lhdf5 -lhdf5_cpp #Library Dependecies
INT_INCLUDE = -Iinclude
EXT_INCLUDE = -I/usr/include/python3.8 -I/usr/include/cuda -I/usr/include/hdf5/serial

REFS = $(INT_INCLUDE) $(EXT_INCLUDE)
H5FLAGS = -g -Wall  -fopenmp -fPIC#Compiler flags for h5c++
CPPFLAGS= -g  -Wall -MMD -MP  -fopenmp -lm -fPIC #compiler flags for g++
NVFLAGS =  -std=c++11 -G -g -Xcompiler -fopenmp -rdc=true -Xcompiler -fPIC


_MAINOBJ = implSim.o #Main execution files not involved with any individual module
MAINOBJ = $(patsubst %,$(ODIR)/%,$(_MAINOBJ))

_CBETOBJ = fillCbetArrays.o map_fields.o cbet2.o# CBET Module source files
CBETOBJ = $(patsubst %,$(ODIR)/%,$(_CBETOBJ))

_DEVMANOBJ = gpuInterface.o #Device management Module source files
DEVMANOBJ = $(patsubst %,$(ODIR)/%,$(_CBETOBJ))

_TRACEOBJ = fillTraceArrays.o Launch_Ray_XZ.o RayLaunch.o #Ray tracking module source files
TRACEOBJ = $(patsubst %,$(ODIR)/%,$(_TRACEOBJ))

#_FIELDOBJ =  #Core C++ files being used
#FIELDOBJ = $(patsubst %,$(ODIR)/%,$(_FIELDOBJ))

_IOOBJ =  hdf5writer.o importConst.o#Core C++ IO module source files
IOOBJ = $(patsubst %,$(ODIR)/%,$(_IOOBJ))

_LIBOBJ =  customMath.o #Core IO module source files
LIBOBJ = $(patsubst %,$(ODIR)/%,$(_LIBOBJ))

_CUOBJ = trackray.o cudahelper.o cbetMemOps.o traceMemOps.o
CUOBJ = $(patsubst %,$(CU_ODIR)/%,$(_CUOBJ))

$(ODIR)/%.o: $(SRCDIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c $(CPPFLAGS) $(REFS) -o $@ $^ $(LIBS)
$(CU_ODIR)/%.o: $(SRCDIR)/%.cu#$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(NV) -c $(NVFLAGS) $(REFS) $^ -o $@ $(LIBS)

implSim: $(MAINOBJ) $(LIBOBJ) $(CBETOBJ) $(TRACEOBJ) $(DEVMANOBJ) $(IOOBJ) $(CUOBJ)#Program compile
	$(NV) -dlink $(NVFLAGS) $(CUOBJ) $(REFS) -o $(CU_ODIR)/cumulativeGPU.o $(LIBS)
	$(CPP)  $(CPPFLAGS)   $(REFS) $^ $(CU_ODIR)/cumulativeGPU.o -o $@ $(LIBS)

.phony: clean

#removes all object and dependency files, must be run when a change is made to .h files
clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
	rm -f $(ODIR)/*.d *~ core $(INCDIR)/*~
	rm -f $(CU_ODIR)/*.o *~ core $(INCDIR)/*~
.phony: reset

reset:
	make clean
	make


.phony: deb install

deb install:
	apt-get install libhdf5-serial-dev
	apt-get install hdf5-helpers
	apt-get install python3
	apt-get install python3-matplotlib
	apt-get install python3-numpy
	apt-get install python3-tk
	apt-get install libomp-dev

.phony:run

run:
	make
	./implSim

.phony:plot

plot:
	python3 matplotting.py
