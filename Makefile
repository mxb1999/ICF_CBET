CPP=g++#compilers being used
H5=h5c++
NV=nvcc
#directories
IDIR=include
ODIR=Bin
SRCDIR=src
CBET_DIR=src/CBET_Mod
FIELD_DIR=src/FieldSolve_Mod
INIT_DIR=src/Init_Mod
IO_DIR=src/IO_Mod
TRACE_DIR=src/Trace_Mod
HDIR=lib/h5
LDIR=lib
OPDIR=output
CU_ODIR=Bin/cuda


H5FLAGS = -g -Wall -Werror -fopenmp -Iinclude#Compiler flags for h5c++
CPPFLAGS= -g  -Wall -MMD -MP -Werror -fopenmp -lm #compiler flags for g++
LIBS = 	-Iinclude -Iinclude/GPU -I/src/include -I/usr/include/python3.8 -lpython3.8 -I/usr/include/cuda  -L/usr/local/cuda/lib64 -lcudart -lcuda -I/usr/include/hdf5/ -L/usr/lib/hdf5 -lhdf5 #Library Dependecies
HLIBS =  -I/usr/include/hdf5/serial -L/usr/include/hdf5/serial#hdf5 libraries
NVFLAGS =  -std=c++11 -Xcompiler -fopenmp -Xcompiler -fPIC
_MAINOBJ = implSim.o #Main execution files not involved with any individual module
MAINOBJ = $(patsubst %,$(ODIR)/%,$(_MAINOBJ))

_CBETOBJ = cbet.o #CBET Module source files
CBETOBJ = $(patsubst %,$(ODIR)/%,$(_CBETOBJ))

_TRACEOBJ = Launch_Ray_XZ.o RayLaunch.o#Ray tracking module source files
TRACEOBJ = $(patsubst %,$(ODIR)/%,$(_TRACEOBJ))

#_FIELDOBJ =  #Core C++ files being used
#FIELDOBJ = $(patsubst %,$(ODIR)/%,$(_FIELDOBJ))

_INITOBJ = Initialize.o  #Initialization module source files
INITOBJ = $(patsubst %,$(ODIR)/%,$(_INITOBJ))

_IOOBJ =  hdf5writer.o gpuInterface.o cudahelper.o#Core IO module source files
IOOBJ = $(patsubst %,$(ODIR)/%,$(_IOOBJ))

_LIBOBJ =  customMath.o #Core IO module source files
LIBOBJ = $(patsubst %,$(ODIR)/%,$(_LIBOBJ))

_CUOBJ = trackray.o
CUOBJ = $(patsubst %,$(CU_ODIR)/%,$(_CUOBJ))

$(ODIR)/%.o: $(SRCDIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c -fPIC -fopenmp -g -o $@ $^ $(LIBS)
$(ODIR)/%.o: $(CBET_DIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c -fPIC -fopenmp -g -o $@ $^ $(LIBS)
$(ODIR)/%.o: $(FIELD_DIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c -fPIC -fopenmp -g -o $@ $^ $(LIBS)
$(ODIR)/%.o: $(INIT_DIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c -fPIC -fopenmp -g -o $@ $^ $(LIBS)
$(ODIR)/%.o: $(TRACE_DIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c -fPIC -fopenmp -g -o $@ $^ $(LIBS)
$(ODIR)/%.o: $(LDIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
		$(CPP) -c -fPIC -fopenmp -g -o $@ $^ $(LIBS)
$(ODIR)/%.o:  $(IO_DIR)/%.cpp#Compile instructions for I/O files
	$(H5) -c -fPIC -fopenmp -g -o $@ $^ $(LIBS)
$(CU_ODIR)/%.o: $(TRACE_DIR)/%.cu  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(NV) -c $(NVFLAGS) -rdc=true -g $^ -o $@  $(LIBS)


implSim: $(INITOBJ) $(MAINOBJ) $(LIBOBJ) $(CBETOBJ) $(TRACEOBJ)  $(IOOBJ) #Program compile
	$(NV) -dlink $(NVFLAGS) $(CU_OBJ) -o Bin/gpu.o $(LIBS)
	$(H5)  $(CPPFLAGS)   $^ Bin/gpu.o -o $@ $(LIBS)

.phony: clean

#removes all object and dependency files, must be run when a change is made to .h files
clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
	rm -f $(ODIR)/*.d *~ core $(INCDIR)/*~

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
