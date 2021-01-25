CPP=g++#compilers being used
H5=h5c++

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

H5FLAGS = -g -Wall -Werror -fopenmp -Iinclude#Compiler flags for h5c++
CPPFLAGS= -g  -Wall -std=c++20 -Werror -fopenmp #compiler flags for g++
LIBS = 	-Iinclude -MMD -MP -lm -I/src/include -I/usr/include/python3.8 -lpython3.8 -lcuda#Library Dependecies
HLIBS =  -I/usr/include/hdf5/serial -L/usr/include/hdf5/serial#hdf5 libraries

_MAINOBJ = implSim.o #Main execution files not involved with any individual module
MAINOBJ = $(patsubst %,$(ODIR)/%,$(_MAINOBJ))

_CBETOBJ = cbet.o #CBET Module source files
CBETOBJ = $(patsubst %,$(ODIR)/%,$(_CBETOBJ))

_TRACEOBJ = Launch_Ray_XZ.o RayLaunch.o #Ray tracking module source files
TRACEOBJ = $(patsubst %,$(ODIR)/%,$(_TRACEOBJ))

#_FIELDOBJ =  #Core C++ files being used
#FIELDOBJ = $(patsubst %,$(ODIR)/%,$(_FIELDOBJ))

_INITOBJ = Initialize.o  #Initialization module source files
INITOBJ = $(patsubst %,$(ODIR)/%,$(_INITOBJ))

_IOOBJ =  hdf5writer.o #Core IO module source files
IOOBJ = $(patsubst %,$(ODIR)/%,$(_IOOBJ))

_LIBOBJ =  customMath.o #Core IO module source files
LIBOBJ = $(patsubst %,$(ODIR)/%,$(_LIBOBJ))

$(ODIR)/%.o: $(SRCDIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c -fopenmp -g -o $@ $^ $(LIBS)

$(ODIR)/%.o: $(CBET_DIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c -fopenmp -g -o $@ $^ $(LIBS)
$(ODIR)/%.o: $(FIELD_DIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c -fopenmp -g -o $@ $^ $(LIBS)
$(ODIR)/%.o: $(INIT_DIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c -fopenmp -g -o $@ $^ $(LIBS)
$(ODIR)/%.o: $(TRACE_DIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
	$(CPP) -c -fopenmp -g -o $@ $^ $(LIBS)
$(ODIR)/%.o: $(LDIR)/%.cpp  #$(CBET_DIR)/%.cpp $(FIELD_DIR)/%.cpp $(INIT_DIR)/%.cpp $(TRACE_DIR)/%.cpp#Compile instructions for individual C++ source files
		$(CPP) -c -fopenmp -g -o $@ $^ $(LIBS)
$(ODIR)/%.o:  $(IO_DIR)/%.cpp#Compile instructions for I/O files
	$(H5) -c -fopenmp -g -o $@ $^ $(LIBS)



implSim: $(INITOBJ) $(MAINOBJ) $(LIBOBJ) $(CBETOBJ) $(TRACEOBJ)  $(IOOBJ) #Program compile
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
