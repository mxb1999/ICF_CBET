CC=h5c++ #compiler being used
IDIR=include
ODIR=Bin
LDIR=src
vpath %.cpp  src
vpath %.h  src/includes

CFLAGS= -g -I/usr/include/hdf5/serial -L/usr/include/hdf5/serial -Wall -Werror -fopenmp #compiler flags
_DEPS = implSim.h declarations.h simConst.h customMath.h #.h Dependecies
DEBS = $(patsubst %, $(IDIR)/%,$(_DEPS))

LIBS = 	-lm #Library Dependecies
_OBJ = Initialize.o cbet.o customMath.o hdf5writer.o implSim.o Launch_Ray_XZ.o RayLaunch.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o:  $(LDIR)/%.cpp $(DEPS)
	$(CC) -c -fopenmp -g -o $@ $^

implSim: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^  $(LIBS)

.phony: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
