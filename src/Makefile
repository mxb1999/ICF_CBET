CC=h5c++ #compiler being used
IDIR=../include
ODIR=obj
LDIR=../src

CFLAGS=-Wall -Werror -O3 -fopenmp #compiler flags
_DEPS = implSim.h declarations.h simConst.h customMath.h #.h Dependecies
DEBS = $(patsubst %, $(IDIR)/%,$(_DEPS))

LIBS = 	-lm #Library Dependecies

_OBJ = Initialize.o cbet.o customMath.o hdf5writer.o implSim.o Launch_Ray_XZ.o RayLaunch.o
OBJ = $(patsubst %, $(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -o $@ $< $(CFLAGS)

implSim: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.phony: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
