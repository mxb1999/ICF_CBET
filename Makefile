CPP=g++#compiler being used
H5=h5c++
MATDIR=Matlab
YORDIR=yorick

IDIR=include
ODIR=Bin
SRC_DIR=src
HDIR=lib/h5
OPDIR=output

H5FLAGS = -g -Wall -Werror -fopenmp -Iinclude
CPPFLAGS= -g  -Wall -std=c++11 -Werror -fopenmp #compiler flags
LIBS = 	-Iinclude -MMD -MP -lm -I/src/include -I/usr/include/python3.8 -lpython3.8#Library Dependecies
HLIBS =  -I/usr/include/hdf5/serial -L/usr/include/hdf5/serial

_COBJ = Initialize.o cbet.o customMath.o hdf5writer.o implSim.o Launch_Ray_XZ.o RayLaunch.o
COBJ = $(patsubst %,$(ODIR)/%,$(_COBJ))

_HOBJ = hdf5writer.o
HOBJ = $(patsubst %,$(ODIR)/%,$(_HOBJ))

$(ODIR)/%.o:  $(SRC_DIR)/%.cpp
	$(CPP) -c -fopenmp -g -o $@ $^ $(LIBS)

$(ODIR)/%.o:  $(HDIR)/%.cpp
	$(H5) -c -fopenmp -g -o $@ $^ $(LIBS)



implSim: $(COBJ) $(HOBJ)
	$(H5) $(CPPFLAGS) -o $@ $^  $(LIBS)

.phony: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
	rm -f $(ODIR)/*.d *~ core $(INCDIR)/*~
	rm -f $(ODIR)/implSim.hdf *~ core $(INCDIR)/*~

compare:
	matlab -nodisplay -nosplash -nodesktop -r "run('m201705_rayTraceTest_v06.m'); exit;" | tail -n +11
	cd yorick && yorick -i cbet.i
	make
