#Simple Makefile
Open Terminal and type commands:
	h5c++ -c -fopenmp -O3  -std=c++17  implSim.cpp -I./hdf5-1.12.0/c++/src -I./hdf5-1.12.0/src
	h5c++ -c -fopenmp -O3  -std=c++17 -I./hdf5-1.12.0/c++/src -I./hdf5-1.12.0/src Initialize.cpp
	h5c++ -c -fopenmp -O3  -std=c++17 -I./hdf5-1.12.0/c++/src -I./hdf5-1.12.0/src cbet.cpp
	h5c++ -c -fopenmp -O3  -std=c++17 -I./hdf5-1.12.0/c++/src -I./hdf5-1.12.0/src RayLaunch.cpp
	h5c++ -c -fopenmp -O3  -std=c++17 -I./hdf5-1.12.0/c++/src -I./hdf5-1.12.0/src Launch_Ray_XZ.cpp
	h5c++ -c -fopenmp -O3 -std=c++17 customMath.cpp
	h5c++ -c -fopenmp -O3  -std=c++17 -I./hdf5-1.12.0/c++/src hdf5writer.cpp
	h5c++ -Wall -fopenmp -O3 -Werror -std=c++17 -I./hdf5-1.12.0/c++/src -I./hdf5-1.12.0/src -g -o implSim implSim.o hdf5writer.o cbet.o  customMath.o RayLaunch.o Launch_Ray_XZ.o Initialize.o
	./implSim
	#python3 matplotting.py
