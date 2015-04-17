CXXFLAGS = -std=c++0x -c -g -O2
OMP = g++
OMPFLAGS = -fopenmp -O2
MPI = mpic++
MPIFLAGS = -O2
VPATH = ./QuadTreeCustom/

all: icing

icing : IcingDriver.o PLOT3D.o Bucket.o
	$(CXX) $^ -o $@

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend
