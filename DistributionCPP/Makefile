CXXFLAGS = -c -g -Wall -O2
OMP = g++
OMPFLAGS = -fopenmp -O2
MPI = mpic++
MPIFLAGS = -O2

all: test

test : PLOT3Ddriver.cc
	$(CXX) PLOT3Ddriver.cc PLOT3D.cc -o test

clean:
	$(RM) *.o
	$(RM) .depend

depend:
	$(CXX) -MM $(CXXFLAGS) *.cc > .depend

-include .depend