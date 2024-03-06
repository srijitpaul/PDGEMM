CXX = g++
MPICXX = mpic++
CXXFLAGS = -std=c++11
ifdef OPENMP
	CXXFLAGS += -fopenmp
endif
ENVVARS = 

VERBOSE ?= 0

matmul:
	$(CXX) -o $@.bin $@.cpp $(CXXFLAGS) -DVERBOSE=$(VERBOSE)

mpi-matmul:
	$(ENVVARS) $(MPICXX) -o $@.bin $@.cpp $(CXXFLAGS) -DVERBOSE=$(VERBOSE)

all: matmul mpi-matmul
