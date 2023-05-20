SRCS = main.cc dft.cc
DEPS = dft.h fft_utils.h fft_types.h my_timer.h

CXXFLAGS = -O2 -g -Wall

all: $(SRCS) $(DEPS)
	g++ $(CXXFLAGS) $(SRCS) -lpthread -fopenmp -o main

omp: $(SRCS) $(DEPS)
	g++ $(CXXFLAGS) -DOPEN_MP $(SRCS) -lpthread -fopenmp -o main