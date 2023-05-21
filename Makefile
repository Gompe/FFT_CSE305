CORE_SRCS = $(wildcard ./core/*.cc)
CORE_H = $(wildcard ./core/*.h)

COMPRESSOR_SRCS = $(wildcard ./compressor/*.cc)
COMPRESSOR_DEPS = $(wildcard ./compressor/*.h)

CXX = g++
CXXFLAGS = -O2 -g -Wall -Wno-psabi -DOPEN_MP

LIBS = -lpthread -fopenmp

INC=-I./

all: $(SRCS) $(DEPS)
	g++ $(CXXFLAGS) $(INC) $(CORE_SRCS) main.cc -lpthread -fopenmp -o main

compressor: $(CORE_SRCS) $(CORE_H) $(COMPRESSOR_SRCS) $(COMPRESSOR_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) $(COMPRESSOR_SRCS) $(CORE_SRCS) $(LIBS) -o compressor.exe