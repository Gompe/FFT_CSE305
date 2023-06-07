CORE_SRCS = $(wildcard ./core/*.cc)
CORE_DEPS = $(wildcard ./core/*.h)

COMPRESSOR_SRCS = $(wildcard ./compressor/*.cc)
COMPRESSOR_DEPS = $(wildcard ./compressor/*.h)

NT_SRCS = $(wildcard ./number_theory/*.cc)
NT_DEPS = $(wildcard ./number_theory/*.cc)


CXX = g++
CXXFLAGS = -O2 -g -Wall -Wno-psabi -std=c++17

LIBS = -lpthread -fopenmp

INC=-I./


compressor: $(CORE_SRCS) $(CORE_DEPS) $(COMPRESSOR_SRCS) $(COMPRESSOR_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) $(COMPRESSOR_SRCS) $(CORE_SRCS) $(LIBS) -o compressor.exe


# For the experiments:

TEST_SRCS = $(CORE_SRCS) $(NT_SRCS)
TEST_DEPS = $(CORE_DEPS) $(NT_DEPS)
TEST_CMD = $(CXX) $(CXXFLAGS) $(INC) $(CORE_SRCS) $(LIBS) 

test_nt: $(TEST_SRCS) $(TEST)
