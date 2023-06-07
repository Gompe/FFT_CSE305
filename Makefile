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


polynomial: $(CORE_SRCS) $(CORE_DEPS) $(NT_SRCS) $(NT_DEPS)
	$(CXX) $(CXXFLAGS) $(INC) $(NT_SRCS) $(CORE_SRCS) $(LIBS) main_polynomial.cc -o polynomial.exe

# For the experiments:

TEST_SRCS = $(CORE_SRCS) $(NT_SRCS)
TEST_DEPS = $(CORE_DEPS) $(NT_DEPS) ./tests/benchmark_timer.h
TEST_CMD = $(CXX) $(CXXFLAGS) $(INC) $(TEST_SRCS) $(LIBS) 

test_nt: $(TEST_SRCS) $(TEST) ./tests/test_nt.cc
	$(TEST_CMD) ./tests/test_nt.cc -o test_nt.exe

test_polynomial: $(TEST_SRCS) $(TEST) ./tests/test_polynomial.cc
	$(TEST_CMD) ./tests/test_polynomial.cc -o test_polynomial.exe


test_parallel: $(TEST_SRCS) $(TEST) ./tests/test_parallel.cc
	$(TEST_CMD) ./tests/test_parallel.cc -o test_parallel.exe

test_dft: $(TEST_SRCS) $(TEST) ./tests/test_dft.cc
	$(TEST_CMD) ./tests/test_dft.cc -o test_dft.exe