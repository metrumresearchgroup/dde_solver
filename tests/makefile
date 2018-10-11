FC = gfortran
CXX = clang++

GTEST ?= gtest_1.7.0
GTEST_MAIN = $(GTEST)/src/gtest_main.cc
GTEST_CXXFLAGS = -isystem $(GTEST)/include -isystem $(GTEST) -O$O -DGTEST_USE_OWN_TR1_TUPLE

DDE_SOLVER_DIR := ..

CXXFLAGS = -I$(DDE_SOLVER_DIR)/include -I$(GTEST)/include -std=c++14
LDFLAGS += -L/usr/local/Cellar/gcc/7.3.0_1/lib/gcc/7 -lquadmath -lgfortran
LDFLAGS += -L/usr/local/Cellar//gcc/7.3.0_1/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.3.0 -lgcc
LDFLAGS += -framework accelerate # for macOS to link BLAS/LAPACK
LDFLAGS += -Wl,-no_compact_unwind

gtest : $(GTEST)/src/gtest-all.o
all : sol1_test sol2_test sol3_test sol4_test
debug: CXXLAGS += -g -DDEBUG -O0
debug : sol1_test sol2_test sol3_test sol4_test

$(GTEST)/src/gtest-all.o: CXXFLAGS = $(GTEST_CXXFLAGS) -g

%.o : %.cpp $(DDE_SOLVER_DIR)/include/dde_solver.hpp
	$(COMPILE.cc) -o $@ $<
sol1_test : sol1_test.o $(GTEST_MAIN) $(GTEST)/src/gtest-all.o $(DDE_SOLVER_DIR)/lib/libdde_solver.a
	$(LINK.cpp) -o $@ $^
sol2_test : sol2_test.o $(GTEST_MAIN) $(GTEST)/src/gtest-all.o $(DDE_SOLVER_DIR)/lib/libdde_solver.a
	$(LINK.cpp) -o $@ $^
sol3_test : sol3_test.o $(GTEST_MAIN) $(GTEST)/src/gtest-all.o $(DDE_SOLVER_DIR)/lib/libdde_solver.a
	$(LINK.cpp) -o $@ $^
sol4_test : sol4_test.o $(GTEST_MAIN) $(GTEST)/src/gtest-all.o $(DDE_SOLVER_DIR)/lib/libdde_solver.a
	$(LINK.cpp) -o $@ $^

clean :
	@find . -name "*.o" -exec rm {} \;
	@find . -name "*_test" -exec rm {} \;
	@find . -name "*.dat" -exec rm {} \;
