ARFLAGS = -rs
LIB_DDE_DIR ?= .

FC = gfortran
CXX = clang++

FFLAGS = -cpp -O2 -ffree-line-length-none

debug: FFLAGS += -g -DDEBUG -O0
debug: $(LIB_DDE_DIR)/lib/libdde_solver.a

BLAS=-framework Accelerate

all: $(LIB_DDE_DIR)/lib/libdde_solver.a

$(LIB_DDE_DIR)/src/%.o : $(LIB_DDE_DIR)/src/%.f90
	$(FC) $(FFLAGS) $(CPPFLAGS) -J$(dir $@) -c $< -o $@

CXXFLAGS = -I$(LIB_DDE_DIR)/include -std=c++14

$(LIB_DDE_DIR)/lib/libdde_solver.a : $(LIB_DDE_DIR)/src/dde_solver_m.o $(LIB_DDE_DIR)/src/dde_solver_bind.o $(LIB_DDE_DIR)/src/dde_solver_bind_cc.o
	$(AR) $(ARFLAGS) $@ $^

LDFLAGS += $(LIB_DDE_DIR)/libdde_solver.a

clean:
	@find src -name "*.o" -exec rm {} \;
	@find src -name "*.mod" -exec rm {} \;
	@find lib -name "*.a" -exec rm {} \;
