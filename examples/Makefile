FC = gfortran
DDE_SOLVER_DIR = ..
PROGS = b2 b2g c2g ex4p4p1 ex4p4p2 ex4p4p3 ex4p4p4 ex4p4p5 ex4p7 ex4p8 exsd1 \
        infection laserex neutvan secdelay volinteq

DDE_LIB = $(DDE_SOLVER_DIR)/lib/libdde_solver.a

all: $(PROGS)

%.o: %.f90
	$(FC) -c -I$(DDE_SOLVER_DIR)/src $(FFLAGS) $< -o $@

b2: b2.o
	$(FC) $< $(DDE_LIB) -o $@
b2g: b2g.o
	$(FC) $< $(DDE_LIB) -o $@
c2g: c2g.o
	$(FC) $< $(DDE_LIB) -o $@
ex4p4p1: ex4p4p1.o
	$(FC) $< $(DDE_LIB) -o $@
ex4p4p2: ex4p4p2.o
	$(FC) $< $(DDE_LIB) -o $@
ex4p4p3: ex4p4p3.o
	$(FC) $< $(DDE_LIB) -o $@
ex4p4p4: ex4p4p4.o
	$(FC) $< $(DDE_LIB) -o $@
ex4p4p5: ex4p4p5.o
	$(FC) $< $(DDE_LIB) -o $@
ex4p7: ex4p7.o
	$(FC) $< $(DDE_LIB) -o $@
ex4p8: ex4p8.o
	$(FC) $< $(DDE_LIB) -o $@
exsd1: exsd1.o
	$(FC) $< $(DDE_LIB) -o $@
infection: infection.o
	$(FC) $< $(DDE_LIB) -o $@
laserex: laserex.o
	$(FC) $< $(DDE_LIB) -o $@
neutvan: neutvan.o
	$(FC) $< $(DDE_LIB) -o $@
secdelay: secdelay.o
	$(FC) $< $(DDE_LIB) -o $@
q1damodule.o: q1damodule.f90
	$(FC) -c $(FFLAGS) $< -o $@
volinteq: q1damodule.o volinteq.o 
	$(FC) $^ $(DDE_LIB) -o $@

clean:
	rm -f $(PROGS) *.o *.mod
	rm -f *.dat
