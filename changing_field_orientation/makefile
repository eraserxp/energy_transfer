##  Makefile to make the programtest2.exe#FFLAGS = -O3 -132 -w 
F90 = ifort 
FFLAGSB = -fast -132 -w 
FFLAGSLIB =/global/software/intel/mkl/10.0.4.023/lib/em64t/libmkl_lapack.a -lguide -lpthread

OBJECTSB = wavepacket.f module.f solve_ode2.f zvode_with_deps.f 
#zvode.f zgbfa.f zgbsl.f zgefa.f zgesl.f  

test: $(OBJECTSB)
	$(F90) $(FFLAGSB) -o test  $(OBJECTSB) $(FFLAGSLIB)
# Tab before $(F90)$

clean:
	rm -f *.o test
# Tab before "rm"
