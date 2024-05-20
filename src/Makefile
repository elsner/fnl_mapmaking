
FC      = mpif90
FC     += -O0 -fpp -qopenmp -traceback -warn all -check all
#FC     += -O3 -xHost -fpp -qopenmp

MODDIR  = ${HOME}/healpix/healpix_mod

INCDIR  = $(patsubst %,-I%,$(subst :, ,${INCLUDE}))
LIBDIR  = $(patsubst %,-L%,$(subst :, ,${LIBRARY_PATH}))

FFLAGS  = $(INCDIR)
LDFLAGS = $(LIBDIR)\
          -lmkl_lapack95_lp64 -lmkl_blas95_lp64\
          -Wl,--start-group\
           ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a\
           ${MKLROOT}/lib/intel64/libmkl_intel_thread.a\
           ${MKLROOT}/lib/intel64/libmkl_core.a\
          -Wl,--end-group\
          -lpthread\
          -lhealpix -lcfitsio

exec    = run
src     = map_main.f90
obj     = map_mod_classes.o map_mod_misc.o\
	  map_mod_sample.o map_main.o


default: $(src) $(exec)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

$(exec): $(obj)
	$(FC) -o $(exec) $(obj) $(LDFLAGS)

clean:
	rm -f *~ *.o *.mod ${MODDIR}/*.o batch_*

cleanall: clean
	rm -f $(exec)
