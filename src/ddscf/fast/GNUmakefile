
 LIBRARY = libddscf.a

  # TESTING := 1

  OBJ = potential.o

  OPTIM := xlm_poles.o fmm.o cheby.o dposv.o fastj.o anlfit.o cart_trans.o \
           gauleg.o solver.o vsint.o dgels.o pot_shell.o  fastj_direct.o \
           fastj_range.o

  USES_BLAS = potential.F cart_trans.f fmm.F pot_shell.F 

ifdef TESTING
  OBJ += $(OPTIM)
else
  OBJ_OPTIMIZE := $(OPTIM)
endif

 LIB_TARGETS = solver testcheby testcarttrans


include ../../config/makefile.h
include ../../config/makelib.h
    
solver:	testsolver.o $(LIBRARY_PATH)
	$(LINK.f) -p -o $@ $^ fftpack/libfftpack.a $(LIBS)

testcheby:	testcheby.o $(LIBRARY_PATH)
	$(LINK.f) -o $@ $^ fftpack/libfftpack.a $(LIBS)

testcarttrans:	testcarttrans.o $(LIBRARY_PATH)
	$(LINK.f) -o $@ $^ fftpack/libfftpack.a $(LIBS)


