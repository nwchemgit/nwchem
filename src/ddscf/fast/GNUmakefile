
 LIBRARY = libddscf.a

  # TESTING := 1

  OPTIM := xlm_poles.o fmm.o cheby.o dposv.o fastj.o anlfit.o cart_trans.o \
           gauleg.o solver.o vsint.o

ifdef TESTING
  OBJ += $(OPTIM)
else
  OBJ_OPTIMIZE := $(OPTIM)
endif

 LIB_TARGETS = solver


include ../../config/makefile.h
include ../../config/makelib.h
    
solver:	testsolver.o $(LIBRARY_PATH)
	$(LINK.f) -o $@ $^ fftpack/libfftpack.a $(LIBS)


