
 LIBRARY = libddscf.a

#   TESTING := 1

  OBJ = 

  OPTIM := xlm_poles.o cheby.o dposv.o fastj.o anlfit.o cart_trans.o \
           gauleg.o solver.o vsint.o dgels.o pot_shell.o  fastj_direct.o \
           fastj_range.o newfmm.o potential.o

LIB_DEFINES = -DFMM_LMAX=30

ifdef TESTING
  OBJ += $(OPTIM)
else
  OBJ_OPTIMIZE := $(OPTIM)
endif

 LIB_TARGETS = solver testcheby testcarttrans testfmm


include ../../config/makefile.h
include ../../config/makelib.h

#FOPTIONS += -pg

# Vast did not speed up pot_shell or xlm_poles
#FOPTIONS += -Pv -Wp,-eaj78,-o
# We don't have KAP
#FOPTIONS += -Pk -Wp,-r=3,-s=3,-o=5,-lc=essl,-l,-f
#FOPTIMIZE := -O3 -qfloat=rsqrt:fltint -NQ40000 -NT80000 
#-lmass
    
solver:	testsolver.o $(LIBRARY_PATH)
	$(LINK.f) -p -o $@ $^  $(LIBS)

testcheby:	testcheby.o $(LIBRARY_PATH)
	$(LINK.f) -o $@ $^ $(LIBS)

testcarttrans:	testcarttrans.o $(LIBRARY_PATH)
	$(LINK.f) -o $@ $^ $(LIBS)

testfmm:	testfmm.o $(LIBRARY_PATH)
	$(LINK.f) -p -o $@ $^  $(LIBS)

newtestfmm:	newtestfmm.o $(LIBRARY_PATH)
	$(LINK.f) -p -o $@ $^  $(LIBS)


