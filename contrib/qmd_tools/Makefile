
FC=gfortran
FFLAGS=-O2
OBJ=dsyev.o lsame.o xerbla.o ilaenv.o dlamch.o dlansy.o dlascl.o\
           dsytrd.o dsterf.o dorgtr.o dsteqr.o dscal.o dsyr2k.o\
           iparmq.o ieeeck.o disnan.o dlassq.o dlatrd.o dsytd2.o\
           dlanst.o dlapy2.o dlae2.o dlasrt.o dorgql.o dorgqr.o\
           dlartg.o dlasr.o dlaev2.o dswap.o dlaset.o\
           ddot.o daxpy.o dgemv.o dlarfg.o dsymv.o\
           dsyr2.o dorg2l.o dlarft.o dlarfb.o dorg2r.o\
           dlaisnan.o dnrm2.o dlarf.o dtrmv.o iladlr.o iladlc.o\
           dcopy.o dtrmm.o dgemm.o dger.o

qmd_analysis: qmd_analysis.o $(OBJ)
	$(FC) $(FFLAGS) -o qmd_analysis qmd_analysis.o $(OBJ)

qmd_analysis.o: qmd_analysis.f90
	$(FC) $(FFLAGS) -c qmd_analysis.f90

%.o: $(NWCHEM_TOP)/src/lapack/double/%.f
	$(FC) $(FFLAGS) -c $^

%.o: $(NWCHEM_TOP)/src/blas/double/%.f
	$(FC) $(FFLAGS) -c $^

dgemm.o: $(NWCHEM_TOP)/src/blas/double/dgemm.F
	$(FC) $(FFLAGS) -c $(NWCHEM_TOP)/src/blas/double/dgemm.F

clean:
	rm -f $(OBJ) qmd_analysis.o
