# Makefile for ersp
FORTRAN = mpif90
#FORTRAN = gfortran
#FLAGS = -O2
#FLAGS =  -fopenmp -g -O3
#FLAGS =  -fopenmp -g -static
FLAGS =  -fopenmp -g
#FLAGS =  -g
SRC1 = rspfn.o mkpol.o testME.o mksecspec.o tetint.o fndgwc.o rdkss.o \
trispline.o spline.o fegrad.o subint.o invertdielf.o \
pspwmatel.o tprojwf.o mkxclda.o mkse.o mkppse.o mksekss.o m_par.o \
par.o mkploss.o mkpwse.o
#BLASDIR = /home/lwcamp/subroutines/Fortran/ATLAS/Build/lib
#LAPACDIR = /home/lwcamp/subroutines/Fortran/ATLAS/Build/lib
#LIBRARIES = -L$(BLASDIR) -lf77blas -latlas -L$(LAPACDIR) -llapack
LIBRARIES = -L/usr/lib64/openmpi/1.4-gcc/lib/
UD = Utilities/
Util = $(UD)linterp.o $(UD)heapsort.o $(UD)grater.o \
$(UD)quadroots.o $(UD)linalg.o $(UD)spherical.o

ersp : $(SRC1) $(Util)
	$(FORTRAN) $(FLAGS) -o ersp $(SRC1) $(Util) $(LIBRARIES)
rspfn.o : rspfn.f90 m_par.o
	$(FORTRAN) $(FLAGS) -c rspfn.f90
mkpol.o : mkpol.f90
	$(FORTRAN) $(FLAGS) -c mkpol.f90
testME.o : testME.f90
	$(FORTRAN) $(FLAGS) -c testME.f90
mksecspec.o : mksecspec.f90
	$(FORTRAN) $(FLAGS) -c mksecspec.f90
mkse.o : mkse.f90
	$(FORTRAN) $(FLAGS) -c mkse.f90
mk3cse.o : mk3cse.f90
	$(FORTRAN) $(FLAGS) -c mk3cse.f90
mkxclda.o : mkxclda.f90
	$(FORTRAN) $(FLAGS) -c mkxclda.f90
mkppse.o : mkppse.f90
	$(FORTRAN) $(FLAGS) -c mkppse.f90
mksekss.o : mksekss.f90
	$(FORTRAN) $(FLAGS) -c mksekss.f90
tetint.o : tetint.f90
	$(FORTRAN) $(FLAGS) -c tetint.f90
fndgwc.o : fndgwc.f90
	$(FORTRAN) $(FLAGS) -c fndgwc.f90
rdkss.o : rdkss.f90
	$(FORTRAN) $(FLAGS) -c rdkss.f90
trispline.o : trispline.f90
	$(FORTRAN) $(FLAGS) -c trispline.f90
spline.o : spline.f90
	$(FORTRAN) $(FLAGS) -c spline.f90
fegrad.o : fegrad.f90
	$(FORTRAN) $(FLAGS) -c fegrad.f90
subint.o : subint.f90
	$(FORTRAN) $(FLAGS) -c subint.f90
invertdielf.o : invertdielf.f90
	$(FORTRAN) $(FLAGS) -c invertdielf.f90
pspwmatel.o : pspwmatel.f90
	$(FORTRAN) $(FLAGS) -c pspwmatel.f90
tprojwf.o : tprojwf.f90
	$(FORTRAN) $(FLAGS) -c tprojwf.f90
par.o : par.f90
	$(FORTRAN) $(FLAGS) -c par.f90
m_par.o : m_par.f90
	$(FORTRAN) $(FLAGS) -c m_par.f90
mkploss.o : mkploss.f90
	$(FORTRAN) $(FLAGS) -c mkploss.f90
mkpwse.o : mkpwse.f90
	$(FORTRAN) $(FLAGS) -c mkpwse.f90
$(UD)linterp.o : $(UD)linterp.f90
	cd $(UD); $(FORTRAN) $(FLAGS) -c linterp.f90
$(UD)heapsort.o : $(UD)heapsort.f90
	cd $(UD); $(FORTRAN) $(FLAGS) -c heapsort.f90
$(UD)grater.o : $(UD)grater.f90
	cd $(UD); $(FORTRAN) $(FLAGS) -c grater.f90
$(UD)quadroots.o : $(UD)quadroots.f90
	cd $(UD); $(FORTRAN) $(FLAGS) -c quadroots.f90
$(UD)linalg.o : $(UD)linalg.f90
	cd $(UD); $(FORTRAN) $(FLAGS) -c linalg.f90
$(UD)spherical.o : $(UD)spherical.f90
	cd $(UD); $(FORTRAN) $(FLAGS) -c spherical.f90
clean	:
	rm -f *.o $(Util) *.mod
