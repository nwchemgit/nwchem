# Makefile for movecs <--> ASCII translators
# $Id$

# Tested on KSR, Cray-T3D
# Requires NWCHEM_TARGET to be set.
# Requires NWChem includes and library for MA (in libnwcutil.a)
ifndef NWCHEM_TOP
NWCHEM_TOP= $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST)))| \
 sed -e 's/\/contrib.*//' )
endif
ifndef NWCHEM_TARGET
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S),Linux)
    NWCHEM_TARGET=LINUX64
  else ifeq ($(UNAME_S),Darwin)
    NWCHEM_TARGET=MACX64
  else
error2:
$(info     )
$(info You must define NWCHEM_TARGET in your environment to be the name)
$(info of the machine you wish to build for ... for example)
$(info     setenv NWCHEM_TARGET SOLARIS)
$(info Known targets are SOLARIS, ...)
$(info See the INSTALL instructions for a complete list)
$(error )
  endif
endif

ifeq ($(FC),f77)
  FC=gfortran
endif
ifeq ($(NWCHEM_TARGET),LINUX)
   FFLAGS = -fno-second-underscore
endif
ifeq ($(NWCHEM_TARGET),$(findstring $(NWCHEM_TARGET), LINUX64 MACX64))
  ifeq ($(FC),gfortran)
   FFLAGS = -fdefault-integer-8
  else
   FFLAGS = -i8
#   FFLAGS =
  endif	
endif


ifeq ($(NWCHEM_TARGET),IBM)
   FFLAGS = -qEXTNAME
endif

ifeq ($(NWCHEM_TARGET),IBM)
   CPPFLAGS = -WF,"-I../../src/tools/install/include"		# Needs MA includes
else
   CPPFLAGS = -I../../src/tools/install/include			# Needs MA includes
endif
LDFLAGS = -L../../src/tools/install/lib -lga -larmci	# Needs MA library
ifdef USE_MPI 
#ifeq ($(FC),$(findstring $(FC),mpifrt mpfort mpif77 mpxlf mpif90 ftn))
ifeq ($(FC),$(findstring $(FC), ftn))
  LIBMPI =
  MPI_INCLUDE =
  MPI_LIB =
else
ifndef MPI_INCLUDE
# check if mpif90 is present
         MPIF90YN = $(shell $(NWCHEM_TOP)/src/tools/guess-mpidefs --mpi_include)
  ifeq ($(MPIF90YN),mpif90notfound)
  errormpif90:
	@echo " "
	@echo "mpif90 not found. Please add its location to PATH"
	@echo "e.g. export PATH=/usr/local/bin:/usr/lib64/openmpi/bin:..."
	@echo " "
  endif
        MPI_INCLUDE = $(shell $(NWCHEM_TOP)/src/tools/guess-mpidefs --mpi_include)
endif
ifndef MPI_LIB
        MPI_LIB     = $(shell $(NWCHEM_TOP)/src/tools/guess-mpidefs --mpi_lib)
endif
ifndef LIBMPI
        LIBMPI      = $(shell $(NWCHEM_TOP)/src/tools/guess-mpidefs --libmpi)
endif
endif
ifdef MPI_LIB 
      LDFLAGS += $(patsubst -L-L%,-L%,-L$(MPI_LIB))
endif 
ifdef OLD_GA
  LDFLAGS += -ltcgmsg-mpi $(LIBMPI) 
else
  LDFLAGS += $(LIBMPI)
endif
endif
ifneq ($(USE_GPROF),)
  FFLAGS += -pg
  LDFLAGS += -pg
endif
COMM_LIBS=  $(shell grep ARMCI_NETWORK_LIBS\ = ${NWCHEM_TOP}/src/tools/build/Makefile | cut -b 22-)
LDFLAGS += $(COMM_LIBS)

#ifdef BLASOPT
#  LDFLAGS += $(BLASOPT)
#endif
#ifdef BLAS_LIB
#  LDFLAGS += $(BLAS_LIB)
#endif
#ifdef BUILD_OPENBLAS
#  LDFLAGS += -L$(NWCHEM_TOP)/src/libext/lib -lnwc_openblas -lpthread
#endif
ifneq ("$(wildcard ${NWCHEM_TOP}/src/tools/install/bin/ga-config)","")
  LDFLAGS+=  $(shell ${NWCHEM_TOP}/src/tools/install/bin/ga-config --ldflags  )
  LDFLAGS+=  $(shell ${NWCHEM_TOP}/src/tools/install/bin/ga-config --libs  )
endif

all:    mov2asc asc2mov

mov2asc:	mov2asc.o
	${FC} ${FFLAGS} ${CPPFLAGS} -o $@ $^ ${LDFLAGS} 

asc2mov:	asc2mov.o
	${FC} ${FFLAGS} ${CPPFLAGS} -o $@ $^ ${LDFLAGS} 


clean realclean:
	$(RM) *.o mov2asc asc2mov



%.o:	%.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) $<

%.o:	%.f
	$(FC) -c $(FFLAGS) $<

%.o:	%.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<


