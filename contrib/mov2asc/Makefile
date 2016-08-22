# Makefile for movecs <--> ASCII translators
# $Id$

# Tested on KSR, Cray-T3D
# Requires NWCHEM_TARGET to be set.
# Requires NWChem includes and library for MA (in libnwcutil.a)

ifeq ($(NWCHEM_TARGET),CRAY-T3D) # (must run on MPP not host)
   FFLAGS = -Wf"-dp" -Wp"-F" -Ccray-t3d 
endif

ifeq ($(NWCHEM_TARGET),CRAY-T3E)
   FFLAGS = -O 1 -d p -F
   FC     = f90
endif

ifeq ($(NWCHEM_TARGET),LINUX)
   FFLAGS = -fno-second-underscore
endif
ifeq ($(NWCHEM_TARGET),LINUX64)
  ifeq ($(FC),gfortran)
   FFLAGS = -fdefault-integer-8
  else
   FFLAGS = -i8
#   FFLAGS =
  endif	
endif

ifeq ($(NWCHEM_TARGET),KSR)
   FFLAGS = -r8
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

ifdef BLASOPT
  LDFLAGS += $(BLASOPT)
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


