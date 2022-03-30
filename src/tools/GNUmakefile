#
# This makefile has to be processed by GNU make
#
include ../config/makefile.h

SHELL = /bin/sh

#
# GNU Autotools-based build of GA
#
	GET_TOOLS=./get-tools

# Conversion of ARMCI env vars to configure options
ifeq ($(ARMCI_NETWORK),OPENIB)
    OPENIB_COMBINED =
    ifneq ($(IB_INCLUDE),)
        OPENIB_COMBINED += $(IB_INCLUDE)
    endif
    ifneq ($(IB_LIB),)
        OPENIB_COMBINED += $(IB_LIB)
    endif
    ifneq ($(IB_LIB_NAME),)
        OPENIB_COMBINED += $(IB_LIB_NAME)
    endif
    ifeq ($(strip $(OPENIB_COMBINED)),)
        MAYBE_ARMCI = --with-openib
    else
        MAYBE_ARMCI = --with-openib="$(strip $(OPENIB_COMBINED))"
    endif
endif # OPENIB

ifneq ($(findstring BGML,$(ARMCI_NETWORK)),)
    BGML_COMBINED =
    ifneq ($(BGML_INCLUDE),)
        BGML_COMBINED += $(BGML_INCLUDE)
    endif
    ifneq ($(BGMLMPI_INCLUDE),)
        BGML_COMBINED += $(BGMLMPI_INCLUDE)
    endif
    ifneq ($(BGML_LIB),)
        BGML_COMBINED += $(BGML_LIB)
    endif
    ifneq ($(BGMLMPI_LIB),)
        BGML_COMBINED += $(BGMLMPI_LIB)
    endif
    ifneq ($(BGML_LIB_NAME),)
        BGML_COMBINED += $(BGML_LIB_NAME)
    endif
    ifneq ($(BGMLMPI_LIB_NAME),)
        BGML_COMBINED += $(BGMLMPI_LIB_NAME)
    endif
    BGML_COMBINED := $(strip $(BGML_COMBINED))
    ifeq ($(BGML_COMBINED),)
        MAYBE_ARMCI = --with-bgml
    else
        MAYBE_ARMCI = --with-bgml="$(BGML_COMBINED)"
    endif
endif # BGML

ifneq ($(findstring DCMF,$(ARMCI_NETWORK)),)
    DCMF_COMBINED =
    ifneq ($(DCMF_INCLUDE),)
        DCMF_COMBINED += $(DCMF_INCLUDE)
    endif
    ifneq ($(DCMFMPI_INCLUDE),)
        DCMF_COMBINED += $(DCMFMPI_INCLUDE)
    endif
    ifneq ($(DCMF_LIB),)
        DCMF_COMBINED += $(DCMF_LIB)
    endif
    ifneq ($(DCMFMPI_LIB),)
        DCMF_COMBINED += $(DCMFMPI_LIB)
    endif
    ifneq ($(DCMF_LIB_NAME),)
        DCMF_COMBINED += $(DCMF_LIB_NAME)
    endif
    ifneq ($(DCMFMPI_LIB_NAME),)
        DCMF_COMBINED += $(DCMFMPI_LIB_NAME)
    endif
    DCMF_COMBINED := $(strip $(DCMF_COMBINED))
    ifeq ($(DCMF_COMBINED),)
        MAYBE_ARMCI = --with-dcmf
    else
        MAYBE_ARMCI = --with-dcmf="$(DCMF_COMBINED)"
    endif
endif # DCMF

ifeq ($(ARMCI_NETWORK),GEMINI)
    GEMINI_COMBINED =
    ifneq ($(GEMINI_INCLUDE),)
        GEMINI_COMBINED += $(GEMINI_INCLUDE)
    endif
    ifneq ($(GEMINI_LIB),)
        GEMINI_COMBINED += $(GEMINI_LIB)
    endif
    ifneq ($(GEMINI_LIB_NAME),)
        GEMINI_COMBINED += $(GEMINI_LIB_NAME)
    endif
    GEMINI_COMBINED := $(strip $(GEMINI_COMBINED))
    ifeq ($(GEMINI_COMBINED),)
        MAYBE_ARMCI = --with-gemini
    else
        MAYBE_ARMCI = --with-gemini="$(GEMINI_COMBINED)"
    endif
endif # GEMINI

ifeq ($(ARMCI_NETWORK),DMAPP)
    nogacrayerror:
	@echo "  "
	@echo "   we discourage the use of ARMCI_NETWORK=DMAPP "
	@echo "   please use ARMCI_NETWORK=MPI-PR, instead "
	@echo "  " 
	@echo
	@exit 1

    DMAPP_COMBINED =
    ifneq ($(CRAY_DMAPP_INCLUDE_OPTS),)
        DMAPP_COMBINED += $(CRAY_DMAPP_INCLUDE_OPTS)
    endif
    ifneq ($(CRAY_DMAPP_POST_LINK_OPTS),)
        DMAPP_COMBINED += $(CRAY_DMAPP_POST_LINK_OPTS)
    endif
    ifneq ($(DMAPP_INCLUDE),)
        DMAPP_COMBINED += $(DMAPP_INCLUDE)
    endif
    ifneq ($(DMAPP_LIB),)
        DMAPP_COMBINED += $(DMAPP_LIB)
    endif
    ifneq ($(DMAPP_LIB_NAME),)
        DMAPP_COMBINED += $(DMAPP_LIB_NAME)
    endif
    DMAPP_COMBINED := $(strip $(DMAPP_COMBINED))
    ifeq ($(DMAPP_COMBINED),)
        MAYBE_ARMCI = --with-dmapp
    else
        MAYBE_ARMCI = --with-dmapp="$(DMAPP_COMBINED)"
    endif
endif # DMAPP used for Gemini and Aries Cray Networks

ifeq ($(ARMCI_NETWORK),PORTALS)
    PORTALS_COMBINED =
    ifneq ($(PORTALS_INCLUDE),)
        PORTALS_COMBINED += $(PORTALS_INCLUDE)
    endif
    ifneq ($(PORTALS_LIB),)
        PORTALS_COMBINED += $(PORTALS_LIB)
    endif
    ifneq ($(PORTALS_LIB_NAME),)
        PORTALS_COMBINED += $(PORTALS_LIB_NAME)
    endif
    PORTALS_COMBINED := $(strip $(PORTALS_COMBINED))
    ifeq ($(PORTALS_COMBINED),)
        MAYBE_ARMCI = --with-portals
    else
        MAYBE_ARMCI = --with-portals="$(PORTALS_COMBINED)"
    endif
endif # PORTALS

ifeq ($(ARMCI_NETWORK),GM)
    GM_COMBINED =
    ifneq ($(GM_INCLUDE),)
        GM_COMBINED += $(GM_INCLUDE)
    endif
    ifneq ($(GM_LIB),)
        GM_COMBINED += $(GM_LIB)
    endif
    ifneq ($(GM_LIB_NAME),)
        GM_COMBINED += $(GM_LIB_NAME)
    endif
    GM_COMBINED := $(strip $(GM_COMBINED))
    ifeq ($(GM_COMBINED),)
        MAYBE_ARMCI = --with-gm
    else
        MAYBE_ARMCI = --with-gm="$(GM_COMBINED)"
    endif
endif # GM

ifeq ($(ARMCI_NETWORK),VIA)
    VIA_COMBINED =
    ifneq ($(VIA_INCLUDE),)
        VIA_COMBINED += $(VIA_INCLUDE)
    endif
    ifneq ($(VIA_LIB),)
        VIA_COMBINED += $(VIA_LIB)
    endif
    ifneq ($(VIA_LIB_NAME),)
        VIA_COMBINED += $(VIA_LIB_NAME)
    endif
    VIA_COMBINED := $(strip $(VIA_COMBINED))
    ifeq ($(VIA_COMBINED),)
        MAYBE_ARMCI = --with-via
    else
        MAYBE_ARMCI = --with-via="$(VIA_COMBINED)"
    endif
endif # VIA

ifeq ($(ARMCI_NETWORK),MELLANOX)
    MELLANOX_COMBINED =
ifneq ($(IB_INCLUDE),)
        MELLANOX_COMBINED += $(IB_INCLUDE)
endif
ifneq ($(IB_LIB),)
        MELLANOX_COMBINED += $(IB_LIB)
endif
ifneq ($(IB_LIB_NAME),)
        MELLANOX_COMBINED += $(IB_LIB_NAME)
endif
MELLANOX_COMBINED := $(strip $(MELLANOX_COMBINED))
ifeq ($(MELLANOX_COMBINED),)
    MAYBE_ARMCI = --with-openib
else
    MAYBE_ARMCI = --with-openib="$(MELLANOX_COMBINED)"
endif
endif # MELLANOX

ifeq ($(ARMCI_NETWORK),QUADRICS)
endif # QUADRICS
ifeq ($(ARMCI_NETWORK),ELAN3)
endif # ELAN3
ifeq ($(ARMCI_NETWORK),ELAN4)
endif # ELAN4

ifeq ($(ARMCI_NETWORK),LAPI)
    LAPI_COMBINED =
ifneq ($(LAPI_INCLUDE),)
        LAPI_COMBINED += $(LAPI_INCLUDE)
endif
ifneq ($(LAPI_LIB),)
        LAPI_COMBINED += $(LAPI_LIB)
endif
ifneq ($(LAPI_LIB_NAME),)
        LAPI_COMBINED += $(LAPI_LIB_NAME)
endif
LAPI_COMBINED := $(strip $(LAPI_COMBINED))
ifeq ($(LAPI_COMBINED),)
    MAYBE_ARMCI = --with-lapi
else
    MAYBE_ARMCI = --with-lapi="$(LAPI_COMBINED)"
endif
endif # LAPI

# ARMCI over MPI using MPI-Spawn to create the extra data server process
ifeq ($(ARMCI_NETWORK),MPI-SPAWN)
    MAYBE_ARMCI = --with-mpi-spawn
endif # MPI-SPAWN

# ARMCI with progress thread
ifeq ($(ARMCI_NETWORK),MPI-PT) 
    MAYBE_ARMCI = --with-mpi-pt
endif
ifeq ($(ARMCI_NETWORK),MPI_PT) 
    MAYBE_ARMCI = --with-mpi-pt
endif # MPI-PT

# ARMCI over Multi-threaded MPI
ifeq ($(ARMCI_NETWORK),MPI-MT) 
    MAYBE_ARMCI = --with-mpi-mt
endif
ifeq ($(ARMCI_NETWORK),MPI_MT) 
    MAYBE_ARMCI = --with-mpi-mt
endif # MPI-MT

# ARMCI over MPI-1 with progress rank
ifeq ($(ARMCI_NETWORK),MPI-PR) 
    MAYBE_ARMCI = --with-mpi-pr
endif
ifeq ($(ARMCI_NETWORK),MPI_PR) 
    MAYBE_ARMCI = --with-mpi-pr
endif # MPI-PR
ifeq ($(ARMCI_NETWORK),MPI3) 
    MAYBE_ARMCI = --with-mpi3
endif # MPI3

# ARMCI over MPI without data server
ifeq ($(ARMCI_NETWORK),MPI-TS) 
        MAYBE_ARMCI = --with-mpi-ts
endif
ifeq ($(ARMCI_NETWORK),MPI_TS) 
        MAYBE_ARMCI = --with-mpi-ts
endif # MPI_TS

ifeq ($(ARMCI_NETWORK),ARMCI) 
    ifdef EXTERNAL_ARMCI_PATH
        MAYBE_ARMCI = --with-armci=$(EXTERNAL_ARMCI_PATH)
    else
        MAYBE_ARMCI = --with-armci
    endif
endif # ARMCIARMCI

# libfabric port for Omnipath
ifeq ($(ARMCI_NETWORK),OFI) 
    OFI_COMBINED =
    ifneq ($(LIBFABRIC_INCLUDE),)
        OFI_COMBINED += $(LIBFABRIC_INCLUDE)
    endif
    ifneq ($(LIBFABRIC_LIB),)
        OFI_COMBINED += $(LIBFABRIC_LIB)
    endif
    ifneq ($(LIBFABRIC_LIB_NAME),)
        OFI_COMBINED += $(LIBFABRIC_LIB_NAME)
    endif
    ifeq ($(strip $(OFI_COMBINED)),)
        MAYBE_ARMCI = --with-ofi
    else
        MAYBE_ARMCI = --with-ofi="$(strip $(OFI_COMBINED))"
    endif
endif # OFIOFI

ifeq ($(ARMCI_NETWORK),OFA) 
        MAYBE_ARMCI = --with-ofa
endif # OFAOFA

ifeq ($(ARMCI_NETWORK),SOCKETS) 
        MAYBE_ARMCI = --with-sockets
endif # SOCKETS


ifdef USE_MPI

    MPI_COMBINED =
    ifneq ($(MPI_INCLUDE),)
        MPI_COMBINED += $(patsubst -I-I%,-I%,-I${MPI_INCLUDE})
    endif
    ifneq ($(MPI_LIB),)
        MPI_COMBINED += $(patsubst -L-L%,-L%,-L${MPI_LIB})
    endif
    ifneq ($(LIBMPI),)
        MPI_COMBINED += $(LIBMPI)
        ifeq ($(ARMCI_NETWORK),ARMCI)
           MPI_COMBINED += -lpthread
        endif
    endif
    ifneq ($(MPI_LOC),)
        MPI_COMBINED += $(MPI_LOC)
    endif
    MPI_COMBINED := $(strip $(MPI_COMBINED))
    ifeq ($(MPI_COMBINED),)
        MAYBE_MPI = --with-mpi
        ifeq ($(MPICC),)
            ifneq ($(CC),)
                MAYBE_MPICC = MPICC="$(strip $(CC))"
            endif
        else
            MAYBE_MPICC = MPICC="$(strip $(MPICC))"
        endif
        ifeq ($(FC),ftn)
            CXX=CC
        endif
        ifeq ($(MPICXX),)
            ifneq ($(CXX),)
                MAYBE_MPICXX = MPICXX="$(strip $(CXX))"
            endif
        else
            MAYBE_MPICXX = MPICXX="$(strip $(MPICXX))"
        endif
        ifeq ($(MPIF77),)
            ifneq ($(FC),)
                MAYBE_MPIF77 = MPIF77="$(strip $(FC))"
            else
                ifneq ($(F77),)
                    MAYBE_MPIF77 = MPIF77="$(strip $(F77))"
                endif
            endif
        else
            MAYBE_MPIF77 = MPIF77="$(strip $(MPIF77))"
        endif
        ifneq ($(MPIFC),)
            MAYBE_MPIF77 = MPIF77="$(strip $(MPIFC))"
        endif
    else
        MAYBE_MPI = --with-mpi="$(MPI_COMBINED)"
    endif
else # USE_MPI
#old    MAYBE_MPI = --without-mpi
nompierror:
	@echo "  "
	@echo "  " You must set USE_MPI and
	@echo "  " have a working MPI installation
	@echo "  "  to compile NWChem
	@echo
	@exit 1
endif # USE_MPI

ifneq ($(CC),)
    MAYBE_CC = CC="$(strip $(CC))"
endif

ifneq ($(CXX),)
    MAYBE_CXX = CXX="$(strip $(CXX))"
endif

ifneq ($(FC),)
    MAYBE_F77 = F77="$(strip $(FC))"
else
    ifneq ($(F77),)
        MAYBE_F77 = F77="$(strip $(F77))"
    endif
endif

# There really is too much logic here. Overall the aim is that if
# BLASOPT is set to something sensible, i.e. not "" or " " then
# that value should be used, otherwise if BLAS_LIB is set to something
# sensible that value should be used, and otherwise do without an
# external BLAS library and build the internal source code version
# instead. HvD.
ifndef LAPACK_LIB
  ifdef BLAS_LIB
     LAPACK_LIB=$(BLAS_LIB)
  endif
  ifdef BLASOPT
     LAPACK_LIB=$(BLASOPT)
  endif
endif
#enforce 32-bit ints for 32-bit arch
ifeq ($(TARGET),$(findstring $(TARGET),LINUX MACX))
  BLAS_SIZE=4
endif
ifneq ($(BLASOPT),)
    ifneq ($(BLAS_LIB),)
        ifeq ($(BLAS_SIZE),4)
            MAYBE_BLAS = --with-blas4="$(strip $(BLAS_LIB))"
        else
            MAYBE_BLAS = --with-blas8="$(strip $(BLAS_LIB))"
        endif
    else
        ifeq ($(strip $(BLASOPT)),)
            MAYBE_BLAS = --without-blas
        else
            ifeq ($(BLAS_SIZE),4)
                MAYBE_BLAS = --with-blas4="$(strip $(BLASOPT))"
            else
                MAYBE_BLAS = --with-blas8="$(strip $(BLASOPT))"
            endif
        endif
    endif
else
    ifneq ($(BLAS_LIB),)
        ifeq ($(BLAS_SIZE),4)
            MAYBE_BLAS = --with-blas4="$(strip $(BLAS_LIB))"
        else
            MAYBE_BLAS = --with-blas8="$(strip $(BLAS_LIB))"
        endif
    else
        MAYBE_BLAS = --without-blas
    endif
endif

ifneq ($(LAPACK_LIB),)
        MAYBE_LAPACK = --with-lapack="$(strip $(LAPACK_LIB))"
else
    MAYBE_LAPACK = --without-lapack
endif

ifneq ($(SCALAPACK_LIB),)
    ifeq ($(SCALAPACK_SIZE),4)
        MAYBE_SCALAPACK = --with-scalapack="$(strip $(SCALAPACK_LIB))"
    else
        MAYBE_SCALAPACK = --with-scalapack8="$(strip $(SCALAPACK_LIB))"
    endif
endif
ifneq ($(ELPA),)
    ifeq ($(ELPA_SIZE),8)
        MAYBE_ELPA = --with-elpa8="$(strip $(ELPA))"
    else
        MAYBE_ELPA = --with-elpa="$(strip $(ELPA))"
    endif
endif

ifneq ($(SCALAPACK),)
    ifeq ($(SCALAPACK_SIZE),8)
        MAYBE_SCALAPACK = --with-scalapack8="$(strip $(SCALAPACK))"
    else
        MAYBE_SCALAPACK = --with-scalapack="$(strip $(SCALAPACK))"
    endif
endif

ifeq ($(SCALAPACK_LIB),)
    ifeq ($(SCALAPACK),)
        MAYBE_SCALAPACK = --without-scalapack
    endif
endif

ifeq ($(TARGET),$(findstring $(TARGET),MACX MACX64))
#gtimeout from brew coreutils is required
#  GTIMEOUT_YES=$(shell rm -f gtimeout.tmp ;/usr/local/bin/gtimeout  --version > gtimeout.tmp 2>&1; head -1 gtimeout.tmp | awk ' /GNU/ {print "Y";exit};{print "N"}')
  GTIMEOUT_YES=$(shell /usr/local/bin/zgtimeout  --version 2>&1 | head -1  | awk ' /GNU/ {print "Y";exit};{print "N"}')
  ifdef USE_GAGITHUB
  ifeq ($(GTIMEOUT_YES),N)
    gtimeout_fail:
	@echo " "
	@echo "   gtimeout missing"
	@echo "   Please install coreutils from Homebrew "
	@echo " "
	exit 1
  endif
  endif
endif

ifeq ($(TARGET),BGQ)
    FFLAGS_FORGA+=-qsmp
endif

ifeq ($(TARGET),BGP)
    MAYBE_HOST = --host=powerpc-bgp-linux 
    FFLAGS="-g"
    CFLAGS="-g"
    CXXFLAGS="-g"
endif

ifdef USE_GPROF
    FFLAGS="-pg"
    CFLAGS="-pg"
    CXXFLAGS="-pg"
    FFLAGS_FORGA+="-pg"
    CFLAGS_FORGA+="-pg"
    CXXFLAGS_FORGA+="-pg"
endif

ifdef USE_KNL
    MAYBE_OFFLOAD= INTEL_64ALIGN=1
    ifeq ($(_FC),ifort)
    FFLAGS_FORGA+="-no-vec"
    endif
    ifeq ($(_CC),icc)
    CFLAGS_FORGA+="-no-vec"
    endif
    ifeq ($(FC),ftn)
    MAYBE_HOST=--host=x86_64-pc-linux-gnu 
    endif
endif
ifdef USE_OFFLOAD
    MAYBE_OFFLOAD= INTEL_64ALIGN=1
endif
ifdef USE_EISPACK
    MAYBE_EISPACK= --enable-eispack
endif
    MAYBE_SYSVSHMEM= "ARMCI_DEFAULT_SHMMAX_UBOUND=131072"
GOTCGWIN32= $(shell $(CC) -dM -E - </dev/null 2> /dev/null |grep CYGWIN32|cut -c22)
GOTMINGW32= $(shell $(CC) -dM -E - </dev/null 2> /dev/null |grep MINGW32|cut -c21)
GOTMINGW64= $(shell gcc -dM -E -   </dev/null 2> /dev/null |grep MINGW64|cut -c21)
ifeq ($(GOTCGWIN32),1)
# weak pragmas not working with cygwin
   MAYBE_OFFLOAD += --disable-weak
endif
ifeq ($(GOTMINGW32),1)
# weak pragmas not working with cygwin
#MSMPI madness
   FFLAGS_FORGA += -fno-range-check 
   MAYBE_OFFLOAD += --disable-weak
endif
ifeq ($(GOTMINGW64),1)
#MSMPI madness
   FFLAGS_FORGA += -fno-range-check -m64
   CFLAGS_FORGA +=  -m64
   CC = gcc
endif
ifeq ($(GNU_GE_8),true)
  FFLAGS_FORGA += -std=legacy
endif
#detect clang
GOTCLANG= $(shell $(CC) -dM -E - </dev/null 2> /dev/null |grep __clang__|head -1|cut -c19)
GOTFREEBSD= $(shell uname -o 2>&1|awk ' /FreeBSD/ {print "1";exit}')
ifeq ($(_FC),ifort)
  USE_FPICF=1
endif
ifeq ($(_FC),ifx)
  USE_FPICF=1
endif
ifeq ($(_FC),pgf90)
  USE_FPICF=1
endif
ifdef USE_FPICF
  FFLAGS_FORGA += "-fPIC"
endif
ifdef USE_CRAYSHASTA
  CFLAGS_FORGA += -D__CRAYXE
endif
ifeq ($(GOTCLANG),1)
  CFLAGS_FORGA += "-fPIC"
  ifeq ($(TARGET),MACX64)
    CFLAGS_FORGA +=-DMPICH_NO_ATTR_TYPE_TAGS
  endif
endif
ifeq ($(PE_ENV),CRAY)
   CFLAGS_FORGA +=-O0 -g -L/opt/cray/pe/mpt/default/gni/mpich-cray/90/lib
   FFLAGS_FORGA +=-O0 -g
endif

ifeq ($(GOTFREEBSD),1)
    UNSETFLAGS_FORGA = LDFLAGS=" " CPPFLAGS=" " LIBS=" "
endif

ifdef FFLAGS_FORGA
    MAYBE_FFLAGS=FFLAGS="$(FFLAGS_FORGA)"
endif
ifdef CFLAGS_FORGA
    MAYBE_CFLAGS=CFLAGS="$(CFLAGS_FORGA)"
endif
ifdef CXXFLAGS_FORGA
    MAYBE_CXXFLAGS=CXXFLAGS="$(CXXFLAGS_FORGA)"
endif
# End conversion of env vars to configure options
USE_GAGITHUB=1
ifdef USE_GASVN
 USE_GAGITHUB=
endif
USE_GA_RELEASE=1
ifdef GA_BRANCH
 USE_GA_RELEASE=
endif
ifdef DEV_GA
 USE_GA_RELEASE=
endif

ifdef USE_GAGITHUB
  MAKEFLAGS = -j 1 --no-print-directory
	GET_TOOLS=./get-tools-github
  ifdef USE_GA_RELEASE
    ifndef GA_RELEASE_NO
      GA_RELEASE_NO=5.8.1
    endif
    GA_DIR=ga-$(GA_RELEASE_NO)
  else
#DEV_GA=1 gives the current "best" dev branch, not develop
#to get develop, set GA_BRANCH-develop
    ifdef DEV_GA
      GA_DIR = ga-develop
      GA_BRANCH = develop
    else
      GA_DIR0=ga-$(GA_BRANCH)
      GA_DIR=$(shell echo $(GA_DIR0) | sed -e 's/\//_/g')
    endif
  endif
else
ifdef DEV_GA
    GA_DIR = ga-dev
else
    ifdef EXP_GA
        GA_DIR = ga-exp1
    else
        ifdef GA_STABLE
           GA_DIR = $(GA_STABLE)
        else
           GA_DIR = ga-5-5
        endif
    endif
endif
endif
ifndef ARMCI_NETWORK
    ARMCI_NETWORK=MPI-TS
    MAYBE_ARMCI = --with-mpi-ts
endif
# CUDA UM support
ifdef NWCHEM_LINK_CUDA
    MAYBE_ARMCI +=  --enable-cuda-mem
endif
#
# Apparently weak bindings do not work with CYGWIN64 at the moment. There seems
# to be an issue with the COFF object format that gets in the way (with ELF
# this is apparently not a problem). Some details can be found at
# http://cygwin.com/ml/cygwin-developers/2013-04/msg00008.html
# Anyway, I found when using weak bindings the code fails to link as the weak
# symbols cannot be found. HVD (2014 Nov 4).
#
ifeq ($(TARGET),CYGWIN64) 
        MAYBE_ARMCI += --disable-weak
endif # CYGWIN64
#SCORE-P
ifeq ($(CC),$(findstring $(CC),scorep scorep-wrapper scorep-gcc scorep-mpicc scorep-cc scorep-nvcc scorep-xlc scorep-xlc_r))
#fix for configure or cmake
     MAYBE_HOST += SCOREP_WRAPPER=OFF
endif
ifeq ($(FC),$(findstring $(FC),scorep scorep-wrapper scorep-gfortran scorep-mpif90 scorep-ftn scorep-xlf scorep-xlf_r scorep-mpif77))
#fix for configure or cmake
     MAYBE_HOST += SCOREP_WRAPPER=OFF
endif
BUILDDIR = build
INSTALLDIR = install
CONFIGURE_PATH = ../$(GA_DIR)/configure
TOOLDIR = $(shell pwd)
CONFIGURE_ARGS = --prefix=$(TOOLDIR)/$(INSTALLDIR) $(MAYBE_HOST) --with-tcgmsg $(MAYBE_MPI) --enable-peigs --enable-underscoring --disable-mpi-tests $(MAYBE_SCALAPACK) $(MAYBE_ELPA) $(MAYBE_LAPACK) $(MAYBE_BLAS) $(MAYBE_ARMCI) $(MAYBE_CC) $(MAYBE_MPICC) $(MAYBE_CXX) $(MAYBE_MPICXX) $(MAYBE_F77) $(MAYBE_MPIF77) $(MAYBE_FFLAGS) $(MAYBE_CFLAGS) $(MAYBE_CXXFLAGS) $(MAYBE_OFFLOAD) $(MAYBE_SYSVSHMEM) $(MAYBE_EISPACK) $(UNSETFLAGS_FORGA)
SPACE :=
SPACE +=
STAMP_FC = $(subst /,-,$(subst $(SPACE),-,$(strip $(FC))).FC.stamp)
STAMP_CC = $(subst /,-,$(subst $(SPACE),-,$(strip $(CC))).CC.stamp)
ifdef ARMCI_NETWORK
STAMP_DIR = $(subst $(SPACE),-,$(strip $(GA_DIR))).$(subst $(SPACE),-,$(strip $(ARMCI_NETWORK))).DIR.stamp
else
STAMP_DIR = $(subst $(SPACE),-,$(strip $(GA_DIR))).DIR.stamp
endif

$(INSTALLDIR)/lib/libga.a: $(BUILDDIR)/.libs/libga.a
	@echo ''
	@echo '*** Installing Parallel Tools ****'
	@echo ''
	@(cd $(BUILDDIR) && $(MAKE) install) || exit 1
	@(../util/util_ga_version.bash)
	${GA_PATH}/bin/ga-config --cppflags  > $(NWCHEM_TOP)/src/ga_cppflags.txt
	${GA_PATH}/bin/ga-config --ldflags  > $(NWCHEM_TOP)/src/ga_ldflags.txt
	${GA_PATH}/bin/ga-config --use_scalapack | awk ' /1/ {print "Y";exit};{print "N"}'  > $(NWCHEM_TOP)/src/ga_use_scalapack.txt

$(BUILDDIR)/.libs/libga.a:  $(BUILDDIR)/Makefile $(BUILDDIR)/armci/Makefile
	@echo ''
	@echo '*** Building Parallel Tools ****'
	@echo ''
	@(cd $(BUILDDIR) && $(MAKE)) || exit 1

configure_ga: $(BUILDDIR)/Makefile

$(BUILDDIR)/Makefile: $(GA_DIR)/Makefile.in $(BUILDDIR)/config.status
	cd $(BUILDDIR) && $(SHELL) ./config.status

$(BUILDDIR)/armci/Makefile: $(GA_DIR)/armci/Makefile.in $(BUILDDIR)/armci/config.status
	cd $(BUILDDIR)/armci && $(SHELL) ./config.status

$(BUILDDIR)/armci/config.status: $(BUILDDIR)/config.status
$(BUILDDIR)/config.status: $(GA_DIR)/configure $(STAMP_FC) $(STAMP_CC) $(STAMP_DIR)
	@echo ''
	@echo '*** Configuring Parallel Tools ****'
	@echo ''
	@(test -d $(BUILDDIR)) || mkdir $(BUILDDIR);
	@(cd $(BUILDDIR) && echo $(CONFIGURE_PATH) $(CONFIGURE_ARGS) && $(CONFIGURE_PATH) $(CONFIGURE_ARGS)) || exit 1
$(GA_DIR)/Makefile.in:
	@echo "No Global Arrays found: Attempting running ./get-tools"
	$(GET_TOOLS)
$(GA_DIR)/configure: $(GA_DIR)/Makefile.in
	@echo "No Global Arrays found: Attempting running ./get-tools"
	$(GET_TOOLS)
	touch $(GA_DIR)/configure
ifndef DONTCARE
$(GA_DIR)/armci/Makefile.in: $(GA_DIR)/configure
	@echo "No Global Arrays found: Attempting running ./get-tools"
	$(GET_TOOLS)
	touch $(GA_DIR)/armci/Makefile.in
$(GA_DIR)/armci/configure: $(GA_DIR)/armci/Makefile.in
	@echo "3No Global Arrays found: Attempting running ./get-tools"
	$(GET_TOOLS)
	touch $(GA_DIR)/armci/configure
endif # DONTCARE
$(STAMP_FC): 
	rm -f ./*.FC.stamp
	echo "" > $(STAMP_FC)
$(STAMP_CC): 
	rm -f ./*.CC.stamp
	echo "" > $(STAMP_CC)
$(STAMP_DIR): 
	rm -rf $(BUILDDIR) $(INSTALLDIR)
	rm -f ./*.DIR.stamp
	echo "" > $(STAMP_DIR)

clean:
	@test -e $(BUILDDIR)/Makefile && { cd $(BUILDDIR) && $(MAKE) clean; } || echo "Not configured"

realclean:
	@rm -rf $(BUILDDIR) $(INSTALLDIR)
# remove all files that are not .o objects (eg lib*.a) from libga.a
cleanlibga:
	mkdir -p ob
	cd ob && ar -x ../build/.libs/libga.a && rm -f ../build/.libs/libga.a && ar crv ../build/.libs/libga.a *.o && mkdir -p ../install/lib && cp ../build/.libs/libga.a ../install/lib
	cd ..
	rm -rf ob


#
# Common targets shared between GNUmakefile's and automake Makefiles
#

.PHONY: cleanF
cleanF:
	@echo Not done here
.PHONY: cleanDEP
cleanDEP:
	@echo Not done here
source:

depend:
	@echo No dependency analysis here
sngl_to_dbl dbl_to_sngl 64_to_32 32_to_64:
	@echo No blas conversion done here
#
# This is a convenience target that will make the TAGS file for current 
# checked out source tree.  This is only useful if you know something 
# about emacs.  Note: find, grep and etags must be in your path.
#
tags_file:
	find . \( -name "*.[cfFh]" -o -name "*.fh" \) -print | grep -v "\./include" | grep -v "\./tools/include" | etags -
tags_clean:
	find . -name TAGS -print -exec rm -f "{}" ";"
