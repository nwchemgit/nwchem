#
# This makefile has to be processed by GNU make
#
include ../config/makefile.h

SHELL = /bin/sh

ifdef OLD_GA
#
# GNUmakefile-based build of GA
# 

export ARMCI_DIR_NAME := armci
ifeq ($(ARMCI_NETWORK),PORTALS)
    export ARMCI_DIR_NAME := armci-portals
endif

ifdef NWCHEM_TARGET_CPU
 TARGET_CPU = $(NWCHEM_TARGET_CPU)
endif
MAKESUBDIRS = for dir in $(SUBDIRS); do $(MAKE)  -C $$dir $@ || exit 1 ; done
    SUBDIRS = ma global tcgmsg-mpi LinAlg/lapack+blas tcgmsg $(ARMCI_DIR_NAME)/src pario
ifndef USE_MPI
  EXTRA_INC = "TCG_INCLUDE=../../include"
      EXTRA = "MSG_COMMS=TCGMSG"
     MP_VER = TCG
else
     MP_VER = MPI
endif
     GA_VER = O

ifndef TARGET
info:
	@echo ''
	@echo '*** This is the Global Arrays package ****'
	@echo "*** you need to define TARGET machine, for example"
	@echo "***     setenv TARGET SUN"
	@echo "*** or "
	@echo "***     make TARGET=SUN"
	@echo ''
endif

export TARGET TARGET_CPU BINDIR

tools:  $(STAMP) lib/$(TARGET) include
	@echo ''
	@echo '*** Building Parallel Tools ****'
	@echo ''
	@echo '>>> making Memory Allocator (MA) library <<< '
	(cd ./ma; $(MAKE))
ifdef USE_MPI
	@echo '>>> making TCGMSG library on top of MPI <<< '
	(cd ./tcgmsg-mpi; $(MAKE) MACHINE=$(TARGET) )
else
	@echo '>>> making TCGMSG library <<<'
	(cd ./tcgmsg; $(MAKE) MACHINE=$(TARGET))
endif
	@echo '>>> making Linear Algebra utility libraries <<< '
	(cd ./LinAlg/lapack+blas; $(MAKE))
	@echo '>>> making ARMCI remote memory copy library <<< '
	(cd ./$(ARMCI_DIR_NAME)/src; \
	$(MAKE) LIB_DISTRIB=../../lib $(EXTRA) $(EXTRA_INC) )
	@echo '>>> making Global Arrays library <<< '
	(cd ./global; $(MAKE) DIAG=PAR )
	(cd ./pario; $(MAKE))

lib/$(TARGET):
	$(MKDIR) -p $@

include: 
	$(MKDIR) -p $@

include_stamp: tools

clean:
	$(MAKESUBDIRS) $^
	$(RM) -r include/* *.stamp
	(cd ./$(ARMCI_DIR_NAME)/src; $(MAKE) LIB_DISTRIB=../../lib clean )

STAMP = $(TARGET)-$(MP_VER)-$(GA_VER).stamp

$(STAMP): 
	$(MAKE) clean
	echo "" > $(STAMP)

realclean: clean

else
#
# GNU Autotools-based build of GA
#

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

ifeq ($(ARMCI_NETWORK),MPI-SPAWN)
    MAYBE_ARMCI = --with-mpi-spawn
endif # MPI-SPAWN

# ARMCI over Multi-threaded MPI
ifeq ($(ARMCI_NETWORK),MPI-MT) 
    MAYBE_ARMCI = --with-mpi-mt
endif # MPI-MT

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
    MAYBE_MPI = --without-mpi
endif # USE_MPI

ifneq ($(CC),)
    MAYBE_CC = CC="$(strip $(CC))"
endif

ifneq ($(FC),)
    MAYBE_F77 = F77="$(strip $(FC))"
else
    ifneq ($(F77),)
        MAYBE_F77 = F77="$(strip $(F77))"
    endif
endif

ifneq ($(BLAS_LIB),)
    ifeq ($(BLAS_SIZE),4)
        MAYBE_BLAS = --with-blas4="$(strip $(BLAS_LIB))"
    else
        MAYBE_BLAS = --with-blas8="$(strip $(BLAS_LIB))"
    endif
else
    MAYBE_BLAS = --without-blas
endif
ifneq ($(BLASOPT),)
     ifeq ($(BLAS_SIZE),4)
         MAYBE_BLAS = --with-blas4="$(strip $(BLASOPT))"
     else
         MAYBE_BLAS = --with-blas8="$(strip $(BLASOPT))"
     endif
endif

ifneq ($(LAPACK_LIB),)
    ifeq ($(LAPACK_SIZE),4)
        MAYBE_LAPACK = --with-lapack="$(strip $(LAPACK_LIB))"
    else
        MAYBE_LAPACK = --with-lapack="$(strip $(LAPACK_LIB))"
    endif
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

ifeq ($(TARGET),BGP)
    MAYBE_HOST = --host=powerpc-bgp-linux 
    FFLAGS="-g"
    CFLAGS="-g"
    CXXFLAGS="-g"
endif

ifeq ($(USE_GPROF),y)
    FFLAGS="-pg"
    CFLAGS="-pg"
    CXXFLAGS="-pg"
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

ifdef DEV_GA
GA_DIR = ga-dev
else
GA_DIR = ga-5-1
endif
BUILDDIR = build
INSTALLDIR = install
CONFIGURE_PATH = ../$(GA_DIR)/configure
TOOLDIR = $(shell pwd)
CONFIGURE_ARGS = --prefix=$(TOOLDIR)/$(INSTALLDIR) $(MAYBE_HOST) --with-tcgmsg $(MAYBE_MPI) --enable-peigs --enable-underscoring --disable-mpi-tests $(MAYBE_SCALAPACK) $(MAYBE_LAPACK) $(MAYBE_BLAS) $(MAYBE_ARMCI) $(MAYBE_CC) $(MAYBE_MPICC) $(MAYBE_F77) $(MAYBE_MPIF77) $(MAYBE_FFLAGS) $(MAYBE_CFLAGS) $(MAYBE_CXXFLAGS)
SPACE :=
SPACE +=
STAMP_FC = $(subst /,-,$(subst $(SPACE),-,$(strip $(FC))).FC.stamp)
STAMP_CC = $(subst /,-,$(subst $(SPACE),-,$(strip $(CC))).CC.stamp)
STAMP_DIR = $(subst $(SPACE),-,$(strip $(GA_DIR))).DIR.stamp

$(INSTALLDIR)/lib/libga.a: $(BUILDDIR)/.libs/libga.a
	@echo ''
	@echo '*** Installing Parallel Tools ****'
	@echo ''
	@(cd $(BUILDDIR) && $(MAKE) install) || exit 1

$(BUILDDIR)/.libs/libga.a: $(BUILDDIR)/Makefile $(BUILDDIR)/armci/Makefile
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
	@echo "You must run ./get-tools" && exit 1
$(GA_DIR)/armci/Makefile.in:
	@echo "You must run ./get-tools" && exit 1
$(GA_DIR)/configure:
	@echo "You must run ./get-tools" && exit 1
$(GA_DIR)/armci/configure:
	@echo "You must run ./get-tools" && exit 1
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

endif
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
