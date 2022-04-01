
# $Id$
#

# Common definitions for all makefiles ... these can be overridden
# either in each makefile by putting additional definitions below the
# include statement, or on the command line

#
# TOPDIR points to your top-level directory that contains 
# src, lib, config, ... (SRCDIR, etc., are derived from TOPDIR)
# Do a setenv for NWCHEM_TOP to be the top level directory
#
# RELEASE is empty for the development branch of NWChem and is the
# version number for releases.  If RELEASE is not empty and if NWCHEM_TOP 
# does not already contain the value of RELEASE, it is appended with 
# a hyphen to NWCHEM_TOP in order to derive the directory path TOPDIR.

# For development tree 
#RELEASE := 
# For current release tree
RELEASE := 7.0.0

#
ifndef NWCHEM_TOP
    #error1:
	#$(info     )
	#$(info You must define NWCHEM_TOP in your environment to be the path)
	#$(info of the top level nwchem directory ... something like)
	#$(info     setenv NWCHEM_TOP /msrc/home/elvis/nwchem)
	#$(info     )
	#$(error )
    NWCHEM_TOP= $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST)))| \
    sed -e 's/\/src.*//' )
endif

# Select the old (pre-autotools version of GA) by uncommenting the next line.
# The value of OLD_GA does not matter -- it is detected as an ifdef only.
#OLD_GA = y 

#
# Do a setenv for NWCHEM_TARGET to be the machine and NWCHEM_TARGET_CPU the CPU to build for
#
# NWCHEM_TARGET :  
#                  CYGNUS       (Windows under Cygwin tools)
#                  IBM
#                  LINUX        NWCHEM_TARGET_CPU :
#                                                  nothing for X86 (e.g. do not set this)
#                                                  POWERPC for MkLinux 
#                  SOLARIS      NWCHEM_TARGET_CPU : not defined or ULTRA
#                  LAPI         NWCHEM_TARGET_CPU : P2SC
#                      (uses thread safe libraries and LAPI)
#
#

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

ifneq ($(NWCHEM_TARGET),$(findstring $(NWCHEM_TARGET), LINUX64 LINUX MACX MACX64 CATAMOUNT CYGWIN64 CYGNUYS CYGWIN BGL BGP BGQ HPUX HPUX64 IBM IBM64 LAPI LAPI64 PURESOLARIS SOLARIS SOLARIS64))
    error20:
	$(info     )
	$(info unrecognized NWCHEM_TARGET value $(NWCHEM_TARGET))
	$(info     )
	$(error )
endif


# trying to fix the BLAS_LIB BLASOPT madness ...
# they have extactly the same purpose (why did we  introduce BLAS_LIB?)
# BLASOPT has higher priority
ifdef BLASOPT
    BLAS_LIB = $(BLASOPT)
else
    ifdef BLAS_LIB
        BLASOPT = $(BLAS_LIB)
    endif
endif


TARGET := $(NWCHEM_TARGET)
TOPDIR := $(NWCHEM_TOP)


ifeq (,$(RELEASE))
    CODE_BRANCH := Development
else
    CODE_BRANCH := $(RELEASE)
endif

#dummy:
#   @echo NWCHEM_TOP=$(NWCHEM_TOP) RELEASE=$(RELEASE) TOPDIR=$(TOPDIR)


SRCDIR := $(TOPDIR)/src
ifndef NWCHEM_TARGET_CPU
    LIBDIR := $(TOPDIR)/lib/$(NWCHEM_TARGET)
    BINDIR := $(TOPDIR)/bin/$(NWCHEM_TARGET)
else
    LIBDIR := $(TOPDIR)/lib/$(NWCHEM_TARGET)_$(NWCHEM_TARGET_CPU)
    BINDIR := $(TOPDIR)/bin/$(NWCHEM_TARGET)_$(NWCHEM_TARGET_CPU)
endif

INCDIR := $(TOPDIR)/src/include
CNFDIR := $(TOPDIR)/src/config


ifdef EXTERNAL_GA_PATH
   #check if ga-config is there
    ifeq ("$(wildcard ${EXTERNAL_GA_PATH}/bin/ga-config)","")
        $(info  )
        $(info invalid EXTERNAL_GA_PATH)
        $(info ga-config not found)
        $(info  )
        $(error )
    endif

    #check if f77 was enabled
    GA_HAS_F77 = $(shell ${EXTERNAL_GA_PATH}/bin/ga-config --enable-f77 | awk '/yes/ {print "Y"}')
    ifneq ($(GA_HAS_F77),Y)
        $(info NWChem requires Global Arrays built with Fortran support)
        $(error )
    endif

    #check peigs interface       
    GA_HAS_PEIGS = $(shell ${EXTERNAL_GA_PATH}/bin/ga-config --use_peigs | awk '/1/ {print "Y"}')
    GA_HAS_SCALAPACK = $(shell ${EXTERNAL_GA_PATH}/bin/ga-config --use_scalapack | awk '/1/ {print "Y"}')
    ifneq ($(GA_HAS_PEIGS),Y)
        ifneq ($(GA_HAS_SCALAPACK),Y)
            $(info NWChem requires Global Arrays built with either Peigs or Scalapack support)
            $(error )
        endif
    endif

    #check blas size
    GA_BLAS_SIZE = $(shell ${EXTERNAL_GA_PATH}/bin/ga-config --blas_size)
    ifndef BLAS_SIZE
        BLAS_SIZE=8
    endif

    ifneq ($(BLAS_SIZE),$(GA_BLAS_SIZE))
         $(info )
         $(info NWChem requires Global Arrays built with same BLAS size )
         $(info you asked BLAS_SIZE=${BLAS_SIZE})
         $(info Global Arrays was built with BLAS size=${GA_BLAS_SIZE})
         $(info )
         $(error )
    endif

    GA_PATH=$(EXTERNAL_GA_PATH)

else
    GA_PATH=$(NWCHEM_TOP)/src/tools/install
endif


ifeq ($(BLAS_SIZE),4)
    USE_64TO32=y
endif
      
#
# Define LIBPATH to be paths for libraries that you are linking in
# from precompiled sources and are not building now. These libraries
# will be searched AFTER anything you are building now.
# e.g. LIBPATH = -L/msrc/proj/mss/lib
#
LIBPATH = 

ifdef OLD_GA
    LIBPATH = -L$(SRCDIR)/tools/lib/$(TARGET)
else
   #case guard against case when tools have not been compiled yet
    ifeq ("$(wildcard ${GA_PATH}/bin/ga-config)","")
        LIBPATH = -L$(SRCDIR)/tools/install/lib
    else

        ifneq ("$(wildcard ${NWCHEM_TOP}/src/ga_ldflags.txt)","")
            GA_LDFLAGS= $(shell cat $(NWCHEM_TOP)/src/ga_ldflags.txt)
        endif

        ifeq ($(GA_LDFLAGS),)
            GA_LDFLAGS=  $(shell ${GA_PATH}/bin/ga-config --ldflags  )
        endif

      #extract GA libs location from last word in GA_LDLFLAGS
        LIBPATH :=  $(word $(words ${GA_LDFLAGS}),${GA_LDFLAGS}) 
        ifdef EXTERNAL_GA_PATH
            LIBPATH += -L$(shell $(NWCHEM_TOP)/src/tools/guess-mpidefs --mpi_lib)
        endif
    endif
endif


#
# Define INCPATH to be directories to get includes for
# libraries that you are not building now.  These directories
# will be searched AFTER anything you are building now.
#
INCPATH = 
ifdef OLD_GA
    INCPATH = -I$(SRCDIR)/tools/include
else
    #case guard against case when tools have not been compiled yet
    ifeq ("$(wildcard ${GA_PATH}/bin/ga-config)","")
        INCPATH = -I$(SRCDIR)/tools/install/include
    else
        ifneq ("$(wildcard ${NWCHEM_TOP}/src/ga_cppflags.txt)","")
            GA_CPPFLAGS= $(shell cat $(NWCHEM_TOP)/src/ga_cppflags.txt)
        endif

        ifeq ($(GA_CPPFLAGS),)
            GA_CPPFLAGS=  $(shell ${GA_PATH}/bin/ga-config --cppflags  )
        endif

        INCPATH :=  $(word $(words ${GA_CPPFLAGS}),${GA_CPPFLAGS})

        ifdef EXTERNAL_GA_PATH
            INCPATH += -I$(shell $(NWCHEM_TOP)/src/tools/guess-mpidefs --mpi_include)
        endif
    endif
endif

# These subdirectories will build the core, or supporting libraries
# that are required by all NWChem modules.  The include directory is
# first to insure that all include files will be properly installed
# prior to trying to compile anything.
#
# The core libraries are usually rather platform-dependent and are
# specified below.  Use of MPI requires substituting the tcgmsg-mpi
# wrapper for the normal tcgmsg library.
# the 2 following environmental variables are need for linking
# LIBMPI - represents the name of mpi library (with -l)
# MPI_LIB - represents the path to the mpi library
#LIBMPI =  -lmpich
#MPI_LIB= /usr/local/lib

#JN: under the new structure, tools should be listed first as
# their header files are needed for dependency analysis of
# other NWChem modules
ifdef USE_LIBXC
    NW_CORE_SUBDIRS += libext
endif


ifdef BUILD_OPENBLAS
    ifndef BLAS_SIZE
        BLAS_SIZE=8
    endif

    NW_CORE_SUBDIRS += libext

    #bail out if BLASOPT or LAPACK_LIB or BLAS_LIB are defined by user
    ifneq ($(or $(BLASOPT),$(LAPACK_LIB),$(BLAS_LIB)),)
        $(info     )
        $(info You must unset)
        $(info BLASOPT ,LAPACK_LIB and BLAS_LIB)
        $(info when using BUILD_OPENBLAS )
        $(info )
        $(error )
    endif
    BLASOPT=-L$(NWCHEM_TOP)/src/libext/lib -lnwc_openblas -lpthread
    LAPACK_LIB=$(BLASOPT)      
    BLAS_LIB=$(BLASOPT)
endif


ifdef BUILD_SCALAPACK
    NW_CORE_SUBDIRS += libext
    ifneq ($(or $(SCALAPACK),$(SCALAPACK_LIB)),)
        $(info     )
        $(info You must unset)
        $(info SCALAPACK  and SCALAPACK_LIB)
        $(info when using BUILD_SCALAPACK )
        $(info )
        $(error )
    endif

    ifndef SCALAPACK_SIZE
        SCALAPACK_SIZE=8
    endif
    SCALAPACK=-L$(NWCHEM_TOP)/src/libext/lib -lnwc_scalapack
endif


ifdef BUILD_ELPA
    NW_CORE_SUBDIRS += libext
#   ifeq ($(or $(BUILD_SCALAPACK),$(BUILD_OPENBLAS)),)
#       $(info     )
#       $(info You must set)
#       $(info BUILD_SCALAPACK  and BUILD_OPENBLAS)
#       $(info when using BUILD_ELPA )
#       $(info )
#       $(error )
#   endif

    ifndef SCALAPACK_SIZE
        SCALAPACK_SIZE=8
    endif
    ELPA=-L$(NWCHEM_TOP)/src/libext/lib -lnwc_elpa -I$(NWCHEM_TOP)/src/libext/include/elpa/modules
endif


ifdef BUILD_MPICH
    NW_CORE_SUBDIRS += libext
    PATH := $(NWCHEM_TOP)/src/libext/bin:$(PATH)
    MPI_INCLUDE = $(shell PATH=$(NWCHEM_TOP)/src/libext/bin:$(PATH) $(NWCHEM_TOP)/src/tools/guess-mpidefs --mpi_include)
    MPI_LIB     = $(shell PATH=$(NWCHEM_TOP)/src/libext/bin:$(PATH)  $(NWCHEM_TOP)/src/tools/guess-mpidefs --mpi_lib)
    LIBMPI      = $(shell PATH=$(NWCHEM_TOP)/src/libext/bin:$(PATH) $(NWCHEM_TOP)/src/tools/guess-mpidefs --libmpi)
endif


ifndef EXTERNAL_GA_PATH
    NW_CORE_SUBDIRS += tools
endif

NW_CORE_SUBDIRS += include basis geom inp input  \
                  pstat rtdb task symmetry util peigs perfm bq cons $(CORE_SUBDIRS_EXTRA)



# Include the modules to build defined by 'make nwchem_config' at top level
include $(CNFDIR)/nwchem_config.h
#

# Finally, we can set the full list of interesting directories, which
# is what most makefile will care about.

NWSUBDIRS = $(NW_CORE_SUBDIRS) $(NW_MODULE_SUBDIRS)

BUILDING_PYTHON = $(filter $(NWSUBDIRS),python)


##########################################################
#                                                        #
# Should NOT need to modify below here unless porting to #
# a new machine or changing compiler options             #
#                                                        #
##########################################################

# Each machine dependent section should define the following as necessary.
# (defaults if any in parentheses)
#
#         FC = path to Fortran compiler (f77)
#         CC = path to ANSI C compiler (cc)
#         AS = path to the assembler (as)
#         AR = path to archive builder (ar)
#        CPP = path to ANSI-like C preprocessor (cpp)
#     RANLIB = path to ranlib or something harmless if not required
#      SHELL = path to the Bourne Shell
#       MAKE = DON\'T define this ... it will break the passing of command
#              arguments.  Simply use the correct path to GNU make on the 
#              command line and all will work just dandy.
#  MAKEFLAGS = options to GNU make ... -j controls no. of threads used
#              for parallel builds. --no-print-directory says be quiet about
#              changing directory.
#    INSTALL = command to install an executable when it is built
#    
# C/FOPTIONS = essential compiler options independent of optimization level
#              C/FOPTIONS should not usually be overridden on the command line
#   C/FDEBUG = compiler flags to enable enable debugging and used for
#              all routines not vital for performance (OBJ in makelib.h)
#C/FOPTIMIZE = compiler flags to enable optimization for important routines.
#              (OBJ_OPTIMIZE in makelib.h)
#
#              C/FDEBUG and C/FOPTIMIZE can be overridden on the command 
#              line to change the optimization level for routines normally 
#              compiled with them.
#
#  EXPLICITF = undefined if the Fortran compiler runs .F files thru .f
#              Otherwise set it to anything and define FCONVERT to be a 
#              command to make $< (which will be a .F file) into $*.f
#
#   FCONVERT = command to convert a .F into a .f 
#
#  LDOPTIONS = additional options to be passed to the linker (LDFLAGS is
#              built from this and the library path info).  LDOPTIONS is
#              the best way to add to the link command.
#
#    ARFLAGS = options for AR (ru)
#
#    DEFINES = C preprocessor defines for both C and Fortran
#
#  CORE_LIBS = List of libraries and paths for libraries in addition 
#              to the LIBDIR and LIBPATH options.
#
# CORE_SUBDIRS_EXTRA = List of additional directories (e.g., BLAS) that
#                 are needed on this machine.
# MODULE_SUBDIRS_EXTRA = List of additional directories (e.g., stepper) that
#                        are needed for top-level modules on this machine.
#                        (Should not normally be used)
#
# The following are defined for all machines at the bottom of this file
#
#   C/FFLAGS = all options to the C/Fortran compilers (note CPPFLAGS are 
#              separate).  These comprise C/FOPTIONS and C/FOPT.
#   INCLUDES = C preprocessor include paths for both C and Fortran.
#              In princpile this could be machine dependent but is not yet.
#   CPPFLAGS = options to C preprocessor that both C and Fortran will use.
#              This comprises the includes and defines.
#    LDFLAGS = options for the linker.  Currently paths from LIBDIR and
#              LIBPATH.
#
# NWCHEM_TARGET_CPU environment variable is used to select compiler flags
# optimal for given CPU. The following values are recognized on particular
# platforms:
#
#        NWCHEM_TARGET                NWCHEM_TARGET_CPU          
#           LAPI                         P2SC
#


#
# Establish some required defaults which may need overriding
# for some machines

SHELL = /bin/sh
ARFLAGS = r
FDEBUG = -g
CDEBUG = -g
AR = ar
FCONVERT = $(CPP) $(CPPFLAGS) $< > $*.f


ifdef OLD_GA
    CORE_LIBS = -lnwcutil -lpario -lglobal -lma -lpeigs -lperfm -lcons -lbq -lnwcutil
else
    CORE_LIBS = -lnwcutil -lga -larmci -lpeigs -lperfm -lcons -lbq -lnwcutil
endif


ifdef USE_INTEGER4
    integer4:
	@echo 
	@echo USE_INTEGER4 option no longer supported
	@echo please do not set it
	@echo 
	@exit 1
endif

# -U option needed for introducing timestamps in libraries
# see https://bugzilla.redhat.com/show_bug.cgi?id=1195883  (RedHat backtracked)
# see https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=798913                                                                     
# see https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=798804 (Debian did not backtrack)                                          
#      USE_ARUR = $(shell rm -f aru.tmp;ar -U > aru.tmp 2>&1; head -1 aru.tmp| awk ' /no operation/ {print "Y";exit};{print "N"}'; rm -f aru.tmp)

USE_ARUR = $(shell rm -f aru.tmp;ar --help  > aru.tmp 2>&1; grep U aru.tmp| awk ' /ctual timest/ {print "Y";exit};'; rm -f aru.tmp)

ifeq ($(USE_ARUR),Y)
    ARFLAGS = rU 
endif


# strip long paths when FC and/or CC are set from user
ifneq ($(FC),f77)
    _FC = $(notdir $(FC))
endif

ifneq ($(CC),cc)
    _CC = $(notdir $(CC))
endif

#
# Machine specific stuff
#
ifeq ($(TARGET),SOLARIS)
    SHELL := $(NICE) /bin/sh
    RANLIB = echo
    MAKEFLAGS = -j 2 --no-print-directory
    INSTALL = echo $@ is built
    #
    # You can use either the f77 or f90 compiler BUT if using f90
    # you\'ll need to specify -DINTEGER_1='integer*1' in the selci
    # and util makefiles.
    #
    CC = cc
    FC = f77

    # Don\'t need this if using the SUN performance library
    #  need for BLASOPT business because of lapack and other possible missing entries
    ifndef BLASOPT
        CORE_SUBDIRS_EXTRA = blas lapack 
    endif

    DEFINES = -DSOLARIS  -DNOAIO
    # Note that WS6 does not optimize robustly and if using this you must 
    #   - put "-nodpend -xvector=no" on FOPTIONS after -fast
    #   - remove "-fsimple=2 -depend -xvector=yes" from FOPTIMIZE.
    #   - remove -lmvec from CORELIBS
    # to get link with sunperf, type BLASOPT="-xlic_lib=sunperf"
    #
    # These options are set for WS5

    ifeq ($(CC),fcc)
        # Fujitsu SPARC systems (thanks to Herbert Fruchtl)
        COPTIONS = -Kdalign
        COPTIMIZE = -Kfast_GP=2
        DEFINES += -DFUJITSU_SOLARIS
    endif

    ifeq ($(FC),frt)
        #    Fujitsu SPARC systems (thanks to Herbert Fruchtl)
        # Fujitsu with Parallelnavi compilers
        # If using Fujitsu compilers on Sun hardware, replace -Kfast_GP=2 with
        #  -Kfast
        DEFINES += -DFUJITSU_SOLARIS -DEXTNAME
        FOPTIONS = -Kdalign -w -fw -X9 
        FOPTIMIZE = -Kfast_GP=2 
        FDEBUG=
        LINK.f = $(FC) $(LDFLAGS) $(FOPTIONS) $(FOPTIMIZE)
    else
        FOPTIONS = -stackvar -dalign 
        FOPTIMIZE = -fast -O5 -fsimple=2 -depend -xvector=yes
        FDEBUG = -g -O1 -nodepend
        LINK.f = $(FC) $(LDFLAGS) $(FOPTIONS)
    endif
 
    ifeq ($(FC),frt)
        CORE_LIBS +=  -SSL2
    else
        LDOPTIONS = -xildoff
        CORE_LIBS +=  -lnwclapack $(BLASOPT) -lnwcblas  -lmvec
    endif
 
 
    EXTRA_LIBS = -ldl 
    # this creates a static executable
    #EXTRA_LIBS = -Bdynamic -ldl -lXext -lnsl  -Bstatic  
 
 
    ifeq ($(BUILDING_PYTHON),python)
    # Both tk/tcl and BLT need X11 (common)
        EXTRA_LIBS += -lX11
    endif
 
    #end of solaris
endif



ifeq ($(TARGET),SOLARIS64)
#
# Sun running Solaris 64-bit ... NEEDS WORKSHOP 6.1 or later compilers
# due to bugs in earlier compilers.  Also cannot yet use sunperf due to
# the braindead naming of the 64-bit interface.
#
# You can use either the f77 or f90 compiler BUT if using f90
# you will need to specify -DINTEGER_1='integer*1' in the selci
# and util makefiles.
#

    SHELL := $(NICE) /bin/sh
    CORE_SUBDIRS_EXTRA = blas lapack
    CC = cc
    FC = f77
    DEFINES = -DSOLARIS  -DNOAIO -DSOLARIS64
    DEFINES  +=  -DEXT_INT

    COPTIMIZE = -O
    RANLIB = echo
    MAKEFLAGS = -j 2 --no-print-directory
    INSTALL = echo $@ is built

    ifeq ($(CC),fcc)
        # Fujitsu SPARC systems (thanks to Herbert Fruchtl)
        COPTIONS = -Kdalign -KV9FMADD
        COPTIMIZE = -Kfast_GP=2 -KV9FMADD
        DEFINES += -DFUJITSU_SOLARIS
    else
        # SUN/Solaris options for WS6.1
        COPTIONS = -xarch=v9 -dalign
    endif

    ifeq ($(FC),frt)
        # Fujitsu SPARC systems (thanks to Herbert Fruchtl)
        # Fujitsu with Parallelnavi compilers
        # If using Fujitsu compilers on Sun hardware, replace -Kfast_GP=2 with
        #  -Kfast
        DEFINES += -DFUJITSU_SOLARIS -DEXTNAME
        FOPTIONS = -Kdalign -w -fw -X9  -KV9FMADD
        ifdef USE_I4FLAGS
            FOPTIONS += -CcdLL8
        else
            FOPTIONS += -CcdLL8 -CcdII8
        endif
        FOPTIMIZE = -Kfast_GP=2  -KV9FMADD
        FDEBUG=
    else
        # SUN/Solaris f77 options 
        FOPTIONS = -stackvar -fast -nodepend -xvector=no -xarch=v9a
        ifdef USE_I4FLAGS
            FOPTIONS +=  -xtypemap=real:64,double:64,integer:32
        else
            FOPTIONS +=  -xtypemap=real:64,double:64,integer:64
        endif
        FOPTIMIZE = -g -O5
        FDEBUG = -g -O1
    endif

    LINK.f = $(FC) $(LDFLAGS) $(FOPTIONS)
    ifeq ($(FC),frt)
        LDOPTIONS = -SSL2
        CORE_LIBS +=  -lnwclapack -lnwcblas
    else
        LDOPTIONS = -xs -xildoff
        ifdef BLASOPT
            CORE_LIBS +=   $(BLASOPT)   -lmvec
        else
            CORE_LIBS +=  -lnwclapack -lnwcblas  -lmvec
        endif
        CORE_LIBS += -lsocket -lrpcsvc -lnsl
        EXTRA_LIBS =  -ldl -lfsu
    endif
    #end of solaris
endif


ifeq ($(TARGET),PURESOLARIS)
    @echo DEPRECATED
    @exit
endif


ifeq ($(TARGET),cray-sv2)
    @echo DEPRECATED
    @exit
endif


ifeq ($(TARGET),CRAY-T3E)
    @echo DEPRECATED
    @exit
endif


ifeq ($(TARGET),HPUX)
#
# HPUX 11.0
#
# removed reference to MLIB since 8.3 version of MLIB 
# does not support +ppu
#
    CORE_SUBDIRS_EXTRA = blas lapack
    MAKEFLAGS = -j 1 --no-print-directory
    CPP = /lib/cpp -P
    CC = cc
    FC = f90
    LDOPTIONS = -g -Wl,+vallcompatwarnings  +U77   
    LDOPTIONS +=  +DA2.0 +DS2.0 +O2
    LDOPTIONS +=   +O2
    LINK.f = f90   $(LDFLAGS) 
    CORE_LIBS +=  $(BLASOPT) -lnwclapack -lnwcblas   -lm
    FDEBUG = -g
    FOPTIONS =  +ppu -Wl,-a,archive
    COPTIONS = -Aa -D_HPUX_SOURCE +e 
    FOPTIMIZE = +O2 +Onolimit
    FOPTIMIZE += +DA2.0 +DS2.0a  +Odataprefetch  +Onofltacc +Onoinitcheck
    FOPTIMIZE += +Oprocelim +Oentrysched +Ofastaccess
    FVECTORIZE = +Oall +Onofltacc
    COPTIMIZE = -O
    RANLIB = echo

    DEFINES = -DHPUX -DEXTNAME
    ifeq ($(BUILDING_PYTHON),python)
        # needed if python was compiled with gcc (common)
        EXTRA_LIBS += -L/usr/local/lib/gcc-lib/hppa1.0-hp-hpux11.00/2.8.0 -lgcc
    endif
endif


ifeq ($(TARGET),HPUX64)
#
# HPUX 11.0
#
    CORE_SUBDIRS_EXTRA = blas lapack
    _CPU = $(shell uname -m  )
    MAKEFLAGS = -j 1 --no-print-directory
    CPP = /lib/cpp -P
    CC = cc
    FC = f90
    LDOPTIONS = -Wl,+vallcompatwarnings  +U77   
    CORE_LIBS +=  -lnwclapack $(BLASOPT) -lnwcblas  -lm
    CDEBUG =
    FDEBUG = -g
    FOPTIONS =  +ppu  #+U77  
    COPTIONS = -Aa -D_HPUX_SOURCE +e 
    ifeq ($(_CPU),ia64)
        FOPTIONS += +DD64 +DSitanium2 +Ofltacc=relaxed +Olibcalls +Onolimit +FPD
        COPTIONS += +DD64
        FOPTIMIZE = +O2
        FVECTORIZE = +Ofast  +O3 +Onoptrs_to_globals +Oloopblock 
        FDEBUG = +Ofast 
    else
        FOPTIONS +=  +DA2.0W 
        COPTIONS +=  +DA2.0W 
        FOPTIMIZE = +O2
        FVECTORIZE =  +Oall +Onofltacc
    endif
  
    COPTIMIZE = -O
    RANLIB = echo

    DEFINES = -DHPUX -DEXTNAME -DHPUX64 
    ifdef USE_I4FLAGS
    else
        FOPTIONS +=  +i8 
    endif
    DEFINES +=  -DEXT_INT
endif



ifeq ($(TARGET),IBM)
#
# IBM AIX
#

    CORE_SUBDIRS_EXTRA = lapack blas
    FC = xlf
    ifeq ($(FC),xlf)
        _FC=xlf
    endif

    CC = xlc
    ARFLAGS = urs
    RANLIB = echo
    MAKEFLAGS = -j 5 --no-print-directory
    INSTALL = @echo $@ is built
    CPP = /usr/lib/cpp -P

    FOPTIONS = -qEXTNAME -qnosave -qalign=4k -qxlf77=leadzero
#   -qinitauto=FF
    COPTIONS = 
#   -qstrict required with -O3 (according to Edo)
#   -qfloat=rsqrt gives faster square roots (off by -qstrict)
#   -qfloat=fltint gives faster real-integer conversion (off by -qstrict)
#   -qhot seems to break a lot of things so don\'t ever use it
    FOPTIMIZE = -O3 -qstrict -NQ40000 -NT80000 -qarch=auto -qtune=auto -NS2048

    ifdef RSQRT
        FOPTIMIZE  += -qfloat=rsqrt:fltint
    endif

    COPTIMIZE = -O -qarch=auto -qtune=auto

    ifdef  USE_GPROF
        FOPTIONS += -pg
        LDOPTIONS += -pg
    endif

    ifdef  USE_DEBUG
        FOPTIONS += -g
        LDOPTIONS += -g
    endif

    DEFINES = -DIBM -DAIX -DEXTNAME
    CORE_LIBS +=  $(BLASOPT) 
    ifdef USE_ESSL
        DEFINES += -DESSL
        CORE_LIBS += -lessl
    endif

    LIBPATH += -L/usr/lib 

    LDOPTIONS += -bmaxstack:0x60000000 -bmaxdata:0x60000000 -bloadmap:nwchem.lapi_map
    CORE_LIBS +=  -lnwclapack $(BLASOPT) -lnwcblas \
              -brename:.daxpy_,.daxpy \
              -brename:.dcopy_,.dcopy \
              -brename:.ddot_,.ddot \
              -brename:.dgemm_,.dgemm \
              -brename:.dgemv_,.dgemv \
              -brename:.dgesv_,.dgesv \
              -brename:.dgetrf_,.dgetrf \
              -brename:.dgetrs_,.dgetrs \
              -brename:.dlaset_,.dlaset \
              -brename:.dpotrf_,.dpotrf \
              -brename:.dpotri_,.dpotri \
              -brename:.dscal_,.dscal \
              -brename:.dspsvx_,.dspsvx \
              -brename:.idamax_,.idamax \
              -brename:.dswap_,.dswap \
              -brename:.dger_,.dger \
              -brename:.dtrsm_,.dtrsm \
              -brename:.dnrm2_,.dnrm2 \
              -brename:.dtrmm_,.dtrmm \
              -brename:.drot_,.drot \
              -brename:.dasum_,.dasum \
              -brename:.dtrmv_,.dtrmv \
              -brename:.dspmv_,.dspmv \
              -brename:.dspr_,.dspr \
              -brename:.dsyrk_,.dsyrk \
              -brename:.dsyr2k_,.dsyr2k \
              -brename:.dsymv_,.dsymv \
              -brename:.lsame_,.lsame \
              -brename:.xerbla_,.xerbla \
              -brename:.zgemm_,.zgemm \
              -brename:.dsyr2_,.dsyr2 \
              -brename:.dznrm2_,.dznrm2 \
              -brename:.zaxpy_,.zaxpy \
              -brename:.zcopy_,.zcopy \
              -brename:.zdotc_,.zdotc \
              -brename:.zdscal_,.zdscal \
              -brename:.zgemv_,.zgemv \
              -brename:.zgerc_,.zgerc \
              -brename:.zhemv_,.zhemv \
              -brename:.zher2_,.zher2 \
              -brename:.zher2k_,.zher2k \
              -brename:.zscal_,.zscal \
              -brename:.zswap_,.zswap \
              -brename:.ztrmm_,.ztrmm \
              -brename:.ztrmv_,.ztrmv \
              -brename:.izamax_,.izamax
#              -brename:.zherk_,.zherk \
#              -brename:.zhemm_,.zhemm \
#              -brename:.ztrsm_,.ztrsm \
#              -brename:.dtrsv_,.dtrsv \
#              -brename:.ztrsv_,.ztrsv 
#              -brename:.dsymm_,.dsymm \
#

  ##comment out from dtrmm_ inclusive

    EXPLICITF = TRUE
#
endif



ifeq ($(TARGET),IBM64)
# 
# IBM AIX 64-bit
# tested on ecs1 May 10 2000 AIX 4.3.3.10
# does not run on AIX  4.3.2.1 (skunkworks)
#
    ifeq ($(LAPACK_LIB),)
        CORE_SUBDIRS_EXTRA += lapack
    endif

    ifeq ($(BLASOPT),)
        CORE_SUBDIRS_EXTRA += blas
    endif

    FC = xlf
    ifeq ($(FC),xlf)
        _FC=xlf
    endif

    CC = xlc
    AR = ar -X 64
    RANLIB = echo
    MAKEFLAGS = -j 11 --no-print-directory
    INSTALL = @echo $@ is built
    CPP = /usr/lib/cpp -P

    FOPTIONS = -qEXTNAME -qnosave -qalign=4k -q64 -qxlf77=leadzero
    COPTIONS = -q64
    FOPTIMIZE = -O3 -qstrict -NQ40000 -NT80000  -qarch=auto -qtune=auto
    RSQRT=y
    FDEBUG = -O2 -qmaxmem=8192 -qsuppress=1500-030

    ifdef RSQRT
        FOPTIMIZE  += -qfloat=rsqrt:fltint
    endif

    XLF8= $(shell /usr/bin/lslpp -l xlfcmp  2>&1|grep COMM|head -n 1| awk ' / [8-9]./  {print "Y"};/[ ][1][0-9]./  {print "Y"}')

    ifdef XLF8
        FVECTORIZE= -O3 -qstrict -qtune=auto -qarch=auto -qcache=auto -qalign=natural \
                     -qnozerosize -qlargepage -qnozerosize -qipa=level=2
#       FOPTIMIZE = -O4  -NQ40000 -NT80000  -qarch=auto -qtune=auto
        FOPTIMIZE = -O3   -qarch=auto -qtune=auto
        #adding -qstrict to fix linking problem for v10.1
        FOPTIMIZE  += -qstrict# -qipa -qhot -qlargepage -qessl 
        FOPTIONS += -blpdata  
        FOPTIMIZE  += -qfloat=rsqrt:fltint
        FVECTORIZE  += -qfloat=rsqrt:fltint
    endif
    COPTIMIZE = -O -qmaxmem=8192 -qsuppress=1500-030

    DEFINES = -DIBM -DAIX -DEXTNAME
    DEFINES += -DCHKUNDFLW
    LIBPATH += -L/usr/lib -L/lib 

    ifdef USE_I4FLAGS
        FOPTIONS += -qintsize=4
    else
        FOPTIONS += -qintsize=8 
    endif

    DEFINES  += -DEXT_INT 

    ifdef  USE_GPROF
        FOPTIONS += -pg
        LDOPTIONS += -pg
    endif

    ifdef USE_DEBUG
        FOPTIONS += -g
        LDOPTIONS += -g
    endif

    LDOPTIONS += -bloadmap:nwchem.ibm64map -bbigtoc # bigtoc requires bmaxdata
    LDOPTIONS += -bmaxstack:0x80000000 -bmaxdata:0x200000000 # this limits MA to 8GB
    ifeq ($(LAPACK_LIB),)
        CORE_LIBS += -lnwclapack
    endif
    CORE_LIBS += $(BLASOPT)
    ifeq ($(BLASOPT),)
        CORE_LIBS += -lnwcblas
    endif
    XLFBREN = y
    EXPLICITF = TRUE
#
endif



ifeq ($(TARGET),LAPI)
#
    CORE_SUBDIRS_EXTRA = lapack blas
    FC = mpxlf_r 
    ifdef OLDXLF
        FC += -qnohpf
    endif

    CC = mpcc_r
    ARFLAGS = urs
    RANLIB = echo
    MAKEFLAGS = -j 1 --no-print-directory
    INSTALL = @echo $@ is built
    CPP = /usr/lib/cpp -P
    MPILIB = 
    LARGE_FILES = YES

     LDOPTIONS = -lc_r -lxlf90_r -lm_r -qEXTNAME -qnosave -qalign=4k -g 
     LINK.f   = mpcc_r   $(LDFLAGS)
     FOPTIONS = -qEXTNAME -qnosave -qalign=4k  -qxlf77=leadzero
     ifdef  USE_GPROF
         FOPTIONS += -pg
         LDOPTIONS += -pg
     endif

     ifdef  USE_DEBUG
         FOPTIONS += -g
         LDOPTIONS += -g
     endif

     # -qinitauto=7F # note that grad_force breaks with this option
     COPTIONS = 
     FOPTIMIZE = -O3 -qstrict -NQ40000 -NT80000 -NS10000 -qipa=level=2
     ifdef RSQRT
         FOPTIMIZE  += -qfloat=rsqrt:fltint
     endif
     COPTIMIZE = -O
     FC += -qarch=auto -qtune=auto -qcache=auto -qthreaded
     CC += -qarch=auto -qtune=auto -qcache=auto


     DEFINES = -DLAPI -DSP1 -DAIX -DEXTNAME
     DEFINES += -DCHKUNDFLW

     USE_ESSL = YES
     ifdef USE_ESSL
         CORE_LIBS +=  -lessl
     endif
     ifdef USE_PESSL
         DEFINES += -DESSL
         CORE_LIBS += -lpessl -lblacs 
     endif
     # Need ESSL before our own BLAS library but still need our
     # own stuff for misc. missing routines

     CORE_LIBS +=  -lnwclapack -lnwcblas

     LDOPTIONS += -bloadmap:nwchem.lapimap -bbigtoc
     LDOPTIONS += -bmaxstack:0x60000000 -bmaxdata:0x60000000 # needed because of bigtoc

     EXPLICITF = TRUE
#
endif



ifeq ($(TARGET),LAPI64)
#
    CORE_SUBDIRS_EXTRA = lapack blas
    FC = mpxlf_r 
    CC = mpcc_r
    ARFLAGS = urs
    RANLIB = echo
    MAKEFLAGS = -j 3 --no-print-directory
    INSTALL = @echo $@ is built
    CPP = /usr/lib/cpp -P
    MPILIB = 
    LARGE_FILES = YES

    LDOPTIONS = -lc_r -lxlf90_r -lm_r -qEXTNAME -qnosave -q64  -bloadmap:nwchem.lapi64_map $(LAPI64LIBS)
    LINK.f   = mpxlf_r   $(LDFLAGS)

    FOPTIONS = -qEXTNAME -qnosave -q64 -qalign=4k -qxlf77=leadzero -qthreaded
    AR   = ar -X 64
    COPTIONS = -q64
    FOPTIMIZE = -O3 -qstrict -NQ40000 -NT80000
    FOPTIMIZE += -qarch=auto -qtune=auto -qcache=auto
    ifdef RSQRT
        FOPTIMIZE  += -qfloat=rsqrt:fltint
    endif
    COPTIMIZE = -O
    XLF8= $(shell xlf -qversion  2>&1|grep Version|head -1| awk ' / [8-9]./  {print "Y"};/[ ][1][0-9]./  {print "Y"}')
    XLF10 = $(shell xlf -qversion  2>&1|grep Version|head -1| awk ' / 10./ {print "Y"}')
    XLF11 = $(shell xlf -qversion  2>&1|grep Version|head -1| awk ' / 11./ {print "Y"}')

    ifdef XLF8
        FVECTORIZE= -O3 -qstrict -qtune=auto -qarch=auto -qcache=auto -qalign=natural \
                    -qnozerosize -qlargepage -qnozerosize -qipa=level=2
        FOPTIMIZE = -O4  -NQ40000 -NT80000  -qarch=auto -qtune=auto
        FOPTIMIZE  += -qipa -qhot -qlargepage -qessl 
        FOPTIONS += -blpdata  
        FOPTIMIZE  += -qfloat=rsqrt:fltint
        FVECTORIZE  += -qfloat=rsqrt:fltint
    endif

    DEFINES = -DLAPI64 -DEXTNAME -DLAPI -DSP1 -DAIX
    DEFINES += -DCHKUNDFLW
    ifdef USE_I4FLAGS
        FOPTIONS += -qintsize=4
    else
        FOPTIONS += -qintsize=8
    endif

    DEFINES += -DEXT_INT
    CORE_LIBS +=  $(BLASOPT) -lnwclapack -lnwcblas
    LDOPTIONS += -bloadmap:nwchem.lapi64map -bbigtoc
    LDOPTIONS += -bmaxstack:0x80000000 -bmaxdata:0x80000000 # needed because of bigtoc
    # LDOPTIONS += -bmaxstack:0xe0000000 -bmaxdata:0xe0000000 # this for large memory
    XLFBREN = y

    EXPLICITF = TRUE
#
endif



ifeq ($(TARGET),MACX)
    FC = gfortran
    _FC = gfortran
    #
    # MacOSX 
    #
    ifdef USE_VECLIB
        CORE_SUBDIRS_EXTRA =  blas
    else
        CORE_SUBDIRS_EXTRA =  blas lapack
    endif

    _CPU = $(shell machine  )
    FC = gfortran
    INSTALL = @echo nwchem is built
    RANLIB = ranlib
    MAKEFLAGS = -j 1 --no-print-directory
    DEFINES =-DMACX
    COPTIONS = -m32
    FOPTIONS   = -m32
    CFLAGS_FORGA = -m32
    FFLAGS_FORGA = -m32

    # required for mpich2 3.x and clang
    DEFINES +=-DMPICH_NO_ATTR_TYPE_TAGS
    CFLAGS_FORGA +=-DMPICH_NO_ATTR_TYPE_TAGS

    ifeq ($(FC),xlf)
        _FC=xlf
        XLFMAC=y
        FOPTIONS = -qextname -qfixed -qnosave  -qalign=4k
        FOPTIONS +=  -NQ40000 -NT80000 -NS2048 -qmaxmem=8192 -qsuppress=1500-030 -qxlf77=leadzero
        FOPTIMIZE= -O3 -qstrict  -qarch=auto -qtune=auto -qcache=auto -qcompact

        ifdef RSQRT
            FOPTIMIZE  += -qfloat=rsqrt:fltint
        endif

        FVECTORIZE += $(FOPTIMIZE) -qunroll=yes 
        FDEBUG= -O2 -qcompact 
        DEFINES  +=-DXLFLINUX -DCHKUNDFLW
        FOPTIONS += $(INCLUDES) -WF,"$(DEFINES)" $(shell echo $(LIB_DEFINES) | sed -e "s/-D/-WF,-D/g"   | sed -e 's/\"/\\\"/g'  | sed -e "s/\'/\\\'/g")
    endif

    ifeq ($(FC),g77)
        #g77, only decent one form Fink http://fink.sf.net
        #gcc version 3.4 20031015 (experimental)
        _G77V33= $(shell g77 -v  2>&1|egrep spec|head -n 1|awk ' /3.3/  {print "Y"}')
        FDEBUG= -O1 -g
        FOPTIONS   = -fno-second-underscore -fno-globals -Wno-globals
        FOPTIMIZE  = -O3 -fno-inline-functions -funroll-loops
        FOPTIMIZE += -falign-loops=16 -falign-jumps=16 -falign-functions=16
        FOPTIMIZE += -ffast-math -mpowerpc-gpopt
        FOPTIMIZE += -maltivec

        ifeq ($(_G77V33),Y)
            FOPTIONS += -fno-force-mem 
            FOPTIMIZE += -fno-force-mem 
        endif

        ifeq ($(_CPU),ppc970)
            #G5
            FOPTIMIZE += -mtune=970 -mcpu=970 -mpowerpc64
        endif

        ifeq ($(_CPU),ppc7450)
#G4
            FOPTIMIZE += -mtune=7450 -mcpu=7450
        endif
    endif

    FDEBUG = -g -O1

    ifeq ($(FC),gfortran)
        _FC=gfortran
        #gcc version 4.2.0 200512 (experimental)
        LINK.f = gfortran -m32  $(LDFLAGS) 
        FDEBUG = -O0 -g
        FOPTIMIZE  = -O2 -ffast-math
        FOPTIMIZE  += -Wuninitialized -Wno-maybe-uninitialized 
        DEFINES  += -DGFORTRAN
        GNUMAJOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | grep __GNUC__ |cut -c18-)

        ifdef GNUMAJOR
            GNUMINOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __GNUC_MINOR | cut -c24)
            GNU_GE_4_6 = $(shell [ $(GNUMAJOR) -gt 4 ] || [ $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 6 ] && echo true)
            GNU_GE_4_8 = $(shell [ $(GNUMAJOR) -gt 4 ] || [ $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 8 ] && echo true)
            GNU_GE_6 = $(shell [ $(GNUMAJOR) -ge 6  ] && echo true)
        endif

        ifeq ($(GNU_GE_4_6),true)
            DEFINES  += -DGCC46
        endif

        ifeq ($(GNU_GE_4_8),true)
            FDEBUG += -fno-aggressive-loop-optimizations
            FOPTIMIZE +=-fno-aggressive-loop-optimizations
            FOPTIONS +=-fno-aggressive-loop-optimizations
            FFLAGS_FORGA += -fno-aggressive-loop-optimizations
          
            ifeq ($(V),-1)
                FOPTIONS += -w
            else
                FOPTIONS += -Warray-bounds
            endif
        endif

        ifeq ($(GNU_GE_6),true)
            FOPTIMIZE += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
            FOPTIONS += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
            FDEBUG += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
        endif

        ifdef USE_OPENMP
            FOPTIONS  += -fopenmp
            LDOPTIONS += -fopenmp
        endif

        ifeq ($(_CPU),ppc970)
            #G5
            FVECTORIZE = -ffast-math  -O2 -ftree-vectorize 
            FVECTORIZE += -ftree-vectorizer-verbose=1
            FOPTIMIZE += -mtune=970 -mcpu=970 -mpowerpc64
            FVECTORIZE += -mtune=970 -mcpu=970 -mpowerpc64
        endif

        ifeq ($(_CPU),ppc7450)
#G4
            FVECTORIZE = -ffast-math  -O2 -ftree-vectorize 
            FVECTORIZE += -ftree-vectorizer-verbose=1
            FOPTIMIZE  += -fprefetch-loop-arrays #-ftree-loop-linear
            FOPTIMIZE += -mtune=7450 -mcpu=7450
        endif

        ifeq ($(_CPU),i486)
            #gcc version 4.2.0 200608 (experimental)
#           FOPTIONS= -malign-double# this break with gfort 4.2 and later http://gcc.gnu.org/bugzilla/show_bug.cgi?id=29562
            FOPTIMIZE+= -funroll-all-loops -mtune=native 
            FVECTORIZE=-O3 -ffast-math -mtune=native -mfpmath=sse -msse3 -ftree-vectorize -ftree-vectorizer-verbose=1   -fprefetch-loop-arrays  -funroll-all-loops 
#           FOPTIMIZE=-O1
#           FVECTORIZE=-O1
        endif

        ifdef USE_F2C
            #possible segv with use of zdotc (e.g. with GOTO BLAS)
            #http://gcc.gnu.org/bugzilla/show_bug.cgi?id=20178
            FOPTIONS +=  -ff2c -fno-second-underscore
        endif

        DEFINES  += -DCHKUNDFLW -DGCC4
    endif
    # End of gfortran


    ifeq ($(FC),ifort)
        _FC=ifort
        #ifort 9.1
#       LINK.f = ifort  $(LDFLAGS) 
        FOPTIONS   += -align    -mp1 -w -g #-vec-report1

        ifdef  USE_GPROF
            FOPTIONS += -qp
        endif

        ifdef  USE_DEBUG
            FOPTIONS += -g
        endif

        FOPTIMIZE = -O3 -prefetch  -unroll 
        FDEBUG=-O0 -g
        DEFINES   += -DIFCLINUX
    endif


    ifeq ($(CC),xlc)
        COPTIONS  +=  -qlanglvl=extended
    else
        COPTIONS   += -Wall #-no-cpp-precomp
        COPTIMIZE  = -g -O2
    endif

    ifdef  USE_GPROF
        FOPTIONS += -pg
        LDOPTIONS += -pg
        COPTIONS += -pg
    endif

    ifdef  USE_DEBUG
        FOPTIONS += -g
        LDOPTIONS += -g
        COPTIONS += -g
    endif

    ifdef USE_VECLIB
        CORE_LIBS += $(BLASOPT)  -Wl,-framework -Wl,vecLib -lnwcblas
    else
        CORE_LIBS +=   -lnwclapack $(BLASOPT)  -lnwcblas
    endif

    ifeq ($(FC),xlf) 
        LDOPTIONS = -Wl,-multiply_defined -Wl,warning
    else
#       _GCC4= $(shell gcc -v  2>&1|egrep spec|head -n 1|awk ' / 3./  {print "N";exit}; / 2./ {print "N";exit};{print "Y"}')
        _GCC4= $(shell $(CC) -dM -E - < /dev/null | egrep __VERS | cut -c22|awk ' /3/  {print "N";exit}; /2/ {print "N";exit};{print "Y"}')
        ifeq ($(_GCC4),Y) 
#           EXTRA_LIBS += 
        else
            EXTRA_LIBS += -lm -lcc_dynamic
        endif
    endif
endif
#


_V104=$(shell uname -v 2>&1|awk ' /Version 7./ {print "N";exit}; /Version 8./ {print "Y";exit}')
ifeq ($(_V104),Y)
    EXTRA_LIBS +=-lSystemStubs
endif



ifeq ($(TARGET),MACX64)
#
# MacOSX 64bit
#
    ifeq ($(FC),f77)
        FC = gfortran
        _FC = gfortran
    endif

    ifeq ($(shell $(CNFDIR)/strip_compiler.sh $(FC)),gfortran)
        _FC := gfortran
    endif

    ifdef USE_VECLIB
        ifdef USE_64TO32
            CORE_SUBDIRS_EXTRA =
        else
            vecliberr:
	        @echo The Apple supplied vector math library does not support 8-byte integers
	        @echo You must also set USE_64TO32 and do a "make 64_to_32" to change the source code
	        @exit 1
        endif
    else
        ifeq ($(BLASOPT),)
            CORE_SUBDIRS_EXTRA += blas
        endif
        ifeq ($(LAPACK_LIB),)
            CORE_SUBDIRS_EXTRA += lapack
        endif
    endif

    INSTALL = @echo nwchem is built
    RANLIB = ranlib
#   MAKEFLAGS = -j 1 --no-print-directory
    DEFINES   = -DMACX
    DEFINES  += -DEXT_INT
    LINK.f    = $(FC) $(LDFLAGS) -Wl,-flat_namespace
    GOTCLANG= $(shell $(CC) -dM -E - </dev/null 2> /dev/null |grep __clang__|head -1|cut -c19)
    ifeq ($(GOTCLANG),1)
        COPTIONS   += -fPIC
    endif

    ifeq ($(_FC),gfortran)
        #gcc version 
        FOPTIONS  = -cpp #-Wextra #-Wunused #-ffast-math

        ifdef USE_I4FLAGS
        else
            FOPTIONS += -fdefault-integer-8
        endif

        FOPTIMIZE = -O2 -ffast-math
        ifeq ($(V),-1)
            FOPTIONS += -w
        else
            FOPTIMIZE  += -Wuninitialized -Wno-maybe-uninitialized
        endif

        DEFINES   += -DGFORTRAN -DGCC4
#
        FOPTIMIZE+= -funroll-all-loops
        ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
            FOPTIMIZE+= -mtune=native
        endif

        #FVECTORIZE=-O3 -ffast-math -mtune=native -mfpmath=sse -msse3 -ftree-vectorize -ftree-vectorizer-verbose=1   -fprefetch-loop-arrays  -funroll-all-loops
        FVECTORIZE=-O3 -ffast-math -ftree-vectorize -ftree-vectorizer-verbose=1 -funroll-all-loops
        ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
            FVECTORIZE+= -mtune=native
        endif

        GNUMAJOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | grep __GNUC__ |cut -c18-)

        ifneq ($(strip $(GNUMAJOR)),)
            GNUMINOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __GNUC_MINOR | cut -c24)
            GNU_GE_4_6 = $(shell [ $(GNUMAJOR) -gt 4 -o \( $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 6 \) ] && echo true)
            GNU_GE_4_8 = $(shell [ $(GNUMAJOR) -gt 4 -o \( $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 8 \) ] && echo true)
            GNU_GE_6 = $(shell [ $(GNUMAJOR) -ge 6  ] && echo true)
            GNU_GE_8 = $(shell [ $(GNUMAJOR) -ge 8  ] && echo true)

            ifeq ($(GNU_GE_4_6),true)
                DEFINES  += -DGCC46
            endif

            ifeq ($(GNU_GE_4_8),true)
                FDEBUG += -fno-aggressive-loop-optimizations
                FOPTIMIZE +=-fno-aggressive-loop-optimizations
                FOPTIONS +=-fno-aggressive-loop-optimizations
                FFLAGS_FORGA += -fno-aggressive-loop-optimizations
                ifeq ($(V),-1)
                    FOPTIONS += -w
                else
                    FOPTIONS += -Warray-bounds
                endif
            endif # GNU_GE_4_8

            ifeq ($(GNU_GE_6),true)
                FOPTIMIZE += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
                FOPTIONS += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
                FDEBUG += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
            endif
        endif # GNUMAJOR

        ifeq ($(GNU_GE_8),true)
            FOPTIONS   += -std=legacy
        endif

        ifdef USE_OPENMP
            FOPTIONS  += -fopenmp
            LDOPTIONS += -fopenmp
        endif

        ifdef  USE_ASAN
            FOPTIONS += -fsanitize=address -fsanitize-recover=address
            LDOPTIONS += -fsanitize=address -fsanitize-recover=address
        endif

        ifdef  USE_FPE
            FOPTIONS += -ffpe-trap=invalid,zero,overflow  -fbacktrace
        endif
    endif # gfortran

    ifdef  USE_GPROF
        FOPTIONS += -pg
        LDOPTIONS += -pg
        COPTIONS += -pg
    endif

    ifdef  USE_DEBUG
        FOPTIONS += -g
        LDOPTIONS += -g
        COPTIONS += -g
    endif

    ifdef USE_VECLIB
        CORE_LIBS += $(BLASOPT)  -Wl,-framework -Wl,vecLib -lnwcblas
    else
        ifeq ($(LAPACK_LIB),)
            CORE_LIBS += -lnwclapack
        endif
        CORE_LIBS +=    $(BLASOPT)
        ifeq ($(BLASOPT),)
            CORE_LIBS += -lnwcblas
        endif
    endif


    ifeq ($(FC),ifort)
        _IFCV11= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 11) {print "Y";exit}}')
        _IFCV12= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 12) {print "Y";exit}}')
        _IFCV14= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 14) {print "Y";exit}}')
        _IFCV15ORNEWER=$(shell ifort -logo  2>&1|egrep "Version "|head -n 1 | sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 15) {print "Y";exit}}')
        _IFCV17=$(shell ifort -logo  2>&1|egrep "Version "|head -n 1 | sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 17) {print "Y";exit}}')
        DEFINES  += -DIFCV8 -DIFCLINUX

        ifdef USE_I4FLAGS
        else
            FOPTIONS += -i8 
        endif

        FOPTIONS += -fpp -g -no-save-temps
        FDEBUG    = -O2 -g
        FOPTIMIZE = -O3

        ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
            FOPTIMIZE += -xHost
        endif

        ifdef USE_OPENMP
            ifeq ($(_IFCV15ORNEWER), Y)
                FOPTIONS  += -qopenmp
                LDOPTIONS += -qopenmp
            else
                FOPTIONS  += -openmp
                LDOPTIONS += -openmp
            endif
        endif

        ifdef  USE_FPE
            FOPTIONS += -fpe0 -traceback #-fp-model  precise
        endif

        ifdef USE_NOSIMD
            FOPTIONS  += -no-simd
        endif

        ifeq ($(_IFCV11),Y) 
            # The MKL option is available since Intel 11.1.
            # TODO: Test minor version number above when major=11.
            # https://software.intel.com/en-us/articles/using-mkl-in-intel-compiler-mkl-qmkl-options
            ifndef BLASOPT
                ifndef USE_INTERNALBLAS
                # When user requests OpenMP, MKL should use it, otherwise not (to avoid oversubscription).
                    ifdef USE_OPENMP
                        BLASOPT = -mkl=parallel
                    else
                        BLASOPT = -mkl=sequential
                    endif
                endif
            endif
        endif

        ifeq ($(_IFCV11),Y) 
            #next 2 lines needed for fp accuracy
            FDEBUG += -fp-model source
            ifeq ($(_IFCV12),Y) 
                FOPTIONS += -fimf-arch-consistency=true
            endif
        endif
        ifeq ($(V),-1)
            FOPTIONS += -diag-disable=7713,8291,15009
        endif
    endif #ifort



    ifdef USE_CCDYNAMIC
        EXTRA_LIBS += -lm -lcc_dynamic
    endif

#   required for mpich2 3.x and clang
    COPTIONS +=-DMPICH_NO_ATTR_TYPE_TAGS
#

endif




ifeq ($(TARGET),$(findstring $(TARGET),LINUX CYGNUS CYGWIN))
#
#
# Linux or Cygwin under Windows running on an x86 using g77
#
    NICE = nice -n 2
    SHELL := $(NICE) /bin/sh

    ifeq ($(BLASOPT),)
        CORE_SUBDIRS_EXTRA += blas
    endif

    ifeq ($(LAPACK_LIB),)
        CORE_SUBDIRS_EXTRA += lapack
    endif

    CC = gcc
    RANLIB = ranlib
    MAKEFLAGS = -j 1 --no-print-directory
    INSTALL = @echo $@ is built
    CPP = gcc -E -nostdinc -undef -P
    FCONVERT = (/bin/cp $< /tmp/$$$$.c; \
               $(CPP) $(CPPFLAGS) /tmp/$$$$.c | sed '/^$$/d' > $*.f; \
               /bin/rm -f /tmp/$$$$.c) || exit 1

    FC=gfortran
    ifeq ($(FC),f77)
        FC = gfortran
        _FC = gfortran
    endif

    ifeq ($(shell $(CNFDIR)/strip_compiler.sh $(FC)),gfortran)
        _FC := gfortran
    endif
    
    ifeq ($(FC),$(findstring $(FC),i686-w64-mingw32.static-gfortran))
        _FC = gfortran
    endif

    ifeq ($(shell $(CNFDIR)/strip_compiler.sh $(CC)),gcc)
        ifneq ($(CC),cc)
            _CC = gcc
        endif
    endif

    ifeq ($(CC),$(findstring $(CC),i686-w64-mingw32.static-gcc))
        ifneq ($(CC),cc)
            _CC = gcc
        endif
    endif

    LINUXCPU = $(shell uname -m |\
               awk ' /sparc/ { print "sparc" }; /i*86/ { print "x86" };  /ppc*/ { print "ppc"};  /arm*/ { print "arm"}; /mips*/ { print "mips"} ' )

    GOTMINGW32= $(shell $(CC) -dM -E - </dev/null 2> /dev/null |grep MINGW32|cut -c21)
    ifeq ($(GOTMINGW32),1)
        ifdef USE_OPENMP 
            errorompming:
	        @echo 
	        @echo "  MinGW environment does not support OpenMP"
	        @echo "  Please unset USE_OPENMP"
	        @echo 
        endif
    endif
  

    DEFINES = -DLINUX

    FOPTIMIZE  = -O2 
    COPTIONS   = -Wall
    COPTIMIZE  = -g -O2

    ifeq ($(_FC),gfortran)
        FOPTIONS   = # -Wextra -Wunused  
        FOPTIMIZE  += -ffast-math

        ifeq ($(V),-1)
            FOPTIONS += -w
        else
            FOPTIMIZE  += -Wuninitialized -Wno-maybe-uninitialized
        endif

        DEFINES  += -DGFORTRAN
        GNUMAJOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | grep __GNUC__ |cut -c18-)

        ifdef GNUMAJOR
            GNUMINOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __VERS | cut -c24)
            GNU_GE_4_6 = $(shell [ $(GNUMAJOR) -gt 4 ] || [ $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 6 ] && echo true)
            GNU_GE_4_8 = $(shell [ $(GNUMAJOR) -gt 4 ] || [ $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 8 ] && echo true)
            GNU_GE_6 = $(shell [ $(GNUMAJOR) -ge 6  ] && echo true)
            GNU_GE_8 = $(shell [ $(GNUMAJOR) -ge 8  ] && echo true)

            ifeq ($(GNU_GE_4_6),true)
                ifdef  USE_FPE
                    FOPTIONS += -ffpe-trap=invalid,zero,overflow  -fbacktrace
                endif
                DEFINES  += -DGCC46
            endif

            ifeq ($(GNU_GE_4_8),true)
                FDEBUG +=-O2 -g -fno-aggressive-loop-optimizations
                FOPTIMIZE +=-fno-aggressive-loop-optimizations
                FFLAGS_FORGA += -fno-aggressive-loop-optimizations
            endif

            ifeq ($(GNU_GE_6),true)
                FOPTIMIZE += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
                FOPTIONS += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
                FDEBUG += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
            endif
        endif
    endif

    ifeq ($(GNU_GE_8),true)
        FOPTIONS   += -std=legacy
    endif

    ifdef USE_OPENMP
        FOPTIONS  += -fopenmp
        LDOPTIONS += -fopenmp
    endif


    ifeq ($(LINUXCPU),x86) 
        ifeq ($(TARGET),CYGNUS)
            DEFINES += -DCYGNUS
        endif

        ifeq ($(TARGET),CYGWIN)
            DEFINES += -DCYGWIN -DCYGNUS
        endif
  
        _CPU = $(shell uname -m  )

        ifeq ($(FC),g77)
            _G77V33= $(shell g77 -v  2>&1|egrep spec|head -n 1|awk ' /3.3/  {print "Y"}')
            FOPTIONS   += -fno-second-underscore   
            FOPTIONS   += -fno-f90  -ffixed-line-length-72 -ffixed-form
            FOPTIMIZE  +=  -O2  -malign-double -finline-functions 
            COPTIONS   += -Wall  -malign-double 
            COPTIMIZE  += -g -O2
            FOPTIONS  +=  -malign-double -fno-globals -Wno-globals  -fno-silent #-Wunused  
            FOPTIMIZE += -Wuninitialized
            FOPTIMIZE += -ffast-math -funroll-loops -fstrength-reduce
            FOPTIMIZE += -fno-move-all-movables -fno-reduce-all-givs 
            FOPTIMIZE += -fforce-addr 
#           see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=13037
#           for atomscf/orderd.f  (bug report by Kirill Smelkov)
            ifeq ($(_G77V33),Y)
                FOPTIONS += -fno-force-mem 
                FOPTIMIZE += -fno-force-mem 
            else
                FOPTIMIZE += -fforce-mem 
            endif
        endif


        ifeq ($(GOTMINGW32),1)
            _CPU=i786
        else
            ifeq ($(_CPU),i686)
                _GOTSSE2= $(shell cat /proc/cpuinfo | egrep sse2 | tail -n 1 | awk ' /sse2/  {print "Y"}')
                ifeq ($(_GOTSSE2),Y) 
                    _CPU=i786
                endif
            endif
        endif


        ifeq ($(_CPU),i786)
#           COPTIONS   =  -march=i686 
            ifdef USE_GCC31
                FDEBUG=-O1 -g
                ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                    COPTIMIZE +=-march=pentium4 -mcpu=pentium4 #-msse2 -mfpmath=sse 
                    FOPTIMIZE +=-march=pentium4 -mcpu=pentium4# -msse2 -mfpmath=sse
                endif
            else
                COPTIONS   = -Wall -malign-double 
            endif
        else
        endif

        ifeq ($(_CPU),k7)
            FOPTIONS   = -fno-second-underscore  
            COPTIONS   = -Wall -malign-double
            ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                ifdef  USE_GCC31
                    FOPTIONS += -march=athlon
                    COPTIONS += -march=athlon
                else
                    FOPTIONS += -march=k6
                    COPTIONS += -march=k6
                endif
            endif
        endif

        ifeq ($(FC),pgf77)
            DEFINES   += -DPGLINUX
#           added -Kieee to get dlamc1 to work on pgf77 3.1-3 EA Jun 8th 2000
            FOPTIONS   = -Mdalign -Minform,warn -Mnolist -Minfo=loop -Munixlogical -Kieee
            FOPTIONS   = -Mdalign -Minform,warn -Mnolist -Munixlogical -Kieee
            ifeq ($(V),-1)
                 FOPTIONS += -Minform,fatal
            else
                 FOPTIONS += -Minform,warn
            endif
            ifdef USE_OPTREPORT
                 FOPTIONS += -Minfo=loop
            endif
            ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                ifeq ($(_CPU),i586)
                    FOPTIONS  += -tp p5  
                endif
                ifeq ($(_CPU),i686)
                    FOPTIONS  += -tp p6 -Mvect=prefetch
                endif
                ifeq ($(_CPU),i786)
                    FOPTIONS  += -tp piv  -Mcache_align  -Mvect=prefetch
                endif
            endif
            FOPTIMIZE  = -O2 -Mvect=assoc,cachesize:262144 -Munroll -Mnoframe
        endif

#       _FC=g77
#
        ifeq ($(FC),ifc)
            _FC=ifc
        endif

        ifeq ($(FC),ifort)
            _FC=ifc
        endif

        ifeq ($(_FC),ifc)
            FOPTIONS   =  -align    -mp1 -w -g #-vec-report1
            ifdef  USE_GPROF
                FOPTIONS += -qp
            endif
            ifdef  USE_DEBUG
                FOPTIONS += -g
            endif
            _IFCV7= $(shell ifort -v  2>&1|egrep "Version "|head -n 1|awk '/7./ {print "Y"; exit}')
            ifneq ($(_IFCV7),Y)
                DEFINES+= -DIFCV8
                ifeq ($(FC),ifc)
                    FOPTIONS += -quiet
                endif
            endif
            ifdef  USE_FPE
                FOPTIONS += -fpe-all=0 -traceback #-fp-model  precise
            endif

            FOPTIMIZE = -O3 -prefetch  -unroll

            ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                ifeq ($(_CPU),i586)
                    FOPTIMIZE +=  -tpp5 -xi # this are for PentiumII
                endif

                ifeq ($(_CPU),k7)
                    FOPTIMIZE +=  -xM  # this are for Athlon
                endif
                ifeq ($(_CPU),i686)
                    FOPTIMIZE +=  -tpp6 -xK   # this are for PentiumIII
                endif
                ifeq ($(_CPU),i786)
                    FOPTIMIZE +=  -tpp7 -xW    # this are for PentiumIV
                endif
            endif

            DEFINES   += -DIFCLINUX

            ifneq ($(_IFCV7),Y)
                FOPTIMIZE += -ansi_alias-
            endif
        endif


        ifeq ($(_FC),gfortran)
            LINK.f = gfortran  $(LDFLAGS) 

            ifneq ($(_CPU),arm)
                FOPTIONS  += -m32
                COPTIONS  += -m32
                CFLAGS_FORGA += -m32
                FFLAGS_FORGA += -m32
            endif

            FOPTIMIZE  += -O2 -ffast-math

            ifeq ($(V),-1)
                FOPTIONS += -w
            else
                FOPTIMIZE  += -Wuninitialized -Wno-maybe-uninitialized
            endif

            ifeq ($(_CPU),i786)
                ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                    FOPTIONS += -march=pentium4 -mtune=pentium4
                endif
                FVECTORIZE = $(FOPTIMIZE) -O3 -ftree-vectorize 
                FVECTORIZE += -ftree-vectorizer-verbose=1
#               FOPTIMIZE  += -fprefetch-loop-arrays -ftree-loop-linear
            else
                FOPTIONS +=  -ffloat-store
            endif

            ifdef USE_F2C
                #possible segv with use of zdotc (e.g. with GOTO BLAS)
                #http://gcc.gnu.org/bugzilla/show_bug.cgi?id=20178
                FOPTIONS +=  -ff2c -fno-second-underscore
            endif
            FDEBUG += -g -O0
            DEFINES  += -DCHKUNDFLW -DGCC4
        endif


        ifeq ($(CC),icc)
            COPTIONS   =   -mp1 -w -g #-vec-report1
            COPTIMIZE = -O3   -unroll 
            ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                ifeq ($(_CPU),i586)
                    COPTIMIZE +=  -tpp5 -xi # this are for PentiumII
                endif
                ifeq ($(_CPU),i686)
                    COPTIMIZE +=  -tpp6 -xK   # this are for PentiumIII
                endif
                ifeq ($(_CPU),i786)
                    COPTIMIZE +=  -tpp7 -xW   # this are for PentiumIV
                endif
            endif
        endif

    endif #x86


    ifeq ($(LINUXCPU),ppc)
#       this are for PowerPC
#       Tested on SLES 9
#       Feb 7th 2005
#       xlf v9.1
#       xlc v7.0 
#       gcc-3.2.3-42 
        ifeq ($(FC),xlf)
            _FC=xlf
        endif

        ifeq ($(FC),blrts_xlf)
            _FC=xlf
        endif

        ifeq ($(_FC),xlf)
            FOPTIONS  = -q32  -qextname -qfixed 
            FOPTIONS +=  -NQ40000 -NT80000 -NS2048 -qmaxmem=8192 -qsuppress=1500-030 -qxlf77=leadzero
            FOPTIMIZE= -O3 -qstrict -qfloat=fltint
            ifeq ($(FC),blrts_xlf)
                FOPTIMIZE+= -qarch=440 -qtune=440
            else
                FOPTIMIZE+= -qarch=auto -qtune=auto
            endif

            FDEBUG= -O2 -g
#           EXPLICITF = TRUE
            DEFINES  +=   -DXLFLINUX
            CPP=/usr/bin/cpp  -P -C -traditional
        endif

        ifeq ($(FC),g77)
            FOPTIONS   = -fno-second-underscore -fno-globals -Wno-globals
            FOPTIMIZE  = -g -O2
        endif

        ifeq ($(CC),xlc)
            _CC=xlc
        endif

        ifeq ($(CC),blrts_xlc)
            _CC=xlc  
        endif

        ifeq ($(_CC),xlc)
            COPTIONS  +=  -q32 -qlanglvl=extended
        else
            COPTIONS   = -Wall
            COPTIMIZE  = -g -O2
        endif
        LDOPTIONS += -Wl,--relax #-v
    endif

    LINK.f = $(FC) $(FOPTIONS) $(LDFLAGS) 
    ifeq ($(LINUXCPU),x86)
        ifeq ($(FC),pgf77)
            LDOPTIONS += -g -Wl,--export-dynamic
            EXTRA_LIBS += -lm
        else
            ifeq ($(_FC),ifc)
                ifneq ($(_IFCV7),Y)
                    EXTRA_LIBS +=  
                else
                    EXTRA_LIBS +=   -Vaxlib  
                endif
                ifeq ($(_CPU),i786)
                    EXTRA_LIBS +=  -lsvml
                endif
                EXTRA_LIBS += #-static
                LDOPTIONS = -g -Wl,--export-dynamic
            else
                ifeq ($(GOTMINGW32),1)
                    LDOPTIONS += -g -O0 
                    EXTRA_LIBS += -lwsock32
                else
                    LDOPTIONS = -Xlinker --export-dynamic 
#                   LDOPTIONS = --Xlinker -O -Xlinker -static
                    EXTRA_LIBS += -lm
                endif
            endif
        endif
    endif


#   EXTRA_LIBS +=-lefence # link against Electricfence
#   CORE_LIBS += -lnwclapack $(BLASOPT) -lnwcblas
# end of Linux, Cygnus
endif



ifneq ($(TARGET),LINUX)
    ifeq ($(TARGET),$(findstring $(TARGET),LINUX64 CYGWIN64 CATAMOUNT))

        GOTMINGW64=$(shell $(CC) -dM -E - </dev/null 2> /dev/null |grep MINGW64|cut -c21)
        ifeq ($(GOTMINGW64),1)
            _CPU = x86_64
        else
            _CPU = $(shell uname -m  )
        endif

#       ifeq ($(NWCHEM_TARGET),LINUX64)

        ifeq ($(FC),g77)
            g7764:
	        @echo 
	        @echo g77 not supported on 64-bit Linux
	        @echo please use supported gompilers
	        @echo 
	        @exit 1
        endif

        ifeq ($(FC),ftn)
            _FC=pgf90
            ifeq ($(PE_ENV),PGI)
                _FC=pgf90
                _CC=pgcc
            endif

            ifeq ($(PE_ENV),INTEL)
                _FC=ifort
                _CC=icc
            endif

            ifeq ($(PE_ENV),GNU)
                _FC=gfortran
                _CC=gcc
            endif

            ifeq ($(PE_ENV),CRAY)
                _FC=crayftn
#               _CC=craycc
#               as of 2021 cray cc is derived from clang
                _CC=clang
            endif

            ifeq ($(PE_ENV),AOCC)
                _FC=gfortran
                USE_FLANG=1
                _CC=clang
            endif

            ifeq ($(PE_ENV),NVIDIA)
#               nvfortran same as pgf90
                _FC=pgf90
#               nvcc same as pgcc
                _CC=pgcc
            endif
            DEFINES  += -DCRAYXT -DNOIO
            USE_NOIO=1
        endif


        ifeq ($(shell $(CNFDIR)/strip_compiler.sh $(CC)),gcc)
            _CC=gcc
        endif

        ifeq ($(CC),pgcc)
            _CC=pgcc
        endif

        ifeq ($(CC),nvcc)
            _CC=pgcc
        endif

        ifeq ($(CC),icc)
            _CC=icc
        endif

        ifeq ($(FC),pgf90)
            _FC=pgf90
        endif

        ifeq ($(shell $(CNFDIR)/strip_compiler.sh $(FC)),nvfortran)
            _FC := pgf90
        endif

        ifeq ($(FC),pgf77)
            _FC=pgf90
        endif

        ifeq ($(FC),ifc)
            _FC=ifort
        endif

        ifeq ($(FC),ifort)
            _FC=ifort
        endif

        ifeq ($(FC),ifx)
            USE_IFX=1
            _IFCV8=1
            _FC=ifort
        endif

        ifeq ($(shell $(CNFDIR)/strip_compiler.sh $(FC)),gfortran)
            _FC := gfortran
        endif

        ifeq ($(FC),$(findstring $(FC),i686-w64-mingw32.static-gfortran x86_64-w64-mingw32-gfortran-win32))
            _FC := gfortran
        endif

        ifeq ($(shell $(CNFDIR)/strip_compiler.sh $(CC)),gcc)
            ifneq ($(CC),cc)
                _CC := gcc
            endif
        endif

        ifeq ($(CC),$(findstring $(CC),i686-w64-mingw32.static-gcc x86_64-w64-mingw32-gcc-win32))
            ifneq ($(CC),cc)
                _CC= gcc
            endif
        endif

        ifeq ($(FC),armflang)
            _FC=armflang
            USE_FLANG=1
        endif

        ifeq ($(shell $(CNFDIR)/strip_compiler.sh $(FC)),flang)
            _FC=gfortran
            USE_FLANG=1
        endif

        ifeq ($(shell $(CNFDIR)/strip_compiler.sh $(FC)),amdflang)
            _FC=gfortran
            USE_FLANG=1
        endif

        ifeq ($(CC),clang)
            _CC=gcc
        endif

        ifeq ($(CC),amdclang)
            _CC=gcc
        endif

        ifeq ($(CC),icx)
            _CC=gcc
        endif

        ifeq ($(FC),$(findstring $(FC),xlf2008_r xlf_r xlf xlf90 xlf90_r))
            _FC=xlf
        endif

        ifndef _FC
            FC=gfortran
            _FC=gfortran
        endif

        ifndef _CC
            _CC=gcc
        endif

        FOPTIMIZE  = -O2 

        ifeq ($(_CPU),aarch64)
            DONTHAVEM64OPT=Y
        endif

        ifeq ($(_CPU),mips64)
            DONTHAVEM64OPT=Y
            COPTIONS   = -mabi=64
            FOPTIONS   = -mabi=64
            FFLAGS_FORGA   = -mabi=64
            CFLAGS_FORGA   = -mabi=64
        endif

        ifeq ($(_CPU),riscv64)
            DONTHAVEM64OPT=Y
            COPTIONS   =  -march=rv64gc -mabi=lp64d
            FOPTIONS   =  -march=rv64gc -mabi=lp64d
            FFLAGS_FORGA   = -march=rv64gc -mabi=lp64d
            CFLAGS_FORGA   = -march=rv64gc -mabi=lp64d
        endif

        ifeq ($(_CC),$(findstring $(_CC),gcc clang))
            ifneq ($(DONTHAVEM64OPT),Y)
                COPTIONS   = -m64
            endif
        endif

        GOTCLANG= $(shell $(_CC) -dM -E - </dev/null 2> /dev/null |grep __clang__|head -1|cut -c19)
        ifeq ($(GOTCLANG),1)
            COPTIONS   += -fPIC
        endif

        GOTFREEBSD= $(shell uname -o 2>&1|awk ' /FreeBSD/ {print "1";exit}')
        ifeq ($(GOTFREEBSD),1)
            DEFINES  +=-DMPICH_NO_ATTR_TYPE_TAGS
#	    LDOPTIONS +=-Wl,-rpath=/usr/local/lib/gcc7
            LDOPTIONS += $(shell mpif90  -show 2>&1 |cut -d " " -f 2) 
            ARFLAGS = rU
        endif

        ifeq ($(_FC),gfortran)
            ifneq ($(DONTHAVEM64OPT),Y)
                FOPTIONS   = -m64
            endif
            ifdef  USE_ASAN
                FOPTIONS += -fsanitize=address -fsanitize-recover=address
                LDOPTIONS += -fsanitize=address -fsanitize-recover=address
            endif
            ifdef  USE_FPE
                ifdef USE_FLANG
                    $(info     )
                    $(info     USE_FPE not ready for flang)
                    $(info     )
                    $(error )
                else
                    FOPTIONS += -ffpe-trap=invalid,zero,overflow  -fbacktrace
                endif
            else
                FOPTIONS   += -ffast-math #-Wunused  
            endif
            ifeq ($(V),-1)
                FOPTIONS += -w
                COPTIONS += -w
            else
                FOPTIMIZE  += -Wuninitialized
                ifeq ($(_CC),$(findstring $(_CC),gcc clang))
                    COPTIONS += -Wall
		endif
                ifeq ($(GNU_GE_4_8),true)
                    ifndef USE_FLANG
                        FOPTIMIZE  += -Wno-maybe-uninitialized
                    endif
                else
                    FOPTIONS   += -Wuninitialized
                endif
            endif

            DEFINES  += -DGFORTRAN
            DEFINES  += -DCHKUNDFLW -DGCC4

            ifeq ($(USE_FLANG),1)
                GNU_GE_4_6=true
                FOPTIONS+=-mcmodel=medium
                FOPTIONS+=-mcmodel=medium -fno-backslash
                COPTIONS+=-mcmodel=medium
                CFLAGS_FORGA = -mcmodel=medium
                FFLAGS_FORGA = -mcmodel=medium
            else
                GNUMAJOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | grep __GNUC__ |cut -c18-)
                ifdef GNUMAJOR
                    GNUMINOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __GNUC_MINOR | cut -c24)
                    GNU_GE_4_6 = $(shell [ $(GNUMAJOR) -gt 4 ] || [ $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 6 ] && echo true)
                    GNU_GE_4_8 = $(shell [ $(GNUMAJOR) -gt 4 ] || [ $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 8 ] && echo true)
                endif
                GNU_GE_6 = $(shell [ $(GNUMAJOR) -ge 6  ] && echo true)
                GNU_GE_8 = $(shell [ $(GNUMAJOR) -ge 8  ] && echo true)
                GNU_GE_10 = $(shell [ $(GNUMAJOR) -ge 10  ] && echo true)
            endif

            ifeq ($(GNU_GE_4_6),true)
                DEFINES  += -DGCC46
            endif

            ifeq ($(GNU_GE_4_8),true)
                ifeq ($(_CPU),ppc64le)
                    FDEBUG =-O0 -g 
#               else
#                   FDEBUG =-O2 -g
                endif
                FDEBUG +=-fno-aggressive-loop-optimizations
                FOPTIONS +=-fno-aggressive-loop-optimizations
                FOPTIMIZE +=-fno-aggressive-loop-optimizations
                FFLAGS_FORGA += -fno-aggressive-loop-optimizations
            endif

            ifeq ($(GNU_GE_8),true)
                FOPTIONS   += -std=legacy
            endif

            ifdef USE_OPENMP
                FOPTIONS  += -fopenmp
                LDOPTIONS += -fopenmp
                ifdef USE_OFFLOAD
                    DEFINES +=-DUSE_F90_ALLOCATABLE -DUSE_OMP_TEAMS_DISTRIBUTE
                endif
            endif
        endif


        ifeq ($(_FC),gfortran)
            ifdef USE_I4FLAGS
#does not exists
#               FOPTIONS += -fdefault-integer-4
            else
                FOPTIONS += -fdefault-integer-8
            endif
        else ifeq ($(_FC),crayftn)
            ifdef USE_I4FLAGS
                FOPTIONS += -s integer32
            else
                FOPTIONS += -s integer64
            endif
        else ifeq ($(_FC),frt)
            ifdef USE_I4FLAGS
                FOPTIONS += -CcdLL8
            else
                FOPTIONS += -CcdLL8 -CcdII8
            endif
        else
            ifdef USE_I4FLAGS
                FOPTIONS += -i4
            else
                FOPTIONS += -i8
            endif
        endif


        DEFINES  += -DEXT_INT
#       MAKEFLAGS = -j 1 --no-print-directory
        ifeq ($(BLASOPT),)
            CORE_SUBDIRS_EXTRA += blas
        endif
        ifeq ($(LAPACK_LIB),)
            CORE_SUBDIRS_EXTRA += lapack
        endif
        RANLIB = echo
        DEFINES   +=   -DLINUX -DLINUX64

        ifeq ($(_CPU),ia64)
#           Itanium  
#           g77 not working 
#           i4 not working 
#
            FC=ifort
            CC=gcc
            DEFINES   +=   -DLINUXIA64 
            COPTIMIZE = -O1
            ifdef USE_SHARED
                FOPTIONS+= -fPIC
            endif

            ifeq ($(FC),ifort)
                _IFCV9= $(shell ifort -v  2>&1|egrep "Version "|head -n 1|awk '/9./ {print "Y"}; /10./ {print "Y"; exit}')
                _IFCV81= $(shell ifort -v  2>&1|egrep "Version "|head -n 1|awk ' /8.1/  {print "Y";exit}; /9./ {print "Y"; exit}; /10./ {print "Y"; exit}')
                _IFCV8= $(shell ifort -v  2>&1|egrep "Version "|head -n 1|awk ' /8./  {print "Y";exit}; /9./ {print "Y"; exit}; /10./ {print "Y"; exit}')

                ifeq ($(_IFCV8),Y)
                    DEFINES+= -DIFCV8
#                   FOPTIONS += -quiet
                endif
                ifeq ($(_IFCV81),Y)
                    DEFINES+= -DIFCV81
                endif
                ITANIUMNO = $(shell   cat /proc/cpuinfo | egrep family | head -n 1  2>&1 | awk ' /Itanium 2/ { print "-tpp2"; exit };/Itanium/ { print "-tpp1"}')
                FOPTIONS   += -auto -w -ftz $(ITANIUMNO)
                ifdef  USE_GPROF
                    FOPTIONS += -qp
                endif
                ifdef  USE_DEBUG
                    FOPTIONS += -g
                endif
                ifeq ($(_IFCV8),Y)
                    FOPTIONS+= -check nobounds -align dcommons -fpe1
                    FOPTIONS+= -warn truncated_source
                else
                    FOPTIONS+= -align
                endif
                DEFINES  +=   -DIFCLINUX
                FVECTORIZE =  -noalign -O3 -pad  -mP2OPT_hlo_level=2
                FOPTIMIZE =  -O3 -pad -mP2OPT_hlo_level=2  #-mP2OPT_align_array_to_cache_line=TRUE 
                ifeq ($(_IFCV81),Y)
#                   FOPTIMIZE+= -IPF_fp_relaxed # breaks nwdft/xc/xc_pw91lda
                endif
                ifeq ($(_IFCV8),Y)
#                   EXTRA_LIBS += -quiet
                    FDEBUG = -g -O2
                else
                    EXTRA_LIBS += -Vaxlib 
                    FDEBUG = -g -O2
                endif
                ifeq ($(_IFCV9),Y)
                    FOPTIONS+=  -fltconsistency
                endif
                ifdef  USE_GPROF
                    EXTRA_LIBS += -qp
                endif
                LDOPTIONS =   -Qoption,link,--relax # -Qoption,link,-Bstatic  
                ifeq ($(BUILDING_PYTHON),python)
                    LDOPTIONS +=  -Qoption,link,--export-dynamic
                    EXTRA_LIBS += -lutil
                endif
                LDOPTIONS += $(FDEBUG)
                LINK.f = ifort $(LDFLAGS)  
            else
                noefc:
	            @echo 
	            @echo Please do not set FC on linux/ia64
	            @echo the makefile will use the Intel compiler
	            @echo 
	            @exit 1
            endif
            ifdef USE_SHARED
                COPTIONS += -fPIC
            endif

#           CORE_LIBS +=  $(BLASOPT) -lnwclapack -lnwcblas
        endif # end of ia64 bit


        ifeq ($(_CPU),x86_64)
#
#edo        MAKEFLAGS = -j 2 --no-print-directory
            COPTIMIZE = -O1
            ifeq ($(NWCHEM_TARGET),CYGWIN64)
                DEFINES += -DCYGWIN -DCYGNUS
            endif

            ifeq ($(NWCHEM_TARGET),CATAMOUNT)
                FC=pgf90
                CC=gcc
            endif

            # support for Intel(R) Fortran compiler
            ifeq ($(_FC),ifxold)
                DEFINES += -DIFCV8 -DIFCLINUX
                FOPTIONS += -fpp -align
                FOPTIMIZE = -g -O3 -fimf-arch-consistency=true
                ifdef USE_I4FLAGS
                else
                    FOPTIONS += -i8
                endif
                ifdef USE_OPENMP
                    FOPTIONS += -fopenmp
                    ifdef USE_OFFLOAD
                        FOPTIONS += -fopenmp-targets=spirv64
                    endif
                endif
                ifdef IFX_DEBUG
                    # debugging remove at some point
                    FOPTIONS += -std95 -what
                endif
                FDEBUG = $(FOPTIMIZE)
            endif


            # support for traditional Intel(R) Fortran compiler
            ifeq ($(_FC),ifort)
                ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                    _GOTSSE3= $(shell cat /proc/cpuinfo | egrep sse3 | tail -n 1 | awk ' /sse3/  {print "Y"}')
                    _GOTSSE42= $(shell cat /proc/cpuinfo | egrep sse4_2 | tail -n 1 | awk ' /sse4_2/  {print "Y"}')
                    _GOTAVX= $(shell cat /proc/cpuinfo | egrep avx | tail -n 1 | awk ' /avx/  {print "Y"}')
                    _GOTAVX2= $(shell cat /proc/cpuinfo | egrep fma | tail -n 1 | awk ' /fma/  {print "Y"}')
                    _GOTAVX512F= $(shell cat /proc/cpuinfo | egrep avx512f | tail -n 1 | awk ' /avx512f/  {print "Y"}')
                endif
                _IFCE = $(shell ifort -V  2>&1 |head -1 |awk ' /64/ {print "Y";exit};')
                _IFCV7= $(shell ifort -v  2>&1|egrep "Version "|head -n 1|awk ' /7./  {print "Y";exit}')
                _IFCV11= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 11) {print "Y";exit}}')
                _IFCV12= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 12) {print "Y";exit}}')
                _IFCV14= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 14) {print "Y";exit}}')
                _IFCV15ORNEWER=$(shell ifort -logo  2>&1|egrep "Version "|head -n 1 | sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 15) {print "Y";exit}}')
                _IFCV17=$(shell ifort -logo  2>&1|egrep "Version "|head -n 1 | sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 17) {print "Y";exit}}')
                _IFCV18=$(shell ifort -logo  2>&1|egrep "Version "|head -n 1 | sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 18) {print "Y";exit}}')

#               Intel EM64T is required
                ifneq ($(_IFCE),Y)
                    defineFCE:
	                @echo
	                @echo "   " ifort missing or not suitable x86_64 CPUs
	                @echo
	                @exit 1
                endif

                ifneq ($(_IFCV7),Y)
#                   to get EM64T
#                   Intel 8.1 is required
                else
                    @echo ifort 8.1 is required for x86_64 CPUs
                    @exit 1
                endif

                FDEBUG= -O2 -g
                FOPTIMIZE = -O3  -unroll

                ifndef USE_IFX
                    FOPTIMIZE += -ip
                endif

                FOPTIONS += -align -fpp

#               might be not need and the root cause for https://github.com/nwchemgit/nwchem/issues/255
#               CPP=fpp -P
#
                ifeq ($(_IFCV15ORNEWER), Y)
#                   fpp seems to get lost with ifort 15 in the offload bit
#                   only use EXPLICITF for offload because otherwise we want debugging to be easy
#                   FOPTIONS +=  -Qoption,fpp,-P -Qoption,fpp,-c_com=no  -allow nofpp_comments 
                    ifdef USE_OPTREPORT
                        FOPTIONS += -qopt-report-file=stderr
                        FOPTIONS += -qopt-report=3 -qopt-report-phase=vec,cg,loop,ipo
                        ifeq ($(_IFCV17), Y)
                            FOPTIONS += -qopt-report-annotate-position=both
                        endif
                        FOPTIONS += -qopt-report-file=stderr
                    endif
                    ifeq ($(V),-1)
                        FOPTIONS += -diag-disable=7713,8291,15009
                    endif
#                   to avoid compiler crashes on simd directive. e.g .Version 15.0.2.164 Build 20150121
                    ifdef USE_NOSIMD
                        FOPTIONS  += -no-simd
                    endif
                    ifdef USE_OPENMP
                        ifdef USE_IFX
                            FOPTIONS += -fiopenmp
                            ifdef USE_OFFLOAD
                                FOPTIONS += -fopenmp-targets=spirv64
                            endif
                        else
                            FOPTIONS += -qopenmp
                            ifdef USE_OPTREPORT
                                FOPTIONS += -qopt-report-phase=openmp
                            endif
                        endif
                    else
                        FOPTIONS += -qno-openmp
                    endif
                else
                    FOPTIONS += -vec-report6
                    ifdef USE_OPENMP
                        FOPTIONS += -openmp
                        FOPTIONS += -openmp-report2
                    endif
                endif

                ifdef USE_VTUNE
                    ifeq ($(VTUNE_AMPLIFIER_XE_DIR),)
                        $(info USE_VTUNE requires VTUNE_AMPLIFIER_XE_DIR to be set)
                        $(error )
                    endif
                    COPTIONS += -I$(VTUNE_AMPLIFIER_XE_DIR)/include
                    FOPTIONS += -DUSE_VTUNE
                    LDOPTIONS += -L$(VTUNE_AMPLIFIER_XE_DIR)/lib64 
                    EXTRA_LIBS += -littnotify
                endif

                DEFINES+= -DIFCV8 -DIFCLINUX

                ifeq ($(FC),ifc)
                    FOPTIONS += -quiet
                endif

                ifdef  USE_FPE
                    FOPTIONS += -fpe0 -traceback #-fp-model  precise
                endif

                ifeq ($(_IFCV11),Y) 
#                   next 2 lines needed for fp accuracy
                    FDEBUG += -fp-model source
                    ifeq ($(_IFCV12),Y) 
                        FOPTIONS += -fimf-arch-consistency=true
                    endif
                    ifdef USE_KNL
                        FOPTIMIZE += -xMIC-AVX512 
#illegal instr?         FOPTIONS += -qopt-assume-safe-padding
                        FOPTIONS += -align array64byte
                        DEFINES+= -DINTEL_64ALIGN
                    else
#                       FOPTIMIZE += -xHost
                        ifndef USE_IFX
#                           crazy simd options
                            ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                                ifeq ($(_IFCV17), Y)
                                    ifeq ($(_GOTAVX512F),Y)
                                        FOPTIMIZE += -axCORE-AVX512
                                    else ifeq ($(_GOTAVX2),Y)
                                        FOPTIMIZE += -axCORE-AVX2
                                    else ifeq ($(_GOTAVX),Y)
                                        FOPTIMIZE += -axAVX
                                    else ifeq ($(_GOTSSE42),Y)
                                        FOPTIMIZE += -axSSE4.2
                                    else ifeq ($(_GOTSSE3),Y) 
                                        FOPTIMIZE += -axSSE3
                                    endif
                                endif
                            endif
                            FOPTIONS += -finline-limit=250
                        endif
                    endif
                else
                    ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                        ifeq ($(_GOTSSE3),Y) 
                             FOPTIMIZE += -xP -no-prec-div
                        else
                             FOPTIMIZE +=  -tpp7
                             FOPTIMIZE += -xW
                        endif
                    endif
                    ifndef USE_IFX
                        FOPTIMIZE +=  -ip
                    endif
                endif


                # configuration options for MEMKIND library on Intel Xeon Phi processors
                ifdef USE_FASTMEM
                    ifdef MEMKIND_PATH
                        FASTMEM_OPTIONS_LD += -L$(MEMKIND_PATH)
                    endif
                    ifdef USE_MEMKIND_TBB
                        FASTMEM_OPTIONS_LD += -L$(TBB_PATH)
                    endif
                    FASTMEM_OPTIONS_LD += -lmemkind
                    ifdef USE_MEMKIND_TBB
                        FASTMEM_OPTIONS_LD += -ltbbmalloc
                    endif
                    EXTRA_LIBS += $(FASTMEM_OPTIONS_LD)
#                   we need to use ALLOCATABLE data for MEMKIND
                    DEFINES += -DUSE_F90_ALLOCATABLE -DUSE_FASTMEM
                endif

            endif # _FC = ifort (i think)
#
            ifeq ($(_FC),pgf90)
                ifeq ($(FC),ftn)
                    LINK.f = ftn  $(LDFLAGS) $(FOPTIONS)
                endif
                ifeq ($(NWCHEM_TARGET),CATAMOUNT)
                    LINK.f = ftn  $(LDFLAGS) $(FOPTIONS)
                endif
            endif

            ifeq ($(FC),pathf90)
#               pathscale 1.3 compiler
#               tested Sep 30 2004 on RH AW3
                FOPTIONS   += -cpp -Wp,-P
                FOPTIONS   += -fno-second-underscore -fixedform
                FOPTIONS   += -align64
#               FOPTIONS   += -LANG:heap_allocation_threshold=0
                FOPTIMIZE   = -O3 -OPT:Ofast:IEEE_arith=1:IEEE_NaN_inf=ON:Olimit=12000:ro=1:fold_reassociate=ON#:div_split=OFF:fast_nint=OFF
                FVECTORIZE  = -O3 -OPT:Ofast:ro=1 -fno-math-errno
                DEFINES  += -DCHKUNDFLW -DPSCALE
                FDEBUG = -g -O1
                LDOPTIONS = -Wl,--warn-once   -Wl,--relax
            endif

            ifeq ($(_CC),pgcc)
#               COPTIONS   =   -O
            endif

            ifeq ($(_CC),icc)
                ICCV15ORNEWER=$(shell icc -V  2>&1|egrep "Version "|head -n 1 | sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 15) {print "Y";exit}}')
                ifdef USE_KNL
                    COPTIONS   +=   -xMIC-AVX512 -ftz
                    DEFINES+= -DINTEL_64ALIGN
                else
#                   COPTIONS   +=   -xHost -ftz
                    COPTIONS   +=   -ftz
                endif

                ifeq ($(ICCV15ORNEWER), Y)
                    ifdef USE_OPTREPORT
                        COPTIONS   += -qopt-report-phase=vec  -qopt-report-file=stderr
                    endif
                    ifdef USE_OPENMP
                        COPTIONS += -qopenmp
                        ifdef USE_OPTREPORT
                            COPTIONS += -qopt-report-phase:openmp
                        endif
                    endif
                else
                    COPTIONS   += -vec-report=1
                    ifdef USE_OPENMP
                        COPTIONS += -openmp
                        COPTIONS += -openmp-report=2
                    endif
                endif
                COPTIMIZE =  -O3
                COPTIMIZE += -ip -no-prec-div
            endif


            ifeq ($(_CC),$(findstring $(_CC),gcc clang))
                COPTIONS   +=   -O3 -funroll-loops -ffast-math 
                ifdef USE_OPENMP
                    COPTIONS += -fopenmp
                endif
            endif


            ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                ifdef USE_GCC34
                    COPTIONS  +=   -march=k8 -mtune=k8
                endif
            endif

#           CORE_LIBS +=  $(BLASOPT) -lnwclapack -lnwcblas
            ifdef  USE_GPROF
                ifeq ($(NWCHEM_TARGET),CATAMOUNT)
                    FOPTIONS   += -Mprof=func#,lines
                    LDOPTIONS   += -Mprof=func#,lines
                    COPTIONS   += -O2 -finstrument-functions
                else
                    FOPTIONS += -pg
                    COPTIONS += -pg
                    LDOPTIONS += -pg
                    LDFLAGS += -pg
                endif
            endif

            ifdef  USE_DEBUG
                FOPTIONS += -g
                COPTIONS += -g
                LDOPTIONS += -g
                LDFLAGS += -g
            endif

            ifeq ($(_FC),gfortran)
#               gcc version 4.1.0 20050525 (experimental)
                ifdef  USE_GPROF
                    FOPTIONS += -pg
                    COPTIONS += -pg
                    LDOPTIONS += -pg
                    LDFLAGS += -pg
                endif

                ifdef  USE_DEBUG
                    FOPTIONS += -g
                    COPTIONS += -g
                    LDOPTIONS += -g
                    LDFLAGS += -g
                endif

                LINK.f = $(FC)  $(LDFLAGS) 
                FOPTIMIZE  += -O3 
                FOPTIMIZE  += -mfpmath=sse # 

                ifeq ($(GNU_GE_6),true)
                    FOPTIMIZE += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
                    FOPTIONS += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
                   FDEBUG += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
                endif

                ifndef USE_FPE
                    FOPTIMIZE  += -ffast-math #2nd time
                endif

                ifdef USE_FLANG
                    ifdef USE_OPTREPORT
                        FOPTIMIZE  += -Rpass=loop-vectorize -Rpass-missed=loop-vectorize -Rpass-analysis=loop-vectorize
                    endif
                else

                FOPTIMIZE  += -fprefetch-loop-arrays #-ftree-loop-linear
            endif

            ifeq ($(GNU_GE_4_8),true)
                FOPTIMIZE  += -ftree-vectorize   
                ifdef USE_OPTREPORT
                    FOPTIMIZE  += -fopt-info-vec
                endif
            endif

            ifdef USE_FLANG
#               AOMP flang crashes with -g in source using block data
                FDEBUG =  -O
            else
                FDEBUG += -g -O0
                ifeq ($(GNU_GE_4_8),true)
                    FDEBUG +=-fno-aggressive-loop-optimizations
                endif
            endif

            ifdef USE_F2C
#               possible segv with use of zdotc (e.g. with GOTO BLAS)
#               http://gcc.gnu.org/bugzilla/show_bug.cgi?id=20178
                FOPTIONS +=  -ff2c -fno-second-underscore
            endif

            ifeq ($(GNU_GE_4_6),true)
                ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                    FOPTIMIZE +=  -mtune=native
                endif
#               causes slowdows in mp2/ccsd
#               FOPTIONS += -finline-functions
            endif
#           FVECTORIZE  += -ftree-vectorize -ftree-vectorizer-verbose=1
#
            ifdef  USE_FPE
                ifdef USE_FLANG
                    $(info     )
                    $(info     USE_FPE not ready for flang)
                    $(info     )
                    $(error )
                else
                    FOPTIONS += -ffpe-trap=invalid,zero,overflow  -fbacktrace
                endif
            endif

            ifeq ($(GOTMINGW64),1)
                 EXTRA_LIBS += -lwsock32
            endif
        else

        endif


        ifeq ($(_FC),crayftn)
#           Jeff: Cray Fortran supports preprocessing as of version 8.2.2 (at least)
#           EXPLICITF = FALSE
            FOPTIONS += -hsystem_alloc -hoverindex
#           workaround for vectorization failures with cce 11
            FOPTIONS += -hfp1
            ifdef BUILD_OPENBLAS
#               avoid replacing code with library calls (eg _dgemm_) to avoid clash with openblas symbols
                FOPTIONS += -hnopattern
            endif
#               USE_POSIXF is required because getlog is provided (GNU extension)
                DEFINES += -DCRAYFORTRAN -DUSE_POSIXF
                ifdef  USE_FPE
                    FOPTIONS   +=  -Ktrap=fp
                endif
                ifdef USE_OPENMP
                    FOPTIONS   +=  -homp
                endif
                FDEBUG   =  -O scalar1,vector1,ipa1  -g
                FOPTIMIZE = -O scalar3,vector2,ipa2
            endif

            ifeq ($(_FC),craycc)
                COPTIONS   =   -O
            endif
        endif


        ifeq ($(_CPU),$(findstring $(_CPU),aarch64))

            ifeq ($(_CC),gcc)
                COPTIONS   +=   -O3 -funroll-loops -ffast-math 
                ifdef USE_OPENMP
                    COPTIONS += -fopenmp
                endif
            endif

            ifeq ($(_CC),armclang)
                COPTIONS += -O3 -funroll-loops
                ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                    ifeq ($(BLAS_SIZE),8)
                        COPTIMIZE +=  -armpl=ilp64
                    else
                        COPTIMIZE +=  -armpl=lp64
                    endif

                    ifdef USE_A64FX
                        COPTIMIZE += -mtune=a64fx -mcpu=a64fx 
                    else
                        COPTIMIZE += -mtune=native -march=native
                    endif 
                endif

                ifdef USE_OPENMP
                    COPTIONS += -fopenmp
                endif
            endif


            ifeq ($(_CC),fcc)
                COPTIONS += -O3
                ifdef USE_OPENMP
                    COPTIONS += -Kopenmp
                endif
            endif

            ifeq ($(_FC),gfortran)
                ifdef  USE_GPROF
                    FOPTIONS += -pg
                    COPTIONS += -pg
                    LDOPTIONS += -pg
                    LDFLAGS += -pg
                endif
                LINK.f = $(FC)  $(LDFLAGS) 
                FOPTIMIZE  += -O3 
                ifeq ($(GNU_GE_6),true)
                    FOPTIMIZE += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
                    FOPTIONS += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
                    FDEBUG += -fno-tree-dominator-opts # solvation/hnd_cosmo_lib breaks
                endif

                FOPTIMIZE  += -fprefetch-loop-arrays #-ftree-loop-linear
                ifeq ($(GNU_GE_4_8),true)
                    FOPTIMIZE  += -ftree-vectorize
                         ifdef USE_OPTREPORT
                              FOPTIMIZE += -fopt-info-vec
			 endif
                endif

                FDEBUG += -g -O

                ifeq ($(GNU_GE_4_6),true)
                    ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
                        ifdef USE_A64FX
                            FOPTIMIZE += -mtune=a64fx -mcpu=a64fx
                            FOPTIMIZE += -march=armv8.2-a+sve
                        else
                            FOPTIMIZE += -mtune=native -march=native
                        endif
                        FOPTIMIZE += -ffp-contract=fast
                             ifdef USE_OPTREPORT
			          FOPTIMIZE += -fopt-info-vec
			     endif
                        FOPTIMIZE += -fstack-arrays
                    endif
#                   causes slowdows in mp2/ccsd
#                   FOPTIONS += -finline-functions
                endif
                ifndef USE_FPE
                    FOPTIMIZE  += -ffast-math #2nd time
                endif
                ifdef  USE_FPE
                    FOPTIONS += -ffpe-trap=invalid,zero,overflow  -fbacktrace
                endif
            endif  # end of gfortran

            # A64fx
            ifeq ($(FC),frt)

                DEFINES += -DFUJITSU
                FOPTIONS += -fs

                LINK.f = $(FC)  $(LDFLAGS)
                FOPTIMIZE  = -O3

                ifeq ($(V),1)
                    $(info     FUJITSU FOPTIMIZE = ${FOPTIMIZE})
                endif

                ifdef USE_OPENMP
                    FOPTIONS  += -Kopenmp
                    LDOPTIONS += -Kopenmp
                endif

                FDEBUG += -g -O

            endif


            ifeq ($(FC),armflang)
                ifdef USE_SHARED
                    FOPTIONS+= -fPIC
                endif

                ifdef USE_OPENMP
                  FOPTIONS  += -fopenmp
                  LDOPTIONS += -fopenmp
                endif

                DEFINES   +=   -DARMFLANG
                LINK.f = $(FC)  $(LDFLAGS) 
                FOPTIMIZE  = -O3 -Mfma -ffp-contract=fast -fno-backslash

                ifeq ($(V),1)
                    $(info     ARMFLANG FOPTIMIZE = ${FOPTIMIZE})
                endif
                ifeq ($(V),-1)
                    FOPTIONS += -w
                endif

                FDEBUG += -g -O
                ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1) 
                    ifeq ($(BLAS_SIZE),8)
                        FOPTIMIZE +=  -armpl=ilp64
                    else
                        FOPTIMIZE +=  -armpl=lp64
                    endif

                    ifdef USE_A64FX
#                       mpcu=a64fx breaks integrals	    
                        FOPTIMIZE += -mtune=a64fx #-mcpu=a64fx 
                    else
                        FOPTIMIZE +=  -mcpu=native
                    endif
                endif

                ifndef USE_FPE
                    FOPTIMIZE  += -ffast-math #2nd time
                endif
                ifdef  USE_FPE
                    FOPTIONS += -ffpe-trap=invalid,zero,overflow  -fbacktrace
                endif
            endif
        endif # end of aarch64


        ifeq ($(_CPU),$(findstring $(_CPU), ppc64 ppc64le))
            # Tested on Red Hat Enterprise Linux AS release 3 (Taroon Update 3)
            # Tested on SLES 9
            # Feb 5th 2005
            # xlf v9.1
            # xlc v7.0 
            # gcc-3.2.3-42 

            #gfortran become default      FC=xlf
            #      CC=gcc
            ifeq ($(CC),xlc)
                COPTIONS  +=  -q64 -qlanglvl=extended
            else
                #this for gcc/cc
                COPTIONS  +=  -m64  -O
            endif

            ifeq ($(_FC),xlf)
                #RSQRT=y breaks intchk QA
                FOPTIONS  =  -q64 -qextname -qfixed #-qnosave  #-qalign=4k
                FOPTIONS +=  -NQ40000 -NT80000 -qmaxmem=8192 -qxlf77=leadzero
                ifeq ($(V),-1)
                   FOPTIONS += -qsuppress=cmpmsg -qsuppress=1501-264 -qsuppress=1500-030
		endif
                ifdef  USE_GPROF
                    FOPTIONS += -pg
                    LDOPTIONS += -pg
                endif
                ifdef  USE_DEBUG
                    FOPTIONS += -g
                    LDOPTIONS += -g
                endif

                FOPTIMIZE= -O3 -qstrict -qcache=auto  
#               qarch and qtune break trobsa.F with xlf 15.1.2
#               FOPTIMIZE+= -qarch=auto -qtune=auto 
                FDEBUG= -O2 -g
#               XLFMAC=y
                DEFINES  +=   -DXLFLINUX -DCHKUNDFLW

                ifdef USE_OPENMP
                    FOPTIONS  += -qsmp=omp
		    LDOPTIONS += -qsmp=omp
                    ifdef USE_OFFLOAD
                        DEFINES +=-DUSE_F90_ALLOCATABLE -DOPENMP_OFFLOAD -DUSE_OMP_TEAMS_DISTRIBUTE
                        OFFLOAD_FOPTIONS = -qtgtarch=sm_70 -qoffload
                        LDOPTIONS += -qoffload -lcudart -L$(NWC_CUDAPATH)
                    endif
                    LINK.f   = xlf_r   $(LDFLAGS)
                endif

                ifdef USE_I4FLAGS
                    FOPTIONS += -qintsize=4
                else
                    FOPTIONS += -qintsize=8
                endif
                ifeq ($(V),-1)
                    FOPTIONS += -w
                endif
            endif

            ifdef USE_ESSL
                CORE_SUBDIRS_EXTRA = lapack
                CORE_LIBS +=  -lnwclapack
            endif
#           CORE_LIBS +=  $(BLASOPT) -lnwclapack -lnwcblas
#           EXTRA_LIBS +=  -dynamic-linker /lib64/ld64.so.1 -melf64ppc -lxlf90_r -lxlopt -lxlomp_ser -lxl -lxlfmath -ldl -lm -lc -lgcc -lm
        endif # end of ppc64 arch

	ifeq ($(_CC),pgcc)
            COPTIONS += -fast -O3 -Munroll
            ifdef USE_OPENMP
                COPTIONS += -mp
            endif
	endif


        ifeq ($(_FC),pgf90)
            FOPTIONS   += -Mdalign -Kieee
            ifeq ($(_CPU),x86_64)
                FOPTIONS   +=  -Mllalign
            endif
            FOPTIONS   += -Mbackslash
            FOPTIONS   += -Mcache_align  # -Kieee 
            FOPTIMIZE   =  -fast -O3 -Mvect=simd  -Munroll # -Mipa=fast
            ifdef USE_OPTREPORT
                 FOPTIMIZE += -Minfo=loop
            endif
            FVECTORIZE   = -fast    -O3   #-Mipa=fast
            FDEBUG = -g -O2
            DEFINES  += -DCHKUNDFLW -DPGLINUX
            ifeq ($(shell $(CNFDIR)/check_env.sh $(USE_HWOPT)),1)
              FOPTIONS += -tp host
            else
              ifeq ($(_CPU),x86_64)
                FOPTIONS += -tp px
	      endif
            endif
            ifdef USE_OPENMP
	      ifndef UNSET_OPENMP
                FOPTIONS  += -mp
                ifdef USE_OPTREPORT
                    FOPTIONS  += -Minfo=mp
                endif
                LDOPTIONS += -mp
	      endif
            endif
	    ifdef USE_FPE
	        FOPTIONS += -traceback
		FOPTIONS += -Ktrap=inv,divz,ovf
	    endif
        endif


        ifeq ($(NWCHEM_TARGET),CATAMOUNT)
            DEFINES  += -DCATAMOUNT
        endif


        # Jeff: FreeBSD does not link libm automatically with flang
        ifeq ($(USE_FLANG),1)
            EXTRA_LIBS += -lm
        endif

    endif
endif
#endof of LINUX64



ifeq ($(TARGET),$(findstring $(TARGET),BGL BGP BGQ))
#
    CORE_SUBDIRS_EXTRA = lapack blas
    ARFLAGS = urs
    INSTALL = @echo $@ is built
    CPP=/usr/bin/cpp  -P -C -traditional
    LDOPTIONS =  -Wl,--relax

    DEFINES   = -DLINUX
    FDEBUG    = -g
    COPTIMIZE = -g
    FOPTIMIZE = -g
    FOPTIONS  = -g
    DEFINES += -DCHKUNDFLW

#for BGL
    ifeq ($(TARGET),BGL)
        EXPLICITF = TRUE
        DEFINES   += -DBGML -DBGL
        FC         = blrts_xlf
        CC         = $(BGCOMPILERS)/powerpc-bgl-blrts-gnu-gcc
        AR         = $(BGCOMPILERS)/powerpc-bgl-blrts-gnu-ar
        AS         = $(BGCOMPILERS)/powerpc-bgl-blrts-gnu-as
        RANLIB     = $(BGCOMPILERS)/powerpc-bgl-blrts-gnu-ranlib
        DEFINES   += -DXLFLINUX -DBGL
        FOPTIMIZE += -qarch=440 -qtune=440 -qfloat=rsqrt:fltint
        FOPTIONS   = -qEXTNAME -qxlf77=leadzero -NQ40000 -NT80000 -NS2048 -qmaxmem=8192 -qsuppress=1500-030
    endif

#for BGP
    ifeq ($(TARGET),BGP)
        EXPLICITF = TRUE
        DEFINES += -DDCMF -DBGP
        FC     = mpixlf77_r
        CC     = mpicc
        AR     = powerpc-bgp-linux-ar
        AS     = powerpc-bgp-linux-as
        RANLIB = powerpc-bgp-linux-ranlib
        DEFINES += -DXLFLINUX
        FOPTIONS = -qEXTNAME -qxlf77=leadzero -NQ40000 -NT80000 -NS2048 -qmaxmem=8192 -qsuppress=1500-030
        FOPTIONS += -O3 -qstrict -qthreaded -qnosave -qalign=4k
        FOPTIMIZE += -O3 -qarch=450d -qtune=450 -qcache=auto -qunroll=auto -qfloat=rsqrt:fltint
        XLF11 = $(shell bgxlf -qversion  2>&1|grep Version|head -1| awk ' / 11./ {print "Y"}')
    endif

#for BGQ
    ifeq ($(TARGET),BGQ)
        DEFINES += -DPAMI -DBGQ 
        DEFINES += -DEAFHACK -DNOFSCHECK
        AR       = powerpc64-bgq-linux-ar
        AS       = powerpc64-bgq-linux-as
        RANLIB   = powerpc64-bgq-linux-ranlib

        #FC = mpif77
        FC = mpixlf77_r
        _FC = xlf

        ifeq ($(FC),mpif77)
            CC         = mpicc
            DEFINES   += -DGFORTRAN -DGCC4
            FOPTIONS  += -g -funderscoring
            ifeq ($(V),-1)
                FOPTIONS += -w
            else
                FOPTIONS  += -Wuninitialized -Wno-maybe-uninitialized
            endif
            FOPTIMIZE += -O3 -ffast-math
            FDEBUG    += -O1 -g

            # EXT_INT means 64-bit integers are used
            DEFINES   += -DEXT_INT 
            ifndef USE_I4FLAGS
                FOPTIONS  += -fdefault-integer-8
            endif

            # linking ESSL is painful with gfortran
            CORE_LIBS +=  -lnwclapack $(BLASOPT) -lnwcblas

            # Here is an example for ALCF:
            # IBMCMP_ROOT=${IBM_MAIN_DIR}
            # BLASOPT=/soft/libraries/alcf/current/xl/BLAS/lib
            # LAPACK_LIB=/soft/libraries/alcf/current/xl/LAPACK/lib
            # ESSL_LIB=/soft/libraries/essl/current/essl/5.1/lib64
            # XLF_LIB=${IBMCMP_ROOT}/xlf/bg/14.1/bglib64
            # XLSMP_LIB=${IBMCMP_ROOT}/xlsmp/bg/3.1/bglib64
            # XLMASS_LIB=${IBMCMP_ROOT}/xlmass/bg/7.3/bglib64
            # MATH_LIBS="-L${XLMASS_LIB} -lmass -L${LAPACK_LIB} -lnwclapack \
                     -L${ESSL_LIB} -lesslsmpbg -L${XLF_LIB} -lxlf90_r \
                     -L${XLSMP_LIB} -lxlsmp -lxlopt -lxlfmath -lxl \
                     -Wl,--allow-multiple-definition"
            # Note that ESSL _requires_ USE_64TO32 on Blue Gene
        endif

        ifeq ($(FC),mpixlf77_r)
            EXPLICITF  = TRUE
            _FC = xlf
            CC         = mpixlc_r
            DEFINES   += -DXLFLINUX
            XLF11      = $(shell bgxlf -qversion  2>&1|grep Version|head -1| awk ' / 11./ {print "Y"}')
            XLF14      = $(shell bgxlf -qversion  2>&1|grep Version|head -1| awk ' / 14./ {print "Y"}')

            # EXT_INT means 64-bit integers are used
            DEFINES   += -DEXT_INT 

            ifdef USE_I4FLAGS
                FOPTIONS = -qintsize=4
                ifeq ($(BLAS_SIZE),8)
                    @echo "You cannot use BLAS with 64b integers when"
                    @echo "the compiler generates 32b integers (USE_I4FLAGS)!"
                    @exit 1
                endif # BLAS_SIZE
            else
                FOPTIONS = -qintsize=8 
                ifeq ($(BLAS_SIZE),4)
                    ifneq ($(USE_64TO32),y)
                        @echo "You cannot use BLAS with 32b integers when"
                        @echo "the compiler generates 64b integers unless"
                        @echo "you do the 64-to-32 conversion!"
                        @exit 1
                    endif # USE_64TO32
                endif # BLAS_SIZE
            endif # USE_I4FLAGS

            FDEBUG     = -g -qstrict -O3
            FOPTIONS  += -g -qEXTNAME -qxlf77=leadzero
            FOPTIONS  += -qthreaded -qnosave # -qstrict
#           FOPTIMIZE += -g -O3 -qarch=qp -qtune=qp -qcache=auto -qunroll=auto -qfloat=rsqrt
            FOPTIMIZE += -O3 -qarch=qp -qtune=qp -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes #-qnoipa
            FOPTIMIZE += -qreport -qsource -qlistopt -qlist # verbose compiler output

            # ESSL dependencies should be provided by XLF linker
            CORE_LIBS +=  -lnwclapack $(BLASOPT) -lnwcblas
        endif

    endif # end BGQ


endif


###################################################################
#  Some generic settings about programming models, e.g., OpenMP,  #
#  that are orthogonal to the actual compiler used.               #
###################################################################
ifdef USE_OPENMP
    DEFINES += -DUSE_OPENMP
    ifdef USE_OPENMP_TASKS
      DEFINES += -DUSE_OPENMP_TASKS
    endif
    ifdef USE_OFFLOAD
      DEFINES += -DUSE_OFFLOAD
    endif
endif



###################################################################
#  All machine dependent sections should be above here, otherwise #
#  some of the definitions below will be 'lost'                   #
###################################################################
ifeq ($(BUILDING_PYTHON),python)
    ifdef PYTHONVERSION
        GOTPYTHON := $(shell command -v python$(PYTHONVERSION) 2> /dev/null)
        ifndef GOTPYTHON	      
            errorpythonXY:
	        $(info )
	        $(info python$(PYTHONVERSION) not  found in your PATH)
	        $(info )
	        $(error )
        endif
        # check presence of python?-config  
        GOT_PYTHONCONFIG := $(shell command -v python$(PYTHONVERSION)-config 2> /dev/null)
    else	      
        #try python3 first, then python
        GOTPYTHON3 := $(shell command -v python3 2> /dev/null)
        GOTPYTHON2 := $(shell command -v python2 2> /dev/null)
        GOTPYTHON  := $(shell command -v python 2> /dev/null)
        ifdef GOTPYTHON3
            PYTHONVERSION=$(shell python3 -c 'import sys; print("{}.{}".format(sys.version_info.major, sys.version_info.minor))')
        else ifdef GOTPYTHON2
            PYTHONVERSION=$(shell python2 -c 'import sys; print("{}.{}".format(sys.version_info.major, sys.version_info.minor))')
        else ifdef GOTPYTHON
            #last try at python2
            PYTHONVERSION=$(shell python -c 'import sys; print("{}.{}".format(sys.version_info.major, sys.version_info.minor))')
        else
            errorpython3:
	        $(info )
	        $(info Neither python2 nor python3 found in your PATH)$
	        $(info )
	        $(error )
        endif

        # check presence of python?-config  
        GOT_PYTHONCONFIG := $(shell command -v python$(PYTHONVERSION)-config 2> /dev/null)
        ifndef GOT_PYTHONCONFIG	  
            ifdef GOTPYTHON3
                GOT_PYTHONCONFIG := $(shell command -v python3-config 2> /dev/null)
            else ifdef GOTPYTHON2
                GOT_PYTHONCONFIG := $(shell command -v python2-config 2> /dev/null)
            else ifdef GOTPYTHON
                GOT_PYTHONCONFIG := $(shell command -v python-config 2> /dev/null)
            endif
        endif
    endif

    ifndef GOT_PYTHONCONFIG	  
        PYMAJOR:=$(word 1, $(subst ., ,$(PYTHONVERSION)))
        errorpythonconfig:$
	    $(info )
	    $(info python-config not found in your PATH)
	    $(info Please install the packages)
	    $(info python$(PYMAJOR)-dev  (Ubuntu/Debian) or)
	    $(info python$(PYMAJOR)-devel (Redhat/Fedora/Centos))
	    $(info )
	    $(error )
    endif

    #ifdef USE_PYTHONCONFIG
    PYMAJOR:=$(word 1, $(subst ., ,$(PYTHONVERSION)))
    PYMINOR:=$(word 2, $(subst ., ,$(PYTHONVERSION)))
    PYGE38:=$(shell [ $(PYMAJOR) -ge 3 -a $(PYMINOR) -ge 8 ] && echo true)
    ifeq ($(PYGE38),true)
        PYCFG= python$(PYTHONVERSION)-config --ldflags --embed
        ifeq ($(shell uname -s),Darwin)
            PYCFG += | sed -e "s/-lintl //"
        endif
    else
        PYCFG = python$(PYTHONVERSION)-config --ldflags
    endif
    EXTRA_LIBS += -lnwcutil $(shell $(PYCFG))
else
    ifndef PYTHONLIBTYPE
        PYTHONLIBTYPE=a
    endif
    ifndef PYTHONCONFIGDIR
        PYTHONCONFIGDIR=config
    endif
#endif #USE_PYTHONCONFIG
endif
#
######
#PAPI#
######
   
ifdef USE_PAPI 
    DEFINES += -DUSE_PAPI
    ifndef PAPI_LIB
        papierror1:
	    @echo You must define PAPI_LIB in your environment to be the path
	    @echo of your PAPI library installation ... something like
	    @echo     setenv PAPI_LIB ~/papi/lib
	    @exit 1
    endif

    ifndef PAPI_INCLUDE
        papierror2:
	    @echo You must define PAPI_INCLUDE in your environment to be the path
	    @echo of your PAPI headers installation ... something like
	    @echo     setenv PAPI_INCLUDE ~/papi/include
	    @exit 1
    endif

    ifndef LIBPAPI 
        LIBPAPI = -lpapi 
    endif 
    CORE_LIBS += -L$(PAPI_LIB) $(LIBPAPI) 
    INCPATH += -I$(PAPI_INCLUDE)
endif


ifdef USE_YFLOP
    ifndef USE_64TO32
        yfloperr:
	    @echo You must define USE_64TO32 in your environment 
	    @echo to use the USE_YFLOP env. variable
	    @echo 
	    @exit 1
    endif
    DEFINES += -DUSE_YFLOP
endif


#
ifdef USE_FDIST
    DEFINES += -DFDIST
endif

#_USE_SCALAPACK = $(shell cat ${NWCHEM_TOP}/src/tools/build/config.h | awk ' /HAVE_SCALAPACK 1/ {print "Y"}')
#case guard against case when tools have not been compiled yet
ifeq ("$(wildcard ${GA_PATH}/bin/ga-config)","")
else
    ifneq ("$(wildcard ${NWCHEM_TOP}/src/ga_use_scalapack.txt)","")
        _USE_SCALAPACK= $(shell cat $(NWCHEM_TOP)/src/ga_use_scalapack.txt)
    endif
    ifeq ($(_USE_SCALAPACK),)
        _USE_SCALAPACK = $(shell ${GA_PATH}/bin/ga-config  --use_scalapack| awk ' /1/ {print "Y"}')
    endif
endif

ifeq ($(_USE_SCALAPACK),Y)
    DEFINES += -DSCALAPACK
    ifeq ($(XLFBREN),y) 
        CORE_LIBS +=  -brename:.iceil_,.iceil \
                    -brename:.blacs_pinfo_,.blacs_pinfo \
                    -brename:.blacs_get_,.blacs_get \
                    -brename:.blacs_gridinit_,.blacs_gridinit \
                    -brename:.blacs_gridinfo_,.blacs_gridinfo \
                    -brename:.blacs_gridexit_,.blacs_gridexit \
                    -brename:.indxg2p_,.indxg2p \
                    -brename:.descinit_,.descinit \
                    -brename:.numroc_,.numroc \
                    -brename:.pdlamch_,.pdlamch \
                    -brename:.pdsyev_,.pdsyev \
                    -brename:.pdsyevd_,.pdsyevd \
                    -brename:.pdsyevx_,.pdsyevx \
                    -brename:.pdsygvx_,.pdsygvx \
                    -brename:.pdpotri_,.pdpotri \
                    -brename:.pdpotrf_,.pdpotrf \
                    -brename:.pdpotrs_,.pdpotrs \
                    -brename:.pdgetrf_,.pdgetrf \
                    -brename:.pdgetrs_,.pdgetrs 
    endif
endif
CORE_LIBS += $(ELPA) $(SCALAPACK) $(SCALAPACK_LIB)


ifdef USE_64TO32
    CORE_LIBS +=  -l64to32
    NWSUBDIRS += 64to32blas
endif


ifeq ($(LAPACK_LIB),)
    CORE_LIBS +=  -lnwclapack 
else
    CORE_LIBS += $(LAPACK_LIB)
endif


ifeq ($(BLASOPT),)
    CORE_LIBS +=  -lnwcblas 
else
    CORE_LIBS += $(BLASOPT)
endif

ifdef NWCHEM_LINK_CUDA
    CORE_LIBS += -acc -gpu=managed -cuda -cudalib=cublas
endif


ifdef BLASOPT
    BLAS_SUPPLIED=Y
endif

ifdef BLAS_LIB
    BLAS_SUPPLIED=Y
endif

ifndef BLAS_SUPPLIED
    ifndef USE_INTERNALBLAS
        errordgemm:
	    $(info     )
	    $(info NWChem's Performance is degraded by not setting BLASOPT)
	    $(info Please consider using ATLAS, GotoBLAS2, OpenBLAS, Intel MKL,)
	    $(info IBM ESSL, AMD ACML, etc. to improve performance.)
	    $(info If you decide to not use a fast implementation of BLAS/LAPACK,)
	    $(info please define USE_INTERNALBLAS=y and the internal Netlib will be used.)
    endif
else
    ifndef LAPACK_LIB
        ifndef USE_ESSL
            errorlap1:
	        $(error Please define LAPACK_LIB if you have defined BLASOPT or BLAS_LIB)
        endif
    endif
endif


ifdef USE_NOIO
    DEFINES += -DNOIO -DEAFHACK
endif


ifdef USE_SUBGROUPS
    DEFINES += -DGANXTVAL -DUSE_SUBGROUPS
    #turn off peigs for now
else
    ifneq ($(GOTMINGW64),1)
        DEFINES += -DPARALLEL_DIAG
    endif
endif


###################################################################
#  All machine dependent sections should be above here, otherwise #
#  some of the definitions below will be 'lost'                   #
###################################################################
#the new GA uses ARMCI library
ifdef OLD_GA
    ORE_LIBS += -larmci
else
    ifeq ($(ARMCI_NETWORK),ARMCI)
        ifdef EXTERNAL_ARMCI_PATH
            CORE_LIBS += -L$(EXTERNAL_ARMCI_PATH)/lib -larmci
        else
            CORE_LIBS += -larmci
        endif
    else
        CORE_LIBS +=
    endif
endif


# MPI version requires tcgmsg-mpi library
ifdef USE_MPI 
    #ifeq ($(FC),$(findstring $(FC),mpifrt mpfort mpif77 mpxlf mpif90 ftn scorep-ftn))
    ifeq ($(FC),$(findstring $(FC), ftn scorep-ftn))
        LIBMPI =
        MPI_INCLUDE =
        MPI_LIB =
    else
        ifndef MPI_INCLUDE
            # check if mpif90 is present
            MPIF90YN = $(shell $(NWCHEM_TOP)/src/tools/guess-mpidefs --mpi_include)
            ifeq ($(MPIF90YN),mpif90notfound)
                errormpif90:
	            $(info )
	            $(info mpif90 not found. Please add its location to PATH)
	            $(info e.g. export PATH=/usr/local/bin:/usr/lib64/openmpi/bin:...)
	            $(info )
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
        CORE_LIBS += $(patsubst -L-L%,-L%,-L$(MPI_LIB))
    endif 
    ifdef OLD_GA
        CORE_LIBS += -ltcgmsg-mpi $(LIBMPI) 
    else
        CORE_LIBS += $(LIBMPI) 
    endif
else 
    ifdef OLD_GA
      CORE_LIBS += -ltcgmsg 
    else
      CORE_LIBS += 
    endif
endif 


# FFTW3 library inclusion
ifdef USE_FFTW3
    ifndef LIBFFTW3
        LIBFFTW3 = -lfftw3
    endif
    ifdef FFTW3_LIB
        CORE_LIBS += -L$(FFTW3_LIB)
    endif
    CORE_LIBS += $(LIBFFTW3) 
endif


# FEFF library inclusion
ifdef USE_FEFF
    ifndef LIBFEFF
        LIBFEFF = -lfeff
    endif
    ifdef FEFF_LIB
        CORE_LIBS += -L$(FEFF_LIB)
    endif
    CORE_LIBS += $(LIBFEFF) 
endif



# slurm libraries for remaining wall time under slurm resource manager

ifdef SLURM
    EXTRA_LIBS += $(SLURMOPT)
endif

# we build the libxc library
ifdef USE_LIBXC
    DEFINES += -DUSE_LIBXC
    EXTRA_LIBS += -L$(NWCHEM_TOP)/src/libext/libxc/install/lib
    EXTRA_LIBS += -lxcf03 -lxc
endif

# we use an external libxc library out of LIBXC_DIR
ifdef LIBXC_DIR
    DEFINES += -DUSE_LIBXC
    EXTRA_LIBS += -L$(LIBXC_DIR)/lib
    EXTRA_LIBS += -lxcf03 -lxc
endif

ifdef USE_SIMINT
    ifndef SIMINT_HOME
        SIMINT_HOME=$(NWCHEM_TOP)/src/NWints/simint/libsimint_source/simint_install
    endif
    ifdef SIMINT_LIB64
        EXTRA_LIBS += -L$(SIMINT_HOME)/lib64 -lnwc_simint
    else
        EXTRA_LIBS += -L$(SIMINT_HOME)/lib -lnwc_simint
    endif
endif


ifdef USE_PLUMED
    DEFINES += -DUSE_PLUMED
    #check presence of plumed command. TODO
    GOTPLUMED  := $(shell command -v plumed 2> /dev/null)
    ifndef GOTPLUMED
        errorplumed0:
	    $(info )
	    $(info  PLUMED installation not found.)
	    $(info  Please add to your PATH the directory where the plumed command is found )
	    $(info )
    endif
    PLUMED_HOME = $(shell plumed info --configuration|egrep prefix=|head -1|cut -c 8-)
    PLUMED_DYNAMIC_LIBS = $(shell plumed info --configuration|egrep DYNAMIC_LIBS| cut -c 14-)
    PLUMED_HASMPI = $(plumed info --configuration|grep program_can_run_mpi|cut -c 21-21)
    ifeq ($(PLUMED_HASMPI),y)
        DEFINES += -DPLUMED_HASMPI
    endif
    #PLUMED_LOAD= /home/edo/tahoma/apps/plumed262.intel20u2/lib/libplumed.a -ldl  -lstdc++ -lfftw3 -lz -ldl -llapack -lblas   -rdynamic -Wl,-Bsymbolic -fopenmp 
    ifdef PLUMED_DYNAMIC_LIBS
        EXTRA_LIBS += -L$(PLUMED_HOME)/lib -lplumed $(PLUMED_DYNAMIC_LIBS)
    else
        errorplumed:
	    $(info )
	    $(info  PLUMED info command not returning the expected output)
	    $(info  Please file an issue at https://github.com/nwchemgit/nwchem/issues )
	    $(info )
    endif
endif

# CUDA
ifndef CUDA
    CUDA = nvcc
endif

ifdef TCE_CUDA
    DEFINES += -DTCE_CUDA
    CORE_LIBS += $(CUDA_LIBS)
    EXTRA_LIBS += -lstdc++
    ifdef USE_TTLG
        EXTRA_LIBS += -lcublas
    endif
    ifeq ($(_CC),pgcc)
        COPTIONS += -acc
    endif
endif

ifndef HIP
    HIP = hipcc
endif

ifdef TCE_HIP
    DEFINES += -DTCE_HIP $(shell hipconfig --cpp_config)
    CORE_LIBS += $(HIP_LIBS)
    EXTRA_LIBS += -lstdc++
endif

ifdef USE_F90_ALLOCATABLE
    DEFINES += -DUSE_F90_ALLOCATABLE
endif

ifdef GWCMPLX
  ifdef GWEN
    errorgw:
$(info  GWCMPLX and GWEN are incompatible )
$(error )
  endif
  DEFINES += -DGWCMPLX
endif

ifdef GWEN
  DEFINES += -DGWEN
endif

ifdef GWDEBUG
  DEFINES += -DGWDEBUG
endif

# lower level libs used by communication libraries 
#case guard against case when tools have not been compiled yet
#  ifeq ("$(wildcard ${GA_PATH}/bin/ga-config)","")
#  else
COMM_LIBS=  $(shell ${GA_PATH}/bin/ga-config --network_ldflags)
COMM_LIBS +=  $(shell ${GA_PATH}/bin/ga-config --network_libs)
#comex bit
#COMM_LIBS +=  $(shell [ -e ${NWCHEM_TOP}/src/tools/build/comex/config.h ] && grep LIBS\ = ${NWCHEM_TOP}/src/tools/build/comex/Makefile|grep -v _LIBS| cut -b 8-) -lpthread
COMM_LIBS += $(shell [ -e ${GA_PATH}/bin/comex-config ] && ${GA_PATH}/bin/comex-config --libs) -lpthread
ifdef COMM_LIBS 
    CORE_LIBS += $(COMM_LIBS) 
endif 
#endif
ifdef USE_CRAYSHASTA
    CORE_LIBS += -lpmi2
endif
ifdef USE_LINUXAIO
    CORE_LIBS += -lrt
endif

# g++ GNU compatibility (might go away)
#CORE_LIBS += -lstdc++

EXTRA_LIBS += $(CONFIG_LIBS)
CORE_LIBS += $(EXTRA_LIBS)


ifdef OPTIMIZE
    FFLAGS = $(FOPTIONS) $(FOPTIMIZE)
    CFLAGS =  $(COPTIONS) $(COPTIMIZE)
    FFLAGS += $(EXTRA_FOPTIONS) $(EXTRA_FOPTIMIZE)
    CFLAGS += $(EXTRA_COPTIONS) $(EXTRA_COPTIMIZE)
else
# Need FDEBUG after FOPTIONS on SOLARIS to correctly override optimization
    FFLAGS = $(FOPTIONS) $(FDEBUG) 
    CFLAGS = $(COPTIONS) $(CDEBUG) 
    FFLAGS += $(EXTRA_FOPTIONS) $(EXTRA_FDEBUG)
    CFLAGS += $(EXTRA_COPTIONS) $(EXTRA_CDEBUG)
endif

INCLUDES = -I. $(LIB_INCLUDES) -I$(INCDIR) $(INCPATH)
CPPFLAGS = $(INCLUDES) $(DEFINES) $(LIB_DEFINES)
LDFLAGS = $(LDOPTIONS) -L$(LIBDIR) $(LIBPATH)
LIBS = $(NW_MODULE_LIBS) $(CORE_LIBS) 

# I think this will work everywhere, but it might have to become
# machine-dependent 

MKDIR = mkdir
#extract defines to be used with linear algebra libraries
ifdef USE_INTERNALBLAS
    DEFINES += -DBLAS_NOTHREADS
endif
ifdef BUILD_OPENBLAS
    DEFINES += -DOPENBLAS
endif
ifeq ($(shell echo $(BLASOPT) |awk '/openblas/ {print "Y"; exit}'),Y)
    DEFINES += -DOPENBLAS
endif
# NVHPC compilers are distributed wtih OpenBLAS named as libblas/liblapack
ifeq ($(shell echo $(BLASOPT) |awk '/\/nvidia\/hpc_sdk\// {print "Y"; exit}'),Y)
    DEFINES += -DOPENBLAS
endif
ifeq ($(shell echo $(BLASOPT) |awk '/mkl/ {print "Y"; exit}'),Y)
    DEFINES += -DMKL
endif
ifeq ($(shell echo $(BLASOPT) |awk '/blis/ {print "Y"; exit}'),Y)
    DEFINES += -DBLIS
endif
ifeq ($(shell echo $(BLASOPT) |awk '/flexiblas/ {print "Y"; exit}'),Y)
    DEFINES += -DFLEXIBLAS
endif
ifeq ($(shell echo $(BLASOPT) |awk '/Accelerate/ {print "Y"; exit}'),Y)
    DEFINES += -DACCELERATE
endif
ifeq ($(shell echo $(BLASOPT) |awk '/lsci/ {print "Y"; exit}'),Y)
    DEFINES += -DCRAYBLAS
endif
ifeq ($(shell echo $(BLASOPT) |awk '/larmpl/ {print "Y"; exit}'),Y)
    DEFINES += -DARMPL
endif
ifeq ($(shell echo $(BLASOPT) |awk '/latlas/ {print "Y"; exit}'),Y)
    DEFINES += -DBLAS_NOTHREADS
endif
ifeq ($(shell echo $(BLASOPT) |awk '/[Ss][Ss][Ll]2[Bb][Ll][Aa][Mm][Pp]/ {print "Y"; exit}'),Y)
    DEFINES += -DBLAS_OPENMP
else ifeq ($(shell echo $(BLASOPT) |awk '/[Ss][Ss][Ll]2/ {print "Y"; exit}'),Y)
    DEFINES += -DBLAS_NOTHREADS
endif

ifeq ($(shell echo $(BLASOPT) |awk '/lessl/ {print "Y"; exit}'),Y)
    ifeq ($(shell echo $(BLASOPT) |awk '/smp/ {print "Y"; exit}'),Y)
        erroresslsmp:
	    $(info     )
	    $(info essl smp threaded libraries are deprecated)
	    $(info since they conflict with OpenMP parallelization)
	    $(info please use -lessl6464 or -lessl)
	    $(error )
    endif
#   DEFINES += -DBLAS_OPENMP
    DEFINES += -DBLAS_NOTHREADS
#   essl does not has the full lapack library
    EXTRA_LIBS += -lnwclapack
    CORE_SUBDIRS_EXTRA = lapack
endif


#
# Define known suffixes mostly so that .p files don\'t cause pc to be invoked
#
V = 0
ACTUAL_FC := $(FC)
NWFC_-1 = @echo "Compiling $<..."; $(ACTUAL_FC)
NWFC_0 = @echo "Compiling $<..."; $(ACTUAL_FC)
NWFC_1 = $(ACTUAL_FC)
NWFC = $(NWFC_$(V))
ACTUAL_CC := $(CC)
NWCC_-1 = @echo "Compiling $<..."; $(ACTUAL_CC)
NWCC_0 = @echo "Compiling $<..."; $(ACTUAL_CC)
NWCC_1 = $(ACTUAL_CC)
NWCC = $(NWCC_$(V))

.SUFFIXES:
.SUFFIXES: .o .s .F .f .c .cpp

ifndef FLINT
    ifdef EXPLICITF
#
# Needed on machines where FC does not preprocess .F files
# with CPP to get .f files
#
# These rules apply to make-ing of files in specfic directories
.SUFFIXES:
.SUFFIXES: .o .s .F .f .c

.F.o:
	@echo Converting $*.F '->' $*.f
	@$(FCONVERT)
	$(NWFC) -c $(FFLAGS) $*.f
	@$(RM) $*.f

.F.f:
	@echo Converting $*.F '->' $*.f
	@$(FCONVERT)
.f.o:
	$(FC) -c $(FFLAGS) $<

     endif
# 
# More explicit rules to avoid infinite recursion, to get dependencies, and
# for efficiency.  CRAY does not like -o with -c.
#
# These rules apply to make-ing of files in with respect to library files
# both these rules and the rules above are needed.
(%.o): %.F
    ifdef EXPLICITF
	@echo Converting $< '->' $*.f
	@$(FCONVERT)
	$(FC) -c $(FFLAGS) $*.f
        ifndef NWCHEM_KEEPF
	@/bin/rm -f $*.f
        endif
    else
        ifeq ($(XLFMAC),y)
	$(FC)  -c $(FFLAGS) $(INCLUDES) -WF,"$(DEFINES)" $(shell echo $(LIB_DEFINES) | sed -e "s/-D/-WF,-D/g" | sed -e 's/\"/\\"/g')  $<
        else
	$(NWFC)  -c $(FFLAGS) $(CPPFLAGS)  $<
        endif
    endif

(%.o): %.f
	$(NWFC) -c $(FFLAGS) $<

(%.o): %.c
	$(NWCC) -c $(CPPFLAGS) $(CFLAGS) -o $% $<

    ifdef GPU_ARCH
        CUDA_ARCH =  -arch=$(GPU_ARCH) 
    else
        CUDA_ARCH =  -arch=sm_35
    endif


    ifdef TCE_CUDA
        ifdef USE_TTLG
            CUDA_VERS_GE8=$(shell nvcc --version|egrep rel|  awk '/release 9/ {print "Y";exit}; /release 8/ {print "Y";exit};{print "N"}')
            ifeq ($(CUDA_VERS_GE8),N)
                CUDA_FLAGS = -O3 -Xcompiler -std=c++11 -DNOHTIME -Xptxas --warn-on-spills $(CUDA_ARCH) 
            else
                CUDA_FLAGS = -O3  -std=c++11 -DNOHTIME -Xptxas --warn-on-spills $(CUDA_ARCH) 
            endif
(%.o):  %.cu
	$(CUDA) -c -DTCE_CUDA $(CUDA_FLAGS) $(CUDA_INCLUDE) -I$(NWCHEM_TOP)/src/tce/ttlg/includes -o $% $<

        else
(%.o):  %.cu
	$(CUDA) -c -DTCE_CUDA $(CUDA_FLAGS) $(CUDA_INCLUDE) -o $% $<

        endif
    endif

    ifdef TCE_HIP
(%.o):  %.hip.cpp
	$(HIP) -c -DTCE_HIP -fno-gpu-rdc -o $% $<

    endif

(%.o):  %.o

# Preceding line has a tab to make an empty rule

# a .F.f rule is needed for any system target where the default .F.f rule does not work
# AND the EXPLICITF is not already true.  Right now this is only LINUX with g77
    ifneq ($(_FC),xlf)
        ifeq ($(TARGET),LINUX)
.F.f:
	$(FC) -c $(FFLAGS) -E $(CPPFLAGS) $< -o $*.f

        endif
    endif

# else for ifndef Flint
else
#        -------------

# First time thru you need the -L... option. Next time remove it
# and move the library file name to the last but one argument and add -g
    ifeq ($(FLINT),1)
.F.o:; flint $(CPPFLAGS) -L $(SRCDIR)/nwchem.lbt $<

.f.o:; flint $(CPPFLAGS) -L $(SRCDIR)/nwchem.lbt $<

    else
.F.o:; flint $(CPPFLAGS) -g -f -u $(SRCDIR)/nwchem.lbt $<

.f.o:; flint $(CPPFLAGS) -g -f -u $(SRCDIR)/nwchem.lbt $<

    endif

endif
