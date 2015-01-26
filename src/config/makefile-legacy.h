
# $Id: makefile.h 26564 2014-12-22 17:20:10Z jhammond $
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
RELEASE := 
# For current release tree
#RELEASE := 6.5

#

ifndef NWCHEM_TOP
error1:
$(info     )
$(info You must define NWCHEM_TOP in your environment to be the path)
$(info of the top level nwchem directory ... something like)
$(info     setenv NWCHEM_TOP /msrc/home/elvis/nwchem)
$(info     )
$(error )
endif

# Select the old (pre-autotools version of GA) by uncommenting the next line.
# The value of OLD_GA does not matter -- it is detected as an ifdef only.
#OLD_GA = y 

#
# Do a setenv for NWCHEM_TARGET to be the machine and NWCHEM_TARGET_CPU the CPU to build for
#
# NWCHEM_TARGET :  
#                  CYGNUS       (Windows under Cygwin tools)
#                  DECOSF
#                  IBM
#                  LINUX        NWCHEM_TARGET_CPU :
#                                                  nothing for X86 (e.g. do not set this)
#                                                  ALPHA for AlphaLinux (broke)
#                                                  POWERPC for MkLinux 
#                  SGI
#                  SGI_N32      NWCHEM_TARGET_CPU : R8000 or R10000
#                  SGITFP       NWCHEM_TARGET_CPU : R8000 or R10000
#                  SOLARIS      NWCHEM_TARGET_CPU : not defined or ULTRA
#                  LAPI         NWCHEM_TARGET_CPU : P2SC
#                      (uses thread safe libraries and LAPI)
#
#

ifndef NWCHEM_TARGET
error2:
	@echo You must define NWCHEM_TARGET in your environment to be the name
	@echo of the machine you wish to build for ... for example
	@echo     setenv NWCHEM_TARGET SOLARIS
	@echo Known targets are SOLARIS, SGI_N32, ...
	@echo See the INSTALL instructions for a complete list
	@exit 2
endif

     TARGET := $(NWCHEM_TARGET)
ifeq (,$(RELEASE))
     TOPDIR := $(NWCHEM_TOP)
     CODE_BRANCH := Development
else
     ifeq (,$(findstring $(RELEASE),$(NWCHEM_TOP)))
          TOPDIR := $(NWCHEM_TOP)-$(RELEASE)
     else
          TOPDIR := $(NWCHEM_TOP)
     endif
     CODE_BRANCH := $(RELEASE)
endif

#dummy:
#	echo NWCHEM_TOP=$(NWCHEM_TOP) RELEASE=$(RELEASE) TOPDIR=$(TOPDIR)

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
    LIBPATH = -L$(SRCDIR)/tools/install/lib
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
    INCPATH = -I$(SRCDIR)/tools/install/include
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

NW_CORE_SUBDIRS = tools include basis geom inp input  \
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
#       MAKE = DON'T define this ... it will break the passing of command
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
#                 are needed for top-level modules on this machine.
#		  (Should not normally be used)
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
#           SGITFP                       R10000/R8000
#           SGI_N32                      R10000/R8000
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
# you'll need to specify -DINTEGER_1='integer*1' in the selci
# and util makefiles.
#
         CC = cc
         FC = f77

# Don't need this if using the SUN performance library
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
#    Fujitsu SPARC systems (thanks to Herbert Fruchtl)
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
    CORE_LIBS +=  -llapack $(BLASOPT) -lblas  -lmvec
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
# you'll need to specify -DINTEGER_1='integer*1' in the selci
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
#    Fujitsu SPARC systems (thanks to Herbert Fruchtl)
    COPTIONS = -Kdalign -KV9FMADD
    COPTIMIZE = -Kfast_GP=2 -KV9FMADD
    DEFINES += -DFUJITSU_SOLARIS
  else
# SUN/Solaris options for WS6.1
   COPTIONS = -xarch=v9 -dalign
  endif

  ifeq ($(FC),frt)
#    Fujitsu SPARC systems (thanks to Herbert Fruchtl)
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
#  SUN/Solaris f77 options 
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
    CORE_LIBS +=  -llapack -lblas
  else
    LDOPTIONS = -xs -xildoff
ifdef BLASOPT
    CORE_LIBS +=   $(BLASOPT)   -lmvec
else
    CORE_LIBS +=  -llapack -lblas  -lmvec
endif
    CORE_LIBS += -lsocket -lrpcsvc -lnsl
    EXTRA_LIBS =  -ldl -lfsu
  endif


#end of solaris
endif

ifeq ($(TARGET),PURESOLARIS)
#
# NOT TESTED RECENTLY
#
# Sun running Solaris 2.4 or later and if you want to use purecoverage tool you must
#
      SHELL := $(NICE) /bin/sh
    CORE_SUBDIRS_EXTRA = blas lapack
         CC = purecov gcc
         FC = purecov f77
     RANLIB = echo
  MAKEFLAGS = -j 2 --no-print-directory
    INSTALL = echo $@ is built
# -fast introduces many options that must be applied to all files
# -stackvar puts locals on t
# the stack which seems a good thing
#     but may need to increase the stacksize at runtime using limit
# -xs allows debugging without .o files
   FOPTIONS = -Nl199 -fast -dalign -stackvar
   COPTIONS = -Wall
# Under Solaris -O3 is the default with -fast (was -O2 with SUNOS)
# -fsimple=2 enables more rearranging of floating point expressions
# -depend enables more loop restructuring
  FOPTIMIZE = -O3 -fsimple=2 -depend 
# Under Solaris -g no longer disables optimization ... -O2 seems solid
# but is slow and impairs debug ... use -O1 for speed and debugability
     FDEBUG = -g -O1
  COPTIMIZE = -g -O1
   LIBPATH += -L/usr/ucblib
   LIBPATH += -L/afs/msrc/sun4m_54/apps/purecov
   OPTIONS = -xildoff -Bstatic
   CORE_LIBS +=  -llapack $(BLASOPT) -lblas
# First four needed for parallel stuff, last for linking with profiling
	   EXTRA_LIBS = -lsocket -lrpcsvc -lnsl -lucb -lintl -lc -lc -lpurecov_stubs

#end of puresolaris
endif


ifeq ($(TARGET),CRAY-T3E)
#
#
   CORE_SUBDIRS_EXTRA = blas lapack # Only a couple of routines not in scilib
               RANLIB = echo
            MAKEFLAGS = -j 1 --no-print-directory
              INSTALL = @echo $@ is built
        OUTPUT_OPTION =

                   FC = f90
                  CPP = /opt/ctl/CC/CC/lib/mppcpp -P  -N
             FOPTIONS = -d p -F 
             COPTIONS =
               FDEBUG = -O scalar1
#               FDEBUG = -g
            FOPTIMIZE = -O scalar3,aggress,unroll2,vector3
#,pipeline3
               CDEBUG = -O 1
            COPTIMIZE = -O
#
# to debug code you must remove the -s flag unless you know assembler
#
#            LDOPTIONS = -g -Xm  -lmfastv
            LDOPTIONS = -Wl"-s" -Xm  -lmfastv

              DEFINES = -DCRAY_T3E -DCRAY_T3D -D__F90__ -DUSE_FCD

               LINK.f = f90 $(LDFLAGS)

            CORE_LIBS += -llapack $(BLASOPT) -lblas
#
# 
ifeq ($(BUILDING_PYTHON),python)
#** on the NERSC CRAY-T3E you need to:
#**** % module load python
#**** % module load tcltk
      EXTRA_LIBS += -ltk -ltcl -lX11
endif
#

      FCONVERT      = $(CPP) $(CPPFLAGS)  $< | sed '/^\#/D'  > $*.f
      EXPLICITF     = TRUE
endif

ifeq ($(TARGET),SGITFP)
#
# SGI power challenge
#
# CORE_SUBDIRS_EXTRA are those machine specific libraries required 
#
# TPS 95/11/22:
# Optimization options const_copy_limit=18000, global_limit=18000 and 
# fprop_limit=1200 added to FOPTIMIZE to allow full optimization of the 
# MD module nwArgos on SGI Power Indigo^2
#
# RJH ... note that fprop_limit is not supported by 7.0 compilers
#
# TPS 95/12/12:
# Increased fprop_limit to 1750
# Removed -j 12 from MAKEFLAGS
# Added -lnwcutil to core libraries
#
# TPS 96/01/14:
# Increased const_copy_limit and global_limit to 18500
#
# RJH ... from Roberto ... on the R10k TENV=3 may cause very expensive
#     interrupts (he recommends 1 when we go to 10K, but 3 is good for 8k)
#     ... also going to 10K use -SWP:if_conversion=OFF
#     ... in going to 6.1/2 then should also set -SWP:*ivdep*=ON/OFF
#         (default changed from agressive to conservative and we want
#         the agressive)
#     ... sometimes the default KAP parameters are best (only for critical
#         routines)
#     ... -dr=AKC forces it to recognize all compiler directives (C=CRAY
#         not on by default) ... can put everywhere.
#     ... on -WK also add -r=3 (level of reduction) even with -o=1
#     ... could benefit from -Wk on FOPTIMIZE ... actually have it on now.
#     ... roundoff/ieee only modify pipelining which happens only at O3
#
# TPS 96/06/27:
# Added -lnwcutil to core libraries (again!)
#
# TPS 96/07/26:
# Fortran optimization limits: const_copy_limit=20000 
#                              global_limit=20000
#                              fprop_limit=2000
#
# JN 96/10/02:
# Replaced -DLongInteger with -DEXT_INT for consistency with GA, DRA, PEIGS ...
#
# JN 99/05/26: MA now has its own library -lma
#
# TLW 99/10/08:
# From Gerardo Cisneros
#  - took out obsolete options
#       "-OPT:fold_arith_limit=4000"
#       "-OPT:fprop_limit=2000"
#       "-OPT:global_limit=20000"
#       "-SWP:if_conversion=OFF"


	CPP = /usr/lib/cpp
  CORE_SUBDIRS_EXTRA = blas lapack
         FC = f77
     RANLIB = echo


    INSTALL = @echo nwchem is built
  MAKEFLAGS = -j 4 --no-print-directory

# RJH ... moved -OPT... to the FOPTIMIZE macro since it's an optimization!
# (that breaks things).
  FOPTIONS = -d8 -i8 -mips4 -align64 -64 -r8 -G 0 
  COPTIONS = -fullwarn -mips4 -64

#optimization flags for R8000 (IP21)
# FOPTIMIZE_8K = -O3 -OPT:fold_arith_limit=4000:const_copy_limit=20000:global_limit=20000:fprop_limit=2000 -TENV:X=3 -WK,-so=1,-o=1,-r=3,-dr=AKC
 FOPTIMIZE_8K = -O3 -OPT:const_copy_limit=20000 -TENV:X=3 -WK,-so=1,-o=1,-r=3,-dr=AKC
FVECTORIZE_8K = -O3 -TENV:X=3 -WK,-dr=AKC

#optimization flags for R10000 (IP28)
 FOPTIMIZE_10K = -O3 -OPT:const_copy_limit=20000:Olimit=4800 -TENV:X=1 -WK,-so=1,-o=1,-r=3,-dr=AKC
FVECTORIZE_10K = -O3 -TENV:X=1 -WK,-dr=AKC

#optimization flags for R12000 (IP30)
 FOPTIMIZE_12K = -O3 -OPT:const_copy_limit=20000:Olimit=4800 -TENV:X=1 -WK,-so=1,-o=1,-r=3,-dr=AKC
FVECTORIZE_12K = -O3 -TENV:X=1 -WK,-dr=AKC

 COPTIMIZE = -O
 FOPTIMIZE = -O3

ifeq ($(NWCHEM_TARGET_CPU),R12000)
 FOPTIMIZE = $(FOPTIMIZE_12K)
 FVECTORIZE = $(FVECTORIZE_12K)
endif
ifeq ($(NWCHEM_TARGET_CPU),R10000)
 FOPTIMIZE = $(FOPTIMIZE_10K)
 FVECTORIZE = $(FVECTORIZE_10K)
endif
ifeq ($(NWCHEM_TARGET_CPU),R8000)
 FOPTIMIZE = $(FOPTIMIZE_8K)
 FVECTORIZE = $(FVECTORIZE_8K)
endif
  FOPTIMIZE += -OPT:roundoff=3:IEEE_arithmetic=3
  FVECTORIZE += -OPT:roundoff=3:IEEE_arithmetic=3

  DEFINES = -DSGI -DSGITFP -DEXT_INT
  CORE_LIBS += -llapack $(BLASOPT) -lblas
ifeq ($(BUILDING_PYTHON),python)
#needed for python 2.2.2
      EXTRA_LIBS += -lpthread
endif
endif




ifeq ($(TARGET),SGI_N32)
#
# SGI 64-bit MIPS-4 processors (R5k, R8k, R10k) under IRIX > 6.0  (ABI)
#
# JN, 12.06.96:
# -n32 allows to use 64-bit processor features and 32-bit address space
# 32-bit address space - use SGITFP if 64-bit addresses needed
#
# SGI BLAS can be used directly

    CORE_SUBDIRS_EXTRA = lapack
         FC = f77
         CPP = /usr/lib/cpp
         AR = ar
     RANLIB = echo

    INSTALL = @echo nwchem is built
  MAKEFLAGS = -j 4 --no-print-directory
    DEFINES = -DSGI  -DSGI_N32

# RJH ... moved -OPT... to the FOPTIMIZE macro since it's an optimization!
# (that breaks things).
  FOPTIONS = -n32 -mips4 -G 0 
  COPTIONS = -n32 -mips4 -fullwarn

#optimization flags for R8000 (IP21)
 FOPTIMIZE_8K = -O3 -OPT:const_copy_limit=20000 -TENV:X=3 -WK,-so=1,-o=1,-r=3,-dr=AKC
FVECTORIZE_8K = -O3 -TENV:X=3 -WK,-dr=AKC

#optimization flags for R10000 (IP28)
 FOPTIMIZE_10K = -O3 -OPT:const_copy_limit=20000 -TENV:X=1 -WK,-so=1,-o=1,-r=3,-dr=AKC 
FVECTORIZE_10K = -O3 -TENV:X=1 -WK,-dr=AKC

#optimization flags for R12000 (IP27)
# FOPTIMIZE_12K = -O   -r12000  -TARG:platform=ip27  -LNO:cs2=8M -TENV:X=3
# FOPTIMIZE_12K = -O3 -r12000 -TARG:platform=ip27 -LNO:prefetch=1:cs2=8M:fusion=2:fission=2 -LIST:all_options -OPT:swp=ON:space=ON 
# FVECTORIZE_12K = -O3 -r12000 -TARG:platform=ip27 -LNO:prefetch=1:cs2=8M:fusion=2:fission=2 -LIST:all_options -OPT:swp=ON:space=ON 
# The above options are some Edo used to optimize for a particular machine
 FOPTIMIZE_12K = -O3 -OPT:const_copy_limit=20000 -TENV:X=1 -WK,-so=1,-o=1,-r=3,-dr=AKC 
FVECTORIZE_12K = -O3 -TENV:X=1 -WK,-dr=AKC

 FOPTIMIZE = -O3
 COPTIMIZE = -O2

ifeq ($(NWCHEM_TARGET_CPU),R12000)
 FOPTIMIZE = $(FOPTIMIZE_12K)
 FVECTORIZE = $(FVECTORIZE_12K)
endif
ifeq ($(NWCHEM_TARGET_CPU),R10000)
 FOPTIMIZE = $(FOPTIMIZE_10K)
 FVECTORIZE = $(FVECTORIZE_10K)
endif
ifeq ($(NWCHEM_TARGET_CPU),R8000)
 FOPTIMIZE = $(FOPTIMIZE_8K)
 FVECTORIZE = $(FVECTORIZE_8K)
endif
  FOPTIMIZE += -OPT:roundoff=3:IEEE_arithmetic=3
  FVECTORIZE += -OPT:roundoff=3:IEEE_arithmetic=3

ifeq ($(BUILDING_PYTHON),python)
# needed if python was compiled with gcc (common)
      EXTRA_LIBS += -L/msrc/apps/gcc-2.8.1/lib/gcc-lib/mips-sgi-irix6.5/2.8.1 -lgcc
# needed here if using a python version with tk/tcl extensions  (common)
      EXTRA_LIBS += -L/msrc/apps/lib -ltk8.0 -ltcl8.0 
# needed here if using a python version built with BLT extensions
#     EXTRA_LIBS += -L/msrc/apps/lib -lBLT 
# Both tk/tcl and BLT need X11 (common)
      EXTRA_LIBS += -lX11
endif

       CORE_LIBS += -llapack $(BLASOPT) -lblas 
ifeq ($(BUILDING_PYTHON),python)
#needed for python 2.2.2
      EXTRA_LIBS += -lpthread
endif
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
  CORE_LIBS +=  $(BLASOPT) -llapack -lblas   -lm
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
  CORE_LIBS +=  -llapack $(BLASOPT) -lblas  -lm
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
# -qinitauto=FF
   COPTIONS = 
# -qstrict required with -O3 (according to Edo)
# -qfloat=rsqrt gives faster square roots (off by -qstrict)
# -qfloat=fltint gives faster real-integer conversion (off by -qstrict)
# -qhot seems to break a lot of things so don't ever use it
  FOPTIMIZE = -O3 -qstrict -NQ40000 -NT80000 -qarch=auto -qtune=auto -NS2048
  ifdef RSQRT
    FOPTIMIZE  += -qfloat=rsqrt:fltint
  endif
  COPTIMIZE = -O -qarch=auto -qtune=auto
  ifdef  USE_GPROF
    FOPTIONS += -pg
    LDOPTIONS += -pg
  endif

    DEFINES = -DIBM -DAIX -DEXTNAME
   CORE_LIBS +=  $(BLASOPT) 
ifdef USE_ESSL
   DEFINES += -DESSL
   CORE_LIBS += -lessl
endif

       LIBPATH += -L/usr/lib 

  LDOPTIONS += -bmaxstack:0x60000000 -bmaxdata:0x60000000 -bloadmap:nwchem.lapi_map
       CORE_LIBS +=  -llapack $(BLASOPT) -lblas \
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
   ifeq ($(BLAS_LIB),)
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
    FDEBUG = -O2 -qmaxmem=8192
  ifdef RSQRT
    FOPTIMIZE  += -qfloat=rsqrt:fltint
  endif
    XLF8= $(shell /usr/bin/lslpp -l xlfcmp  2>&1|grep COMM|head -n 1| awk ' / [8-9]./  {print "Y"};/[ ][1][0-9]./  {print "Y"}')
  ifdef XLF8
    FVECTORIZE= -O3 -qstrict -qtune=auto -qarch=auto -qcache=auto -qalign=natural -qnozerosize -qlargepage -qnozerosize -qipa=level=2
#old    FOPTIMIZE = -O4  -NQ40000 -NT80000  -qarch=auto -qtune=auto
    FOPTIMIZE = -O3   -qarch=auto -qtune=auto
#adding -qstrict to fix linking problem for v10.1
    FOPTIMIZE  += -qstrict# -qipa -qhot -qlargepage -qessl 
    FOPTIONS += -blpdata  
    FOPTIMIZE  += -qfloat=rsqrt:fltint
    FVECTORIZE  += -qfloat=rsqrt:fltint
  endif
   COPTIMIZE = -O -qmaxmem=8192

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
  LDOPTIONS += -bloadmap:nwchem.ibm64map -bbigtoc # bigtoc requires bmaxdata
  LDOPTIONS += -bmaxstack:0x80000000 -bmaxdata:0x200000000 # this limits MA to 8GB
  ifeq ($(LAPACK_LIB),)
     CORE_LIBS += -llapack
  endif
  CORE_LIBS += $(BLASOPT)
  ifeq ($(BLAS_LIB),)
     CORE_LIBS += -lblas
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

CORE_LIBS +=  -llapack -lblas

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
    FVECTORIZE= -O3 -qstrict -qtune=auto -qarch=auto -qcache=auto -qalign=natural -qnozerosize -qlargepage -qnozerosize -qipa=level=2
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
  CORE_LIBS +=  $(BLASOPT) -llapack -lblas
  LDOPTIONS += -bloadmap:nwchem.lapi64map -bbigtoc
  LDOPTIONS += -bmaxstack:0x80000000 -bmaxdata:0x80000000 # needed because of bigtoc
#  LDOPTIONS += -bmaxstack:0xe0000000 -bmaxdata:0xe0000000 # this for large memory
  XLFBREN = y

 EXPLICITF = TRUE
#
endif

ifeq ($(TARGET),DECOSF)
#
# DEC AXP OSF1
#
# JN 96/10/02:
# Replaced -DLongInteger with -DEXT_INT for consistency with GA, DRA, PEIGS ...

    CORE_SUBDIRS_EXTRA = blas lapack
#                  NICE = nice
#                SHELL := $(NICE) /bin/sh
                    FC = f77
                    AR = ar
                RANLIB = echo
	 	   CPP = /usr/bin/cpp -P -C	

               INSTALL = @echo nwchem is built
             MAKEFLAGS = -j 1 --no-print-directory

# -fpe2 and call to util/dec_fpe.f from nwchem.F necessary to avoid
# braindead alpha undflows inside texas (c6h6 6-31g)

# assume noaccuracy_sensitive was breaking the code in recent versions (EA)
              FDEBUG = -g -O0
              FOPTIONS = -align dcommons -math_library fast -fpe2 -check nounderflow -check nopower -check nooverflow  -warn argument_checking -warn unused -automatic -math_library fast

             COPTIONS = 
             LDOPTIONS = -O
             LINK.f = f77 $(LDFLAGS)
             FOPTIMIZE =  -O4  -tune host -arch host  
             FVECTORIZE = -fast -O4 -tune host -arch host
             COPTIMIZE = -O

               DEFINES = -DDECOSF
ifdef USE_I4FLAGS
              FOPTIONS +=  -i4
else
              FOPTIONS +=  -i8
endif
               DEFINES +=  -DEXT_INT 
             CORE_LIBS +=  -llapack $(BLASOPT) -lblas 
            EXTRA_LIBS = -laio 
ifeq ($(BUILDING_PYTHON),python)
      EXTRA_LIBS += -lX11
endif
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
    FOPTIONS +=  -NQ40000 -NT80000 -NS2048 -qmaxmem=8192 -qxlf77=leadzero
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
        FOPTIMIZE  = -O2 -ffast-math -Wuninitialized 
        DEFINES  += -DGFORTRAN
        GNUMAJOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __VERS | cut -c22)
        ifdef GNUMAJOR
        GNUMINOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __VERS | cut -c24)
        GNU_GE_4_6 = $(shell [ $(GNUMAJOR) -gt 4 -o \( $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 6 \) ] && echo true)
        GNU_GE_4_8 = $(shell [ $(GNUMAJOR) -gt 4 -o \( $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 8 \) ] && echo true)
        endif
        ifeq ($(GNU_GE_4_6),true)
        DEFINES  += -DGCC46
        endif
        ifeq ($(GNU_GE_4_8),true)
          FDEBUG += -fno-aggressive-loop-optimizations
          FOPTIMIZE +=-fno-aggressive-loop-optimizations
          FFLAGS_FORGA += -fno-aggressive-loop-optimizations
          
          FOPTIONS += -Warray-bounds
        endif
        ifdef USE_OPENMP
           FOPTIONS  += -fopenmp
           LDOPTIONS += -fopenmp
           DEFINES += -DUSE_OPENMP
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
#         FOPTIONS= -malign-double# this break with gfort 4.2 and later http://gcc.gnu.org/bugzilla/show_bug.cgi?id=29562
         FOPTIMIZE+= -funroll-all-loops -mtune=native 
         FVECTORIZE=-O3 -ffast-math -mtune=native -mfpmath=sse -msse3 -ftree-vectorize -ftree-vectorizer-verbose=1   -fprefetch-loop-arrays  -funroll-all-loops 
#         FOPTIMIZE=-O1
#         FVECTORIZE=-O1
       endif
        ifdef USE_F2C
#possible segv with use of zdotc (e.g. with GOTO BLAS)
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=20178
          FOPTIONS +=  -ff2c -fno-second-underscore
        endif
        DEFINES  += -DCHKUNDFLW -DGCC4
      endif
      ifeq ($(FC),ifort)
    _FC=ifort
#ifort 9.1
#        LINK.f = ifort  $(LDFLAGS) 
        FOPTIONS   += -align    -mp1 -w -g -vec-report1
  ifdef  USE_GPROF
    FOPTIONS += -qp
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
ifdef USE_VECLIB
             CORE_LIBS += $(BLASOPT)  -Wl,-framework -Wl,vecLib -lblas
else
             CORE_LIBS +=   -llapack $(BLASOPT)  -lblas
endif
  ifeq ($(FC),xlf) 
     LDOPTIONS = -Wl,-multiply_defined -Wl,warning
  else
#  _GCC4= $(shell gcc -v  2>&1|egrep spec|head -n 1|awk ' / 3./  {print "N";exit}; / 2./ {print "N";exit};{print "Y"}')
  _GCC4= $(shell $(CC) -dM -E - < /dev/null | egrep __VERS | cut -c22|awk ' /3/  {print "N";exit}; /2/ {print "N";exit};{print "Y"}')
    ifeq ($(_GCC4),Y) 
#      EXTRA_LIBS += 
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
    ifndef FC
      @echo Defaulting to FC=gfortran
      FC=gfortran
    endif
    ifndef _FC
      _FC=$(FC)
    endif
#
# MacOSX 64bit
#
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
   ifeq ($(BLAS_LIB),)
      CORE_SUBDIRS_EXTRA += blas
   endif
   ifeq ($(LAPACK_LIB),)
      CORE_SUBDIRS_EXTRA += lapack
   endif
endif
               _CPU = $(shell machine  )
               INSTALL = @echo nwchem is built
               RANLIB = ranlib
             MAKEFLAGS = -j 1 --no-print-directory
             DEFINES   = -DMACX
             DEFINES  += -DEXT_INT

      ifeq ($(FC),gfortran)
#gcc version 
        LINK.f = gfortran  $(LDFLAGS) 
        FOPTIONS   = #-Wextra #-Wunused #-ffast-math
        FOPTIONS += -fdefault-integer-8
        FOPTIMIZE  = -O2 -ffast-math -Wuninitialized 
       DEFINES  += -DGFORTRAN -DGCC4
#
         FOPTIMIZE+= -funroll-all-loops -mtune=native 
         FVECTORIZE=-O3 -ffast-math -mtune=native -mfpmath=sse -msse3 -ftree-vectorize -ftree-vectorizer-verbose=1   -fprefetch-loop-arrays  -funroll-all-loops 
#         FOPTIMIZE=-O1
#         FVECTORIZE=-O1
        GNUMAJOR=$(shell $(FC) -dumpversion | cut -f1 -d.)
        GNUMAJOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __VERS | cut -c22)
        ifdef GNUMAJOR
        GNUMINOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __VERS | cut -c24)
        GNU_GE_4_6 = $(shell [ $(GNUMAJOR) -gt 4 -o \( $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 6 \) ] && echo true)
        GNU_GE_4_8 = $(shell [ $(GNUMAJOR) -gt 4 -o \( $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 8 \) ] && echo true)
        ifeq ($(GNU_GE_4_6),true)
         DEFINES  += -DGCC46
        endif
        ifeq ($(GNU_GE_4_8),true)
          FDEBUG += -fno-aggressive-loop-optimizations
          FOPTIMIZE +=-fno-aggressive-loop-optimizations
          FFLAGS_FORGA += -fno-aggressive-loop-optimizations
          FOPTIONS += -Warray-bounds
        endif # GNU_GE_4_8
        endif # GNUMAJOR

        ifdef USE_OPENMP
           FOPTIONS  += -fopenmp
           LDOPTIONS += -fopenmp
           DEFINES += -DUSE_OPENMP
        endif
       endif # gfortran
    ifdef  USE_GPROF
      FOPTIONS += -pg
      LDOPTIONS += -pg
      COPTIONS += -pg
    endif
ifdef USE_VECLIB
             CORE_LIBS += $(BLASOPT)  -Wl,-framework -Wl,vecLib -lblas
else
   ifeq ($(LAPACK_LIB),)
      CORE_LIBS += -llapack
   endif
   CORE_LIBS +=    $(BLASOPT)
   ifeq ($(BLAS_LIB),)
      CORE_LIBS += -lblas
   endif
endif
      ifeq ($(FC),ifort)
       _IFCV11= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 11) {print "Y";exit}}')
       _IFCV12= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 12) {print "Y";exit}}')
       _IFCV14= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 14) {print "Y";exit}}')
       _IFCV15ORNEWER=$(shell ifort -logo  2>&1|egrep "Version "|head -n 1 | sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 15) {print "Y";exit}}')
        DEFINES  += -DIFCV8 -DIFCLINUX
        FOPTIONS += -i8
        FOPTIONS +=  -g -no-save-temps
        FDEBUG    = -O2 -g
        FOPTIMIZE = -O3 -xHost
        ifdef USE_OPENMP
           FOPTIONS  += -openmp
           LDOPTIONS += -openmp
           DEFINES   += -DUSE_OPENMP
        endif
        ifeq ($(_IFCV11),Y) 
#next 2 lines needed for fp accuracy
          FDEBUG += -fp-model source
          ifeq ($(_IFCV12),Y) 
            FOPTIONS += -fimf-arch-consistency=true
          endif
        endif
      endif

#  _GCC4= $(shell gcc -v  2>&1|egrep spec|head -n 1|awk ' / 3./  {print "N";exit}; / 2./ {print "N";exit};{print "Y"}')
  _GCC4= $(shell $(CC) -dM -E - < /dev/null | egrep __VERS | cut -c22|awk ' /3/  {print "N";exit}; /2/ {print "N";exit};{print "Y"}')
    ifeq ($(_GCC4),Y) 
#      EXTRA_LIBS += 
    else
      EXTRA_LIBS += -lm -lcc_dynamic
    endif

# required for mpich2 3.x and clang
    COPTIONS +=-DMPICH_NO_ATTR_TYPE_TAGS
    CFLAGS_FORGA +=-DMPICH_NO_ATTR_TYPE_TAGS
#

endif


ifeq ($(TARGET),$(findstring $(TARGET),LINUX CYGNUS CYGWIN INTERIX))
#
#
# Linux or Cygwin under Windows running on an x86 using g77
#
       NICE = nice -n 2
      SHELL := $(NICE) /bin/sh
    CORE_SUBDIRS_EXTRA = blas lapack
         CC = gcc
     RANLIB = ranlib
  MAKEFLAGS = -j 1 --no-print-directory
    INSTALL = @echo $@ is built
        CPP = gcc -E -nostdinc -undef -P
   FCONVERT = (/bin/cp $< /tmp/$$$$.c; \
			$(CPP) $(CPPFLAGS) /tmp/$$$$.c | sed '/^$$/d' > $*.f; \
			/bin/rm -f /tmp/$$$$.c) || exit 1

         FC=gfortran
         LINUXCPU = $(shell uname -m |\
                 awk ' /sparc/ { print "sparc" }; /i*86/ { print "x86" };  /ppc*/ { print "ppc"} ' )

     GOTMINGW32= $(shell $(CC) -dM -E - </dev/null 2> /dev/null |grep MINGW32|cut -c21)

ifeq ($(BUILDING_PYTHON),python)
#   EXTRA_LIBS += -ltk -ltcl -L/usr/X11R6/lib -lX11 
#   EXTRA_LIBS += -L/home/edo/tcltk/lib/LINUX -ltk8.3 -ltcl8.3 -L/usr/X11R6/lib -lX11 -ldl
# needed if python was built with pthread support
  ifneq ($(GOTMINGW32),1)
   EXTRA_LIBS += $(shell $(PYTHONHOME)/bin/python-config --libs) -lz -lnwcutil
  endif
endif

  DEFINES = -DLINUX

      FOPTIMIZE  = -O2 
      COPTIONS   = -Wall
      COPTIMIZE  = -g -O2
      ifeq ($(FC),gfortran)
        FOPTIONS   = # -Wextra -Wunused  
        FOPTIMIZE  += -ffast-math -Wuninitialized
        _FC=gfortran
        DEFINES  += -DGFORTRAN
        GNUMAJOR=$(shell $(FC) -dumpversion | cut -f1 -d.)
        GNUMAJOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __VERS | cut -c22)
        ifdef GNUMAJOR
          GNUMINOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __VERS | cut -c24)
          GNU_GE_4_6 = $(shell [ $(GNUMAJOR) -gt 4 -o \( $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 6 \) ] && echo true)
          GNU_GE_4_8 = $(shell [ $(GNUMAJOR) -gt 4 -o \( $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 8 \) ] && echo true)
          ifeq ($(GNU_GE_4_6),true)
            DEFINES  += -DGCC46
          endif
          ifeq ($(GNU_GE_4_8),true)
            FDEBUG +=-O2 -g -fno-aggressive-loop-optimizations
            FOPTIMIZE +=-fno-aggressive-loop-optimizations
            FFLAGS_FORGA += -fno-aggressive-loop-optimizations
          endif
         endif
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
    FOPTIMIZE += -Wuninitialized -ffast-math -funroll-loops -fstrength-reduce 
    FOPTIMIZE += -fno-move-all-movables -fno-reduce-all-givs 
    FOPTIMIZE += -fforce-addr 
# see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=13037
#  for atomscf/orderd.f  (bug report by Kirill Smelkov)
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
      COPTIONS   =  -march=i686 
      ifdef USE_GCC31
        FDEBUG=-O1 -g
        COPTIMIZE +=-march=pentium4 -mcpu=pentium4 #-msse2 -mfpmath=sse 
#        COPTIMIZE +=-fprefetch-loop-arrays -minline-all-stringops -fexpensive-optimizations
        FOPTIMIZE +=-march=pentium4 -mcpu=pentium4# -msse2 -mfpmath=sse 
#        FOPTIMIZE +=-fprefetch-loop-arrays -minline-all-stringops -fexpensive-optimizations
      else
#        FOPTIMIZE  += -march=i686
        COPTIONS   = -Wall -march=i686 -malign-double 
      endif
    else
    ifneq ($(_CPU),x86)
      COPTIONS   +=  -march=$(_CPU)
      FOPTIONS   +=  -march=$(_CPU)
    endif
    endif
    ifeq ($(_CPU),k7)
       FOPTIONS   = -fno-second-underscore  
       COPTIONS   = -Wall -malign-double
       ifdef  USE_GCC31
        FOPTIONS += -march=athlon
        COPTIONS += -march=athlon
       else
        FOPTIONS += -march=k6
        COPTIONS += -march=k6
       endif
    endif

  ifeq ($(FC),pgf77)
    DEFINES   += -DPGLINUX
# added -Kieee to get dlamc1 to work on pgf77 3.1-3 EA Jun 8th 2000
    FOPTIONS   = -Mdalign -Minform,warn -Mnolist -Minfo=loop -Munixlogical -Kieee
    ifeq ($(_CPU),i586)
      FOPTIONS  += -tp p5  
    endif
    ifeq ($(_CPU),i686)
      FOPTIONS  += -tp p6 -Mvect=prefetch
    endif
    ifeq ($(_CPU),i786)
      FOPTIONS  += -tp piv  -Mcache_align  -Mvect=prefetch
    endif
    FOPTIMIZE  = -O2 -Mvect=assoc,cachesize:262144 -Munroll -Mnoframe
  endif
# _FC=g77
 ifeq ($(FC),ifc)
     _FC=ifc
 endif
 ifeq ($(FC),ifort)
     _FC=ifc
 endif
  ifeq ($(_FC),ifc)
  FOPTIONS   =  -align    -mp1 -w -g -vec-report1
  ifdef  USE_GPROF
    FOPTIONS += -qp
  endif
    _IFCV7= $(shell ifort -v  2>&1|egrep "Version "|head -n 1|awk '/7./ {print "Y"; exit}')
    _IFCV10= $(shell ifort -v  2>&1|egrep "Version "|head -n 1|awk '/10./ {print "Y"; exit}')
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
    DEFINES   += -DIFCLINUX
    ifneq ($(_IFCV7),Y)
      FOPTIMIZE += -ansi_alias-
    endif
  endif
      ifeq ($(_FC),gfortran)
        LINK.f = gfortran  $(LDFLAGS) 
        FOPTIONS  += -m32
        COPTIONS  += -m32
        CFLAGS_FORGA += -m32
        FFLAGS_FORGA += -m32
        FOPTIMIZE  += -O2 -ffast-math -Wuninitialized
        ifeq ($(_CPU),i786)
          FOPTIONS += -march=pentium4 -mtune=pentium4
          FVECTORIZE = $(FOPTIMIZE) -O3 -ftree-vectorize 
          FVECTORIZE += -ftree-vectorizer-verbose=1
#        FOPTIMIZE  += -fprefetch-loop-arrays -ftree-loop-linear
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
    COPTIONS   =   -mp1 -w -g -vec-report1
    COPTIMIZE = -O3   -unroll 
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

  ifeq ($(LINUXCPU),ppc)
# this are for PowerPC
# Tested on SLES 9
# Feb 7th 2005
# xlf v9.1
# xlc v7.0 
# gcc-3.2.3-42 
    ifeq ($(FC),xlf)
      _FC=xlf
    endif
    ifeq ($(FC),blrts_xlf)
      _FC=xlf
    endif
    ifeq ($(_FC),xlf)
      FOPTIONS  = -q32  -qextname -qfixed 
      FOPTIONS +=  -NQ40000 -NT80000 -NS2048 -qmaxmem=8192 -qxlf77=leadzero
      FOPTIMIZE= -O3 -qstrict -qfloat=fltint
      ifeq ($(FC),blrts_xlf)
        FOPTIMIZE+= -qarch=440 -qtune=440
      else
        FOPTIMIZE+= -qarch=auto -qtune=auto
      endif
      FDEBUG= -O2 -g
      EXPLICITF = TRUE
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
#  LDOPTIONS = --Xlinker -O -Xlinker -static
      EXTRA_LIBS += -lm
   endif
    endif
  endif
endif
#EXTRA_LIBS +=-lefence # link against Electricfence

CORE_LIBS += -llapack $(BLASOPT) -lblas

# end of Linux, Cygnus
endif

ifneq ($(TARGET),LINUX)
ifeq ($(TARGET),$(findstring $(TARGET),LINUX64 CYGWIN64 CATAMOUNT))
     _CPU = $(shell uname -m  )
#ifeq ($(NWCHEM_TARGET),LINUX64)
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
          _CC=craycc
	  endif
          DEFINES  += -DCRAYXT -DNOIO
          USE_NOIO=1
      endif
      ifeq ($(CC),gcc)
        _CC=gcc
      endif
      ifeq ($(CC),pgcc)
        _CC=pgcc
      endif
      ifeq ($(CC),icc)
        _CC=icc
      endif
      ifeq ($(FC),pgf90)
        _FC=pgf90
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
      ifeq ($(FC),gfortran)
       _FC=gfortran
      endif
      ifndef _FC
      FC=gfortran
      _FC=gfortran
      endif
      ifndef _CC
      _CC=gcc
      endif
      FOPTIMIZE  = -O2 
      ifeq ($(_FC),gfortran)
       ifeq ($(_CPU),aarch64)
         DONTHAVEM64OPT=Y
       endif
       ifneq ($(DONTHAVEM64OPT),Y)
         FOPTIONS   = -m64
         COPTIONS   = -m64
       endif
        COPTIONS += -Wall
        FOPTIONS   += -ffast-math #-Wunused  
        FOPTIMIZE  += -ffast-math -Wuninitialized
        DEFINES  += -DGFORTRAN
        DEFINES  += -DCHKUNDFLW -DGCC4
        GNUMAJOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __VERS | cut -c22)
        ifdef GNUMAJOR
        GNUMINOR=$(shell $(FC) -dM -E - < /dev/null 2> /dev/null | egrep __VERS | cut -c24)
        GNU_GE_4_6 = $(shell [ $(GNUMAJOR) -gt 4 -o \( $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 6 \) ] && echo true)
        GNU_GE_4_8 = $(shell [ $(GNUMAJOR) -gt 4 -o \( $(GNUMAJOR) -eq 4 -a $(GNUMINOR) -ge 8 \) ] && echo true)
        endif
        ifeq ($(GNU_GE_4_6),true)
          DEFINES  += -DGCC46
        endif
        ifeq ($(GNU_GE_4_8),true)
          FDEBUG =-O2 -g -fno-aggressive-loop-optimizations
          FOPTIMIZE +=-fno-aggressive-loop-optimizations
          FFLAGS_FORGA += -fno-aggressive-loop-optimizations
          FOPTIONS += -Warray-bounds
	  else
          FOPTIONS   += -Wuninitialized # -Wextra -Wunused
        endif
        ifdef USE_OPENMP
           FOPTIONS  += -fopenmp
           LDOPTIONS += -fopenmp
           DEFINES += -DUSE_OPENMP
        endif
      endif
      ifeq ($(_FC),gfortran)
        ifdef USE_I4FLAGS
#does not exists
#             FOPTIONS += -fdefault-integer-4
        else
             FOPTIONS += -fdefault-integer-8
        endif
      else ifeq ($(_FC),crayftn)
        ifdef USE_I4FLAGS
             FOPTIONS += -s integer32
        else
             FOPTIONS += -s integer64
        endif
      else
        ifdef USE_I4FLAGS
             FOPTIONS += -i4
        else
             FOPTIONS += -i8
        endif
      endif
      DEFINES  += -DEXT_INT
      MAKEFLAGS = -j 1 --no-print-directory
     ifeq ($(BLAS_LIB),)
       CORE_SUBDIRS_EXTRA += blas
     endif
     ifeq ($(LAPACK_LIB),)
       CORE_SUBDIRS_EXTRA += lapack
     endif
     RANLIB = echo
     DEFINES   +=   -DLINUX -DLINUX64
     ifeq ($(_CPU),alpha)
# using COMPAQ/DEC compilers (EA 3/13/2000)
       FC  = fort
       CC  = ccc      
       FOPTIONS   += -assume no2underscore -align dcommons -check nooverflow -assume accuracy_sensitive -check nopower -check nounderflow  -noautomatic
       DEFINES   +=   -DLINUXALPHA
       FOPTIMIZE =  -O4  -tune host -arch host  -math_library fast
       FVECTORIZE = -fast -O5 -tune host -arch host

       ifdef USE_I4FLAGS
         FOPTIONS +=  -fpe0
# needed: binutils 2.11 for -taso option with some more hacking on bfd/elf.c
       else
         FOPTIONS += -fpe3
       endif
         LINK.f = fort $(LDFLAGS)  
# this creates a static executable
#  LINK.f = fort $(LDFLAGS)   -Wl,-Bstatic
#         CORE_LIBS += -llapack $(BLASOPT) -lblas
     endif

    ifeq ($(_CPU),ia64)
# Itanium  
# g77 not working 
# i4 not working 
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
       _IFCV81= $(shell ifort -v  2>&1|egrep "Version "|head -n 1|awk ' /8\.1/  {print "Y";exit}; /9./ {print "Y"; exit}; /10./ {print "Y"; exit}')
       _IFCV8= $(shell ifort -v  2>&1|egrep "Version "|head -n 1|awk ' /8\./  {print "Y";exit}; /9./ {print "Y"; exit}; /10./ {print "Y"; exit}')
       ifeq ($(_IFCV8),Y)
         DEFINES+= -DIFCV8
#         FOPTIONS += -quiet
       endif	
       ifeq ($(_IFCV81),Y)
         DEFINES+= -DIFCV81
       endif	
        ITANIUMNO = $(shell   cat /proc/cpuinfo | egrep family | head -n 1  2>&1 | awk ' /Itanium 2/ { print "-tpp2"; exit };/Itanium/ { print "-tpp1"}')
        FOPTIONS   += -auto -w -ftz $(ITANIUMNO)
        ifdef  USE_GPROF
          FOPTIONS += -qp
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
#         FOPTIMIZE+= -IPF_fp_relaxed # breaks nwdft/xc/xc_pw91lda
       endif
        ifeq ($(_IFCV8),Y)
#         EXTRA_LIBS += -quiet
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

#     CORE_LIBS +=  $(BLASOPT) -llapack -lblas
endif # end of ia32 bit
    ifeq ($(_CPU),x86_64)
#
      MAKEFLAGS = -j 2 --no-print-directory
      COPTIMIZE = -O1
      ifeq ($(NWCHEM_TARGET),CYGWIN64)
        DEFINES += -DCYGWIN -DCYGNUS
      endif
ifeq ($(NWCHEM_TARGET),CATAMOUNT)
      FC=pgf90
      CC=gcc
endif

      ifeq ($(_FC),ifort)
     _GOTSSE3= $(shell cat /proc/cpuinfo | egrep sse3 | tail -n 1 | awk ' /sse3/  {print "Y"}')
       _IFCE = $(shell ifort -V  2>&1 |head -1 |awk ' /64/ {print "Y";exit};')
       _IFCV7= $(shell ifort -v  2>&1|egrep "Version "|head -n 1|awk ' /7./  {print "Y";exit}')
       _IFCV11= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 11) {print "Y";exit}}')
       _IFCV12= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 12) {print "Y";exit}}')
       _IFCV14= $(shell ifort -logo  2>&1|egrep "Version "|head -n 1|sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 14) {print "Y";exit}}')
       _IFCV15ORNEWER=$(shell ifort -logo  2>&1|egrep "Version "|head -n 1 | sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 15) {print "Y";exit}}')
# Intel EM64T is required
      ifneq ($(_IFCE),Y)
        defineFCE: 
	@echo
	@echo "   " ifort missing or not suitable x86_64 CPUs
	@echo
	@exit 1
      endif
       ifneq ($(_IFCV7),Y)
# to get EM64T
# Intel 8.1 is required
       else
           @echo ifort 8.1 is required for x86_64 CPUs
           @exit 1
       endif
       FDEBUG= -O2 -g
       FOPTIMIZE = -O3  -unroll  -ip
       FOPTIONS += -align
	   ifeq ($(_IFCV15ORNEWER), Y)
           FOPTIONS += -qopt-report-file=stderr
#fpp seems to get lost with ifort 15 in the offload bit
           EXPLICITF = TRUE
           CPP=fpp -P 
#           FOPTIONS +=  -Qoption,fpp,-P -Qoption,fpp,-c_com=no  -allow nofpp_comments 
          ifdef USE_OPTREPORT
	  FOPTIONS += -qopt-report=1 -qopt-report-phase=vec 
          endif
         ifdef USE_OPENMP
           FOPTIONS += -qopenmp
           FOPTIONS += -qopt-report-phase=openmp
           COPTIONS += -qopenmp
           DEFINES+= -DUSE_OPENMP 
         endif		   
	   else
         FOPTIONS += -vec-report6
         ifdef USE_OPENMP
           FOPTIONS += -openmp
           FOPTIONS += -openmp-report2
           COPTIONS += -openmp
           DEFINES+= -DUSE_OPENMP 
         endif
       endif
           
       ifdef USE_OFFLOAD
         ifeq ($(_IFCV14), Y)
          ### extra mic compile stuff; make FC=ifort CC=icc  AR=xiar
          FC = ifort
          _FC = ifort
          CC = icc
          _CC = icc
          AR = xiar
          EXTRA_LIBS += -loffload
          DEFINES+= -DUSE_OFFLOAD
          DEFINES+= -DINTEL_64ALIGN
          ifeq ($(_IFCV15ORNEWER), Y)
            FOPTIONS += -qopt-report-phase=offload
            FOPTIONS += -qoffload-option,mic,compiler,"-align array64byte"
            FOPTIONS += -align array64byte
	        LDOPTIONS += -qoffload-option,mic,compiler," -Wl,-zmuldefs"
            FOPTIONS += -watch=mic_cmd 
            COPTIONS += -qopt-report-phase=offload
		  else
            FOPTIONS += -opt-report-phase=offload
#           FOPTIONS += -offload-option,mic,compiler,"-mP2OPT_hlo_use_const_second_pref_dist=1"
            FOPTIONS += -offload-option,mic,compiler,"-align array64byte"
            FOPTIONS += -align array64byte
            FOPTIONS += -offload-option,mic,compiler,"-opt-report-phase=hlo"
	        LDOPTIONS += -offload-option,mic,compiler," -Wl,-zmuldefs"
            FOPTIONS += -watch=mic_cmd 
            COPTIONS += -opt-report-phase=offload
		  endif
         else
error100:
$(info     )
$(info USE_OFFLOAD requires ifort version 14 and later)
$(info     )
$(error )
         endif
       else
          ifdef USE_OPENMP
          ifeq ($(_IFCV15ORNEWER), Y)
             FOPTIONS += -qno-openmp-offload
          else
          ifeq ($(_IFCV14), Y)
             FOPTIONS += -no-openmp-offload
          endif
          endif
          endif
       endif
       DEFINES+= -DIFCV8 -DIFCLINUX
       ifeq ($(FC),ifc)
         FOPTIONS += -quiet
       endif
       ifdef  USE_FPE
         FOPTIONS += -fpe0 -traceback #-fp-model  precise
       endif
        ifeq ($(_IFCV11),Y) 
#next 2 lines needed for fp accuracy
        FDEBUG += -fp-model source
        ifeq ($(_IFCV12),Y) 
        FOPTIONS += -fimf-arch-consistency=true
        endif
        FOPTIMIZE += -xHost
       else
        ifeq ($(_GOTSSE3),Y) 
         FOPTIMIZE += -xP -no-prec-div
        else
         FOPTIMIZE +=  -tpp7 -ip 
         FOPTIMIZE += -xW
        endif
       endif
     endif	
#      
      ifeq ($(_FC),pgf90)
        FOPTIONS   += -Mdalign -Mllalign -Kieee 
#        FOPTIONS   += -tp k8-64  
#        FOPTIONS   +=    -Ktrap=fp
        FOPTIMIZE   = -O3 -fastsse -Mnounroll -Minfo=loop -Mipa=fast
        FVECTORIZE   = -fast  -fastsse  -O3   -Mipa=fast
        FDEBUG = -g -O2
        DEFINES  += -DCHKUNDFLW -DPGLINUX
        ifdef USE_OPENMP
           FOPTIONS  += -mp -Minfo=mp
           LDOPTIONS += -mp
           DEFINES += -DUSE_OPENMP
        endif
       ifeq ($(FC),ftn)
          LINK.f = ftn  $(LDFLAGS) $(FOPTIONS)
       endif
       ifeq ($(NWCHEM_TARGET),CATAMOUNT)
          LINK.f = ftn  $(LDFLAGS) $(FOPTIONS)
       endif
      endif
      ifeq ($(FC),pathf90)
#pathscale 1.3 compiler
# tested Sep 30 2004 on RH AW3
        FOPTIONS   += -cpp -Wp,-P
        FOPTIONS   += -fno-second-underscore -fixedform
        FOPTIONS   += -align64
#        FOPTIONS   += -LANG:heap_allocation_threshold=0
        FOPTIMIZE   = -O3 -OPT:Ofast:IEEE_arith=1:IEEE_NaN_inf=ON:Olimit=12000:ro=1:fold_reassociate=ON#:div_split=OFF:fast_nint=OFF
        FVECTORIZE  = -O3 -OPT:Ofast:ro=1 -fno-math-errno
        DEFINES  += -DCHKUNDFLW -DPSCALE
        FDEBUG = -g -O1
        LDOPTIONS = -Wl,--warn-once   -Wl,--relax
      endif
      ifeq ($(_CC),pgcc)
#        COPTIONS   =   -O
      endif
      ifeq ($(_CC),icc)
	 ICCV15ORNEWER=$(shell icc -V  2>&1|egrep "Version "|head -n 1 | sed 's/.*Version \([0-9][0-9]\).*/\1/' | awk '{if ($$1 >= 15) {print "Y";exit}}')
         COPTIONS   +=   -xHOST -ftz
         ifeq ($(ICCV15ORNEWER), Y)
   	    COPTIONS   += -qopt-report-phase=vec  -qopt-report-file=stderr
            ifdef USE_OPENMP
              COPTIONS += -qopenmp
   	      COPTIONS   +=  -qopt-report-phase:openmp
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
      ifeq ($(_CC),gcc)
        COPTIONS   =   -O3 -funroll-loops -ffast-math 
        ifdef USE_OPENMP
          COPTIONS += -fopenmp
        endif
      endif
      ifdef USE_GCC34
        COPTIONS  +=   -march=k8 -mtune=k8
      endif
#     CORE_LIBS +=  $(BLASOPT) -llapack -lblas
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

      ifeq ($(_FC),gfortran)
     _GOT3DNOW= $(shell cat /proc/cpuinfo | egrep 3dnowext | tail -n 1 | awk ' /3dnowext/  {print "Y"}')
#gcc version 4.1.0 20050525 (experimental)
       ifdef  USE_GPROF
          FOPTIONS += -pg
          COPTIONS += -pg
          LDOPTIONS += -pg
          LDFLAGS += -pg
        endif
	    LINK.f = $(FC)  $(LDFLAGS) 
        FOPTIMIZE  += -O3 
        FOPTIMIZE  += -mfpmath=sse -ffast-math
        FOPTIMIZE  += -fprefetch-loop-arrays #-ftree-loop-linear
        ifeq ($(GNU_GE_4_8),true)
          FOPTIMIZE  += -ftree-vectorize   -fopt-info-vec
        endif

        FDEBUG += -g -O 
        ifdef USE_F2C
#possible segv with use of zdotc (e.g. with GOTO BLAS)
#http://gcc.gnu.org/bugzilla/show_bug.cgi?id=20178
          FOPTIONS +=  -ff2c -fno-second-underscore
        endif
        ifeq ($(GNU_GE_4_6),true) 
          FOPTIMIZE += -march=native -mtune=native
        else
        ifeq ($(_GOT3DNOW),Y) 
#we guess its an opteron
          FOPTIMIZE += -march=opteron -mtune=opteron
        else
#we guess its a nocona em64t
          FOPTIMIZE += -march=nocona -mtune=nocona
        endif
        endif
#        FVECTORIZE  += -ftree-vectorize -ftree-vectorizer-verbose=1
      endif
      ifeq ($(_FC),crayftn)
	# Jeff: Cray Fortran supports preprocessing as of version 8.2.2 (at least)
        #EXPLICITF = FALSE
        #CPP = /usr/bin/cpp  -P -C -traditional
 	#CPPFLAGS += -DCRAYFORTRAN
        #FCONVERT = $(CPP) $(CPPFLAGS) $< > $*.f
        # USE_POSIXF is required because getlog is provided (GNU extension)
        FOPTIONS   +=  -Ktrap=fp -DCRAYFORTRAN -DUSE_POSIXF
        FDEBUG   =    -g
#       FOPTIMIZE = -O2 -O scalar3,thread0,vector1,ipa0
        FOPTIMIZE = -O2 -O scalar3,thread0,vector2,ipa2 #-rdm
      endif
      ifeq ($(_FC),craycc)
        COPTIONS   =   -O
      endif
endif

    ifeq ($(_CPU),ppc64)
# Tested on Red Hat Enterprise Linux AS release 3 (Taroon Update 3)
# Tested on SLES 9
# Feb 5th 2005
# xlf v9.1
# xlc v7.0 
# gcc-3.2.3-42 

#gfortran become default      FC=xlf
         ifeq ($(FC),xlf)
           _FC=xlf
         endif
      CC=gcc
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
        ifdef  USE_GPROF
          FOPTIONS += -pg
          LDOPTIONS += -pg
        endif
        FOPTIMIZE= -O3 -qstrict -qarch=auto -qtune=auto -qcache=auto -qfloat=fltint 
        FDEBUG= -O2 -g
        EXPLICITF = TRUE
        DEFINES  +=   -DXLFLINUX -DCHKUNDFLW
        CPP=/usr/bin/cpp  -P -C -traditional
        ifdef USE_I4FLAGS
          FOPTIONS += -qintsize=4
        else
          FOPTIONS += -qintsize=8
        endif
      endif
#     CORE_LIBS +=  $(BLASOPT) -llapack -lblas
#     EXTRA_LIBS +=  -dynamic-linker /lib64/ld64.so.1 -melf64ppc -lxlf90_r -lxlopt -lxlomp_ser -lxl -lxlfmath -ldl -lm -lc -lgcc -lm
    endif

     ifeq ($(BUILDING_PYTHON),python)
#   EXTRA_LIBS += -ltk -ltcl -L/usr/X11R6/lib -lX11 -ldl
   EXTRA_LIBS += -lnwcutil $(shell $(PYTHONHOME)/bin/python-config --libs) -lz
  LDOPTIONS += -Wl,--export-dynamic 
     endif
ifeq ($(NWCHEM_TARGET),CATAMOUNT)
        DEFINES  += -DCATAMOUNT
endif

endif
endif
#endof of LINUX64

ifeq ($(TARGET),FUJITSU_VPP)
#
# FUJITSU VPP5000 32-bit
#

         FC = frt
      CPP = /lib/cpp -P -C
     RANLIB = echo
  MAKEFLAGS = 
    INSTALL = @echo $@ is built
                        
    DEFINES = -DFUJITSU_VPP
    USE_MPI = TRUE

 #search include files for tools directories that are not built
   LIB_INCLUDES += -I$(NWCHEM_TOP)/src/tools/include \
                   -I$(NWCHEM_TOP)/src/tools/ma \
                   -I$(NWCHEM_TOP)/src/tools/tcgmsg-mpi \
                   -I$(NWCHEM_TOP)/src/tools/global/src \
                   -I$(NWCHEM_TOP)/src/tools/pario/eaf \
                   -I$(NWCHEM_TOP)/src/tools/pario/elio \
                   -I$(NWCHEM_TOP)/src/tools/pario/sf \
                   -I$(NWCHEM_TOP)/src/tools/pario/dra

 #change DEFINES and LIB_DEFINES so that frt understands them and add them
 #to FOPTIONS
     comma:= ,
     end:=
     space:= $(end) $(end)
   FDEFINES_1:= $(DEFINES) $(LIB_DEFINES)
   FDEFINES:= -Wp,$(subst $(space),$(comma),$(strip $(FDEFINES_1)))
    FOPTIONS = -w -Sw -KA32 $(FDEFINES)
    COPTIONS = -KA32
     FDEBUG = -Ob -g
  FOPTIMIZE = -Kfast
  COPTIMIZE = -K4

# removed global, ma, tcgmsg-mpi, as they are part of the native GA
 NW_CORE_SUBDIRS = include basis geom inp input  \
       pstat rtdb task symmetry util peigs $(CORE_SUBDIRS_EXTRA)

ifdef OLD_GA
        CORE_LIBS = -lnwcutil \
                    -L$(GA_LIBDIR) -lglobal -lpario -lma -lpeigs \
                    -ltcgmsg-mpi -L/usr/lang/mpi2/lib32 -lmpi -lmp
else
        CORE_LIBS = -lnwcutil \
                    -L$(GA_LIBDIR) -lga -lpeigs \
                    -L/usr/lang/mpi2/lib32 -lmpi -lmp
endif
       EXTRA_LIBS = -llapackvp -lblasvp -lsocket -Wl,-J,-P,-t,-dy
#end of FUJITSU_VPP 
endif

ifeq ($(TARGET),FUJITSU_VPP64)
#
# FUJITSU VX/VPP 64-bit
#

         FC = frt
      CPP = /lib/cpp -P -C
     RANLIB = echo
  MAKEFLAGS = 
    INSTALL = @echo $@ is built
                        
    DEFINES = -DFUJITSU_VPP -DEXT_INT
    USE_MPI = TRUE
     CORE_SUBDIRS_EXTRA = blas lapack

 #search include files for tools directories that are not built
   LIB_INCLUDES += -I$(NWCHEM_TOP)/src/tools/include \
                   -I$(NWCHEM_TOP)/src/tools/ma \
                   -I$(NWCHEM_TOP)/src/tools/tcgmsg-mpi \
                   -I$(NWCHEM_TOP)/src/tools/global/src \
                   -I$(NWCHEM_TOP)/src/tools/pario/eaf \
                   -I$(NWCHEM_TOP)/src/tools/pario/elio \
                   -I$(NWCHEM_TOP)/src/tools/pario/sf \
                   -I$(NWCHEM_TOP)/src/tools/pario/dra

 #change DEFINES and LIB_DEFINES so that frt understands them and add them
 #to FOPTIONS
     comma:= ,
     end:=
     space:= $(end) $(end)
   FDEFINES_1:= $(DEFINES) $(LIB_DEFINES)
   FDEFINES:= -Wp,$(subst $(space),$(comma),$(strip $(FDEFINES_1)))
    FOPTIONS = -w -Sw -KA64 -CcdII8 -CcdLL8 $(FDEFINES)
    COPTIONS = -KA64
     FDEBUG = -Ob -g
  FOPTIMIZE = -Kfast
  COPTIMIZE = -K4

# removed global, ma, tcgmsg-mpi, as they are part of the native GA
 NW_CORE_SUBDIRS = include basis geom inp input  \
       pstat rtdb task symmetry util peigs $(CORE_SUBDIRS_EXTRA)

ifdef OLD_GA
        CORE_LIBS = -lnwcutil \
                    -L$(GA_LIBDIR) -lglobal -lpeigs -lpario -lma \
                    -ltcgmsg-mpi -L/usr/lang/mpi2/lib64 -lmpi -lmp
else
        CORE_LIBS = -lnwcutil \
                    -L$(GA_LIBDIR) -lga -lpeigs \
                    -L/usr/lang/mpi2/lib64 -lmpi -lmp
endif
       EXTRA_LIBS = -llapack -lblas -lsocket -Wl,-J,-P,-t,-dy
#end of FUJITSU_VPP64
endif

ifeq ($(TARGET),cray-sv2)
#
# Cray sv2 aka x1
#

         FC = ftn
     RANLIB = echo
  MAKEFLAGS = 
    INSTALL = @echo $@ is built
    DEFINES =  -DEXT_INT  -DUSE_POSIXF  -DUSE_FFIO
    USE_FFIO = y
     CORE_SUBDIRS_EXTRA = blas lapack

   FOPTIONS =  -F -s integer64
   ifdef USE_SSP
      FOPTIONS += -O ssp
      COPTIONS += -h ssp
   endif
   FOPTIMIZE = -O scalar3,aggress,unroll2,vector2
      FDEBUG = -O scalar1,vector1
   COPTIMIZE = -O -h inline2  -h aggress


       EXTRA_LIBS = -lsci64 -llapack  -lblas  # need make dbl_to_sngl for this
#       EXTRA_LIBS =  -llapack  -lblas 

#      EXPLICITF     = TRUE
      FCONVERT      = $(CPP) $(CPPFLAGS)  $< | sed '/^\#/D'  > $*.f
#end of sv2
endif

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
    FOPTIONS   = -qEXTNAME -qxlf77=leadzero -NQ40000 -NT80000 -NS2048 -qmaxmem=8192
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
    FOPTIONS = -qEXTNAME -qxlf77=leadzero -NQ40000 -NT80000 -NS2048 -qmaxmem=8192
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

    ifeq ($(FC),mpif77)
        CC         = mpicc
        DEFINES   += -DGFORTRAN -DGCC4

        FOPTIONS  += -g -funderscoring -Wuninitialized 
        FOPTIMIZE += -O3 -ffast-math
        FDEBUG    += -O1 -g

        # EXT_INT means 64-bit integers are used
        DEFINES   += -DEXT_INT 
ifndef USE_I4FLAGS
        FOPTIONS  += -fdefault-integer-8
endif

        # linking ESSL is painful with gfortran
        CORE_LIBS +=  -llapack $(BLASOPT) -lblas

        # Here is an example for ALCF:
        # IBMCMP_ROOT=${IBM_MAIN_DIR}
        # BLAS_LIB=/soft/libraries/alcf/current/xl/BLAS/lib
        # LAPACK_LIB=/soft/libraries/alcf/current/xl/LAPACK/lib
        # ESSL_LIB=/soft/libraries/essl/current/essl/5.1/lib64
        # XLF_LIB=${IBMCMP_ROOT}/xlf/bg/14.1/bglib64
        # XLSMP_LIB=${IBMCMP_ROOT}/xlsmp/bg/3.1/bglib64
        # XLMASS_LIB=${IBMCMP_ROOT}/xlmass/bg/7.3/bglib64
        # MATH_LIBS="-L${XLMASS_LIB} -lmass -L${LAPACK_LIB} -llapack \
                     -L${ESSL_LIB} -lesslsmpbg -L${XLF_LIB} -lxlf90_r \
                     -L${XLSMP_LIB} -lxlsmp -lxlopt -lxlfmath -lxl \
                     -Wl,--allow-multiple-definition"
        # Note that ESSL _requires_ USE_64TO32 on Blue Gene
    endif

    ifeq ($(FC),mpixlf77_r)
        EXPLICITF  = TRUE
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
#        FOPTIMIZE += -g -O3 -qarch=qp -qtune=qp -qcache=auto -qunroll=auto -qfloat=rsqrt
        FOPTIMIZE += -O3 -qarch=qp -qtune=qp -qsimd=auto -qhot=level=1 -qprefetch -qunroll=yes #-qnoipa
        FOPTIMIZE += -qreport -qsource -qlistopt -qlist # verbose compiler output

        # ESSL dependencies should be provided by XLF linker
        CORE_LIBS +=  -llapack $(BLASOPT) -lblas
    endif

   endif # end BGQ


endif


###################################################################
#  All machine dependent sections should be above here, otherwise #
#  some of the definitions below will be 'lost'                   #
###################################################################
ifeq ($(BUILDING_PYTHON),python)
ifndef PYTHONHOME
errorpython1:
	@echo "For python you must define both PYTHONHOME and PYTHONVERSION"
	@echo "E.g., setenv PYTHONHOME /msrc/home/d3g681/Python-1.5.1"
	@echo "      setenv PYTHONVERSION 1.5"
	@echo " building_python <$(BUILDING_PYTHON)>"
	@echo " subdirs <$(NWSUBDIRS)>"
	@exit 1
endif
ifndef PYTHONVERSION
errorpython2:
	@echo "For python you must define both PYTHONHOME and PYTHONVERSION"
	@echo "E.g., setenv PYTHONHOME /msrc/home/d3g681/Python-1.5.1"
	@echo "      setenv PYTHONVERSION 1.5"
	@echo " building_python <$(BUILDING_PYTHON)>"
	@echo " subdirs <$(NWSUBDIRS)>"
	@exit 1
endif
#
ifndef PYTHONLIBTYPE
	PYTHONLIBTYPE=a
endif
ifndef PYTHONCONFIGDIR
	PYTHONCONFIGDIR=config
endif
ifdef USE_PYTHON64
	CORE_LIBS += $(PYTHONHOME)/lib64/python$(PYTHONVERSION)/$(PYTHONCONFIGDIR)/libpython$(PYTHONVERSION).$(PYTHONLIBTYPE)
else
  ifeq ($(GOTMINGW32),1)
  CORE_LIBS += $(PYTHONHOME)/libs/libpython$(PYTHONVERSION).$(PYTHONLIBTYPE)
  else
  CORE_LIBS += $(PYTHONHOME)/lib/python$(PYTHONVERSION)/$(PYTHONCONFIGDIR)/libpython$(PYTHONVERSION).$(PYTHONLIBTYPE)
  endif
endif
endif
#
######
#PAPI#
######
   
ifdef USE_PAPI 
  DEFINES += -DUSE_PAPI
  ifndef PAPI_LIB
    papierror1:
	echo You must define PAPI_LIB in your environment to be the path
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

#
ifdef USE_FDIST
  DEFINES += -DFDIST
endif

_USE_SCALAPACK = $(shell cat ${NWCHEM_TOP}/src/tools/build/config.h | awk ' /HAVE_SCALAPACK\ 1/ {print "Y"}')

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
  CORE_LIBS += $(SCALAPACK) $(PBLAS) $(BLACS)
endif

ifdef USE_64TO32
      CORE_LIBS +=  -l64to32
NWSUBDIRS += 64to32blas
endif
ifdef BLASOPT
       CORE_LIBS +=  $(BLASOPT) 
endif
ifeq ($(LAPACK_LIB),)
      CORE_LIBS +=  -llapack 
endif
ifeq ($(BLAS_LIB),)
      CORE_LIBS +=  -lblas 
endif

ifdef USE_NOIO
 DEFINES += -DNOIO
endif


ifdef USE_SUBGROUPS
  DEFINES += -DGANXTVAL -DUSE_SUBGROUPS
#turn off peigs for now
else
  DEFINES += -DPARALLEL_DIAG
endif
###################################################################
#  All machine dependent sections should be above here, otherwise #
#  some of the definitions below will be 'lost'                   #
###################################################################
#the new GA uses ARMCI library
ifdef OLD_GA
      CORE_LIBS += -larmci
else
  ifeq ($(ARMCI_NETWORK),ARMCI)
    ifdef EXTERNAL_ARMCI_PATH
      CORE_LIBS += -L$(EXTERNAL_ARMCI_PATH)/lib -larmci
    else
      CORE_LIBS += -L$(NWCHEM_TOP)/src/tools/install/lib -larmci
    endif
  else
      CORE_LIBS +=
  endif
endif

# MPI version requires tcgmsg-mpi library

ifdef USE_MPI 
ifeq ($(FC),$(findstring $(FC),mpifrt mpfort mpif77 mpxlf mpif90 ftn))
  LIBMPI =
  MPI_INCLUDE =
  MPI_LIB =
else
ifndef MPI_INCLUDE
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

# CUDA
#ckbn gpu
ifdef TCE_CUDA
 CORE_LIBS += $(CUDA_LIBS)
endif

# lower level libs used by communication libraries 
 
COMM_LIBS=  $(shell grep ARMCI_NETWORK_LIBS\ = ${NWCHEM_TOP}/src/tools/build/Makefile | cut -b 22-)
COMM_LIBS +=  $(shell grep ARMCI_NETWORK_LDFLAGS\ = ${NWCHEM_TOP}/src/tools/build/Makefile | cut -b 24-)

ifdef COMM_LIBS 
 CORE_LIBS += $(COMM_LIBS) 
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
else
# Need FDEBUG after FOPTIONS on SOLARIS to correctly override optimization
    FFLAGS = $(FOPTIONS) $(FDEBUG) 
    CFLAGS = $(COPTIONS) $(CDEBUG) 
endif
  INCLUDES = -I. $(LIB_INCLUDES) -I$(INCDIR) $(INCPATH)
  CPPFLAGS = $(INCLUDES) $(DEFINES) $(LIB_DEFINES)
   LDFLAGS = $(LDOPTIONS) -L$(LIBDIR) $(LIBPATH)
      LIBS = $(NW_MODULE_LIBS) $(CORE_LIBS) 

# I think this will work everywhere, but it might have to become
# machine-dependent 

MKDIR = mkdir

#
# Define known suffixes mostly so that .p files don't cause pc to be invoked
#

.SUFFIXES:	
.SUFFIXES:	.o .s .F .f .c .cpp

ifndef FLINT

ifdef EXPLICITF
#
# Needed on machines where FC does not preprocess .F files
# with CPP to get .f files
#
# These rules apply to make-ing of files in specfic directories
.SUFFIXES:	
.SUFFIXES:	.o .s .F .f .c

.F.o:	
	@echo Converting $*.F '->' $*.f
	@$(FCONVERT)
	$(FC) -c $(FFLAGS) $*.f
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
(%.o):	%.F
ifdef EXPLICITF
	@echo Converting $< '->' $*.f
	@$(FCONVERT)
	$(FC) -c $(FFLAGS) $*.f
ifndef NWCHEM_KEEPF
	@/bin/rm -f $*.f
endif
else
ifeq ($(XLFMAC),y)
	$(FC) -c $(FFLAGS)     $<
else
	$(FC) -c $(FFLAGS) $(CPPFLAGS)  $<
endif
endif

(%.o):	%.f
	$(FC) -c $(FFLAGS) $<

(%.o):	%.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $% $<

(%.o):  %.cu
	$(CUDA) -c $(CUDA_FLAGS) -c  $(CUDA_LIBS) -o $% $<

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


