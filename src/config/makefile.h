#
# $Id: makefile.h,v 1.294 1999-07-21 18:42:33 d3e129 Exp $
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
# a hypehn to NWCHEM_TOP in order to derive the directory path TOPDIR.

# For development tree 
RELEASE := 
# For current release tree
#RELEASE := 3.1

#

ifndef NWCHEM_TOP
error1:
	@echo You must define NWCHEM_TOP in your environment to be the path
	@echo of the top level nwchem directory ... something like
	@echo     setenv NWCHEM_TOP /msrc/home/elvis/nwchem
	@exit 1
endif

# Select the old (pre-1999 version of GA) by uncommenting the next line
#OLD_GA = y 

#
# Do a setenv for NWCHEM_TARGET to be the machine and NWCHEM_TARGET_CPU the CPU to build for
#
# NWCHEM_TARGET :  CONVEX-SPP
#                  CRAY-T3D
#                  CRAY-T3E
#                  DECOSF
#                  DELTA
#                  IBM
#                  KSR
#                  LINUX        NWCHEM_TARGET_CPU :
#                                                  nothing for X86 (e.g. do not set this)
#                                                  ALPHA for AlphaLinux (broke)
#                                                  POWERPC for MkLinux (broke)
#                  PARAGON
#                  SGI
#                  SGI_N32      NWCHEM_TARGET_CPU : R8000 or R10000
#                  SGITFP       NWCHEM_TARGET_CPU : R8000 or R10000
#                  SOLARIS
#                  SP           NWCHEM_TARGET_CPU : P2SC
#                      (uses non-thread safe libraries and MPL)
#                  LAPI      NWCHEM_TARGET_CPU : P2SC
#                      (uses thread safe libraries and LAPI)
#                  SUN
#
# Note that the preprocessor flags for CRAY-T3D and CRAY-T3E are CRAY_T3D and CRAY_T3E respectively
#

ifndef NWCHEM_TARGET
error2:
	@echo You must define NWCHEM_TARGET in your environment to be the name
	@echo of the machine you wish to build for ... for example
	@echo     setenv NWCHEM_TARGET SUN
	@echo Known targets are SUN, DELTA, ...
	@exit 1
endif

#JN SP1 name is obsolete 
ifeq ($(NWCHEM_TARGET),SP1)
    NWCHEM_TARGET = SP
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
     LIBDIR := $(TOPDIR)/lib/$(TARGET)
     BINDIR := $(TOPDIR)/bin/$(TARGET)
     INCDIR := $(TOPDIR)/src/include
     CNFDIR := $(TOPDIR)/src/config

#
# Define LIBPATH to be paths for libraries that you are linking in
# from precompiled sources and are not building now. These libraries
# will be searched AFTER anything you are building now.
# e.g. LIBPATH = -L/msrc/proj/mss/lib
#
    LIBPATH = 
    LIBPATH = -L$(NWCHEM_TOP)/src/tools/lib/$(TARGET)

#
# Define INCPATH to be directories to get includes for
# libraries that you are not building now.  These directories
# will be searched AFTER anything you are building now.
#
    INCPATH = 
    INCPATH = -I$(NWCHEM_TOP)/src/tools/include

# These subdirectories will build the core, or supporting libraries
# that are required by all NWChem modules.  The include directory is
# first to insure that all include files will be properly installed
# prior to trying to compile anything.
#
# The core libraries are usually rather platform-dependent and are
# specified below.  Use of MPI requires substituting the tcgmsg-mpi
# wrapper for the normal tcgmsg library.
# MPLIB - represents the name of mpi library
MPILIB = -lmpi

#JN: under the new structure, tools should be listed first as
# their header files are needed for dependency analysis of
# other NWChem modules

NW_CORE_SUBDIRS = tools include basis geom inp input  \
	pstat rtdb task symmetry util peigs $(CORE_SUBDIRS_EXTRA)

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

# Only the SUN, SGI, KSR, PARAGON, DELTA and IBM versions are up to date

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
# SCRATCH_DEF_DIR = Site specific default directory for scratch files
# PERM__DEF_DIR   = Site specific default directory for permanent files
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
#           SP                           P2SC
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

ifndef SCRATCH_DEF_DIR
 SCRATCH_DEF_DIR = "'.'"
endif
ifndef PERM_DEF_DIR
 PERM_DEF_DIR   = "'.'"
endif

#
# Machine specific stuff
#

ifeq ($(TARGET),SUN)
#
# Sun running SunOS
#

       NICE = nice
      SHELL := $(NICE) /bin/sh
    CORE_SUBDIRS_EXTRA = blas lapack
         CC = gcc
     RANLIB = ranlib
  MAKEFLAGS = -j 1 --no-print-directory
    INSTALL = @echo $@ is built

   FOPTIONS = -Nl199 -fast -dalign
   COPTIONS = -Wall
# -O4 breaks at least inp_* and seems no faster than -O3
  FOPTIMIZE = -O3
  COPTIMIZE = -g -O2

    DEFINES = -DSUN

       CORE_LIBS =  -lutil -lchemio -lglobal -lma -lpeigs -llapack -lblas
endif

ifeq ($(TARGET),SOLARIS)
#
# Sun running Solaris 2.4 or later
# if you want to use purecoverage tool you must
#
# setenv PURECOV 1
#FLINT = 1
#
      SHELL := $(NICE) /bin/sh
    CORE_SUBDIRS_EXTRA = blas lapack
	CPP = /usr/ccs/lib/cpp
#        CC = gcc
#  COPTIONS = -Wall
# COPTIMIZE = -g -O2
         CC = cc
   COPTIONS = 
  COPTIMIZE = -fast
     RANLIB = echo
  MAKEFLAGS = -j 2 --no-print-directory
    INSTALL = echo $@ is built
# -fast introduces many options that must be applied to all files
# -stackvar puts locals on the stack which seems a good thing
#     but may need to increase the stacksize at runtime using limit
# -xs allows debugging without .o files
   FOPTIONS = -Nl199 -fast -dalign -stackvar
# Under Solaris -O3 is the default with -fast (was -O2 with SUNOS)
# -fsimple=2 enables more rearranging of floating point expressions
# -depend enables more loop restructuring ... now implicit in -fast?
# -xvector requires -mvec library
  FOPTIMIZE = -O3 -fsimple=2 -depend -xvector=yes 
# This for ultra-2 -xarch=v8plusa
# Under Solaris -g no longer disables optimization ... -O2 seems solid
# but is slow and impairs debug ... use -O1 for speed and debugability.
# -fast now turns on -depend so must turn it off
     FDEBUG = -g -O1 -nodepend
   LIBPATH += -L/usr/ucblib
    DEFINES = -DSOLARIS  -DNOAIO

# -DPARALLEL_DIAG

  LDOPTIONS = -xildoff

       CORE_LIBS = -lutil -lchemio -lglobal -lma -lpeigs -llapack -lblas
# First four needed for parallel stuff, last for linking with profiling
      EXTRA_LIBS = -lsocket -lrpcsvc -lnsl -lucb -lmvec -ldl 
ifeq ($(BUILDING_PYTHON),python)
# needed if python was compiled with gcc (common)
      EXTRA_LIBS += -L/msrc/apps/gcc-2.8.1/lib/gcc-lib/sparc-sun-solaris2.6/2.8.1/ -lgcc
# needed here if using a python version with tk/tcl extensions (common)
      EXTRA_LIBS += -L/msrc/apps/lib -ltk8.0 -ltcl8.0 
# needed here if using a python version built with BLT extensions
#     EXTRA_LIBS += -L/msrc/apps/lib -lBLT
# Both tk/tcl and BLT need X11 (common)
      EXTRA_LIBS += -lX11
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
   CORE_LIBS = -lutil -lglobal -lma -lpeigs -llapack -lblas
# First four needed for parallel stuff, last for linking with profiling
	   EXTRA_LIBS = -lsocket -lrpcsvc -lnsl -lucb -lintl -lc -lc -lpurecov_stubs

#end of puresolaris
endif


ifeq ($(TARGET),CRAY-T3D)
#
# CRAY-T3D cross-compiled on YMP (atw)
#
   CORE_SUBDIRS_EXTRA =	blas lapack # Only a couple of routines not in scilib
               RANLIB = echo
            MAKEFLAGS = -j 4 --no-print-directory
              INSTALL = @echo $@ is built
        OUTPUT_OPTION = 

                   FC = cf77 
                  CPP = cpp -P  -N
#                   FC = /mpp/bin/cf77 
#                  CPP = /mpp/lib/cpp -P  -N
# gpp does not eat elif
#                 CPP = /usr/lib/gpp -P  -F
# need jump since with all modules code is too big for branches
# noieeedivide seems safe and should be faster
             FOPTIONS = -dp -Ccray-t3d 
             COPTIONS = -Tcray-t3d -hjump
# To make executable smaller use scalar optimization and no -g on all code.
# symmetry/(dosymops.F,sym_movecs_apply_op) break with scalar
# (it is handled separately in symmetry/makefile)
# !! Note that -O option disables any -Wf"-o options" but we need
# !! jump so cannot use -O.
               FDEBUG = -Wf"-o scalar,jump,noieeedivide"
# Not sure yet if these are fully safe ... aggress,unroll
            FOPTIMIZE = -Wf"-o scalar,jump,noieeedivide,aggress,unroll"
               CDEBUG = -O 1
            COPTIMIZE = -O
# -s eliminates symbol tables to make executable smaller (remove -s
# if you want to debug) ... the T3 does load the symbol table (stupid!)
# No need for forcing of block data here as long as each one is
# referenced by an external statement in the code.
            LDOPTIONS = -s -Drdahead=on -L$(LIBDIR) 

# Compilation also depends on compilers defining CRAY
              DEFINES = -DCRAY_T3D -DPARALLEL_DIAG

#               LINK.f = /mpp/bin/mppldr $(LDOPTIONS)
               LINK.f = mppldr $(LDFLAGS)

            CORE_LIBS =  -lutil -lchemio -lglobal -lma -lpeigs -llapack -lblas 

      EXPLICITF     = TRUE
      FCONVERT      = $(CPP) $(CPPFLAGS)  $< | sed '/^\#/D'  > $*.f
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
            FOPTIMIZE = -O scalar3,aggress,unroll2,vector3
#,pipeline3
               CDEBUG = -O 1
            COPTIMIZE = -O
#
# to debug code you must remove the -s flag unless you know assembler
#
            LDOPTIONS = -Wl"-s" -Xm  -lmfastv

              DEFINES = -DCRAY_T3E -DCRAY_T3D -D__F90__ -DPARALLEL_DIAG

               LINK.f = f90 $(LDFLAGS)

            CORE_LIBS = -lutil -lchemio -lglobal -lma -lpeigs -llapack -lblas
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

ifeq ($(TARGET),KSR)
#
# KSR running OSF
#
# JN 96/10/02:
# Replaced -DLongInteger with -DEXT_INT for consistency with GA, DRA, PEIGS ...

	CPP = /usr/lib/cpp -P -C
    CORE_SUBDIRS_EXTRA = blas
     RANLIB = echo
  MAKEFLAGS = -j 10 --no-print-directory
    INSTALL = @echo $@ is built

   FOPTIONS = -r8
   COPTIONS = 
  FOPTIMIZE = -xfpu3 -O1
  COPTIMIZE = -xfpu3 -O1

    DEFINES = -DKSR -DPARALLEL_DIAG -DEXT_INT

#       LIBPATH += -L/home/d3g681/TCGMSG_DISTRIB
        LIBPATH += -L/home2/d3g270/peigs1.1.1 -L/home/d3g681/TCGMSG_DISTRIB
       CORE_LIBS = -lglobal -lma -lutil -lchemio -lpeigs \
                   -lksrlapk -lksrblas -llapack2 -lblas2  
      EXTRA_LIBS = -para -lrpc
endif

ifeq ($(TARGET),PARAGON)
#
# Intel Paragon running OSF
# (native build, or cross-compiled -- the latter is faster in most cases)
#
    CORE_SUBDIRS_EXTRA = lapack

         FC = if77
         CC = icc
         AR = ar860
     RANLIB = echo
	CPP = /usr/lib/cpp -P -C

  MAKEFLAGS = -j 4  --no-print-directory
    INSTALL = @echo $@ is built

  FOPTIONS = -Knoieee
  COPTIONS = -Knoieee
 FOPTIMIZE = -O2 -Minline=1000
FVECTORIZE = -O2 -Minline=1000 # -Mvect
 COPTIMIZE = -O2

#
# __PARAGON__ is defined by the PGI cpp, but NOT when it is invoked by
# the f77 comiler!!!
#
# Do NOT define -DNX or -DIPSC for paragon -- they get into some code for
# real iPSC and Delta that is not applicable to the paragon, which is more
# unixy since it runs OSF/1.
#
    DEFINES = -D__PARAGON__ -DPARALLEL_DIAG
    ARFLAGS = ru

# CAUTION: PGI's linker thinks of -L as adding to the _beginning_ of the
# search path -- contrary to usual unix usage!!!!!
       LIBPATH  += -L/home/delilah11/gifann/lib
       CORE_LIBS = -lchemio -lglobal -lma -lutil -lpeigs \
	           -llapack $(LIBDIR)/liblapack.a -lkmath 
      EXTRA_LIBS = -nx
endif

ifeq ($(TARGET),DELTA)
#
# DELTA/IPSC running NX
#
    CORE_SUBDIRS_EXTRA = lapack blas # -lkamth not reliable
# blas
        FC = if77
        CC = icc
       CPP = /usr/lib/cpp
        AR = ar860
    RANLIB = echo

   INSTALL = "strip860 $(BINDIR)/nwchem; rcp $(BINDIR)/nwchem delta1:"
 MAKEFLAGS = -j 2  --no-print-directory

  FOPTIONS = -Knoieee
  COPTIONS = -Knoieee
 FOPTIMIZE = -O2 		# -Minline=1000 ## Inlining bombs for dtrtri.f
FVECTORIZE = -O4 		# -Mvect corrupts lapack for large vectors
 COPTIMIZE = -O2

   DEFINES = -DNX -DDELTA -DIPSC -DNO_BCOPY  -D__IPSC__ -DPARALLEL_DIAG
        LIBPATH += -L/home/delilah11/gifann/lib
       CORE_LIBS = -lglobal -lma -lutil -lchemio -lglobal -lpeigs \
		   $(LIBDIR)/liblapack.a -llapack -lblas
      EXTRA_LIBS = -node
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
# Added -lutil to core libraries
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
# Added -lutil to core libraries (again!)
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

	CPP = /usr/lib/cpp
  CORE_SUBDIRS_EXTRA = blas lapack
         FC = f77
     RANLIB = echo


    INSTALL = @echo nwchem is built
  MAKEFLAGS = -j 4 --no-print-directory

  FOPTIONS = -d8 -i8 -mips4 -align64 -64 -r8 -G 0 -OPT:roundoff=3:IEEE_arithmetic=3
  COPTIONS = -fullwarn -mips4 -64

#optimization flags for R8000 (IP21)
 FOPTIMIZE_8K = -O3 -OPT:fold_arith_limit=4000:const_copy_limit=20000:global_limit=20000:fprop_limit=2000 -TENV:X=3 -WK,-so=1,-o=1,-r=3,-dr=AKC
FVECTORIZE_8K = -O3 -OPT:fold_arith_limit=4000 -TENV:X=3 -WK,-dr=AKC

#optimization flags for R10000 (IP28)
 FOPTIMIZE_10K = -O3 -OPT:fold_arith_limit=4000:const_copy_limit=20000:global_limit=20000:Olimit=4000 -TENV:X=1 -WK,-so=1,-o=1,-r=3,-dr=AKC
FVECTORIZE_10K = -O3 -OPT:fold_arith_limit=4000 -TENV:X=1 -WK,-dr=AKC

 COPTIMIZE = -O
 FOPTIMIZE = -O3

ifeq ($(NWCHEM_TARGET_CPU),R10000)
 FOPTIMIZE = $(FOPTIMIZE_10K)
 FVECTORIZE = $(FVECTORIZE_10K)
endif
ifeq ($(NWCHEM_TARGET_CPU),R8000)
 FOPTIMIZE = $(FOPTIMIZE_8K)
 FVECTORIZE = $(FVECTORIZE_8K)
endif

    DEFINES = -DSGI -DSGITFP -DEXT_INT -DPARALLEL_DIAG
  CORE_LIBS = -lutil -lchemio -lglobal -lma -lpeigs -llapack -lblas
endif


ifeq ($(TARGET),SGI)
#
# SGI normal (32-bit platform)
#
# CORE_SUBDIRS_EXTRA are those machine specific libraries required 
#
# JN 12/4/96: 
# removed -lblas from CORE_SUBDIRS_EXTRA and -lmalloc from the EXTRA_LIBS
# replaced -mips2 with -mips3
# On IRIX >= 6.1, SGI recomends using -n32 to utilize internal 64-bit
# registers (up to 50% better floating point performance over -32 flag) 

    CORE_SUBDIRS_EXTRA = lapack
         FC = f77
         AR = ar
     RANLIB = echo
     CPP = /usr/lib/cpp

    INSTALL = @echo nwchem is built
  MAKEFLAGS = -j 4 --no-print-directory
    DEFINES = -DSGI

  FOPTIONS = -32 -Nn10000 # -mips3
  COPTIONS =  -32 -fullwarn #-mips3
 FOPTIMIZE = -O2
 COPTIMIZE = -O2

       CORE_LIBS = -lutil -lchemio -lglobal -lma -lpeigs -llapack -lblas
#     EXTRA_LIBS = -lmalloc 
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

  FOPTIONS = -n32 -mips4 -G 0 -OPT:roundoff=3:IEEE_arithmetic=3
  COPTIONS = -n32 -mips4 -fullwarn

#optimization flags for R8000 (IP21)
 FOPTIMIZE_8K = -O3 -OPT:fold_arith_limit=4000:const_copy_limit=20000:global_limit=20000:fprop_limit=2000 -TENV:X=3 -WK,-so=1,-o=1,-r=3,-dr=AKC
FVECTORIZE_8K = -O3 -OPT:fold_arith_limit=4000 -TENV:X=3 -WK,-dr=AKC

#optimization flags for R10000 (IP28)
 FOPTIMIZE_10K = -O3 -OPT:fold_arith_limit=4000:const_copy_limit=20000:global_limit=20000:fprop_limit=2000 -TENV:X=1 -WK,-so=1,-o=1,-r=3,-dr=AKC -SWP:if_conversion=OFF
FVECTORIZE_10K = -O3 -OPT:fold_arith_limit=4000 -TENV:X=1 -WK,-dr=AKC -SWP:if_conversion=OFF

 FOPTIMIZE = -O3
 COPTIMIZE = -O2

ifeq ($(NWCHEM_TARGET_CPU),R10000)
 FOPTIMIZE = $(FOPTIMIZE_10K)
 FVECTORIZE = $(FVECTORIZE_10K)
endif
ifeq ($(NWCHEM_TARGET_CPU),R8000)
 FOPTIMIZE = $(FOPTIMIZE_8K)
 FVECTORIZE = $(FVECTORIZE_8K)
endif
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

       CORE_LIBS = -lutil -lchemio -lglobal -lma -lpeigs -llapack -lblas
endif

ifeq ($(TARGET),HPUX)
#
# HPUX 10.1
#

  CORE_SUBDIRS_EXTRA = blas lapack
  CPP = /lib/cpp -P
# CC = gcc
  CC = cc
  FC = f77
  LDOPTIONS = -g -L/usr/lib
  LINK.f = fort77   $(LDFLAGS)
  CORE_LIBS = -lutil -lchemio -lglobal -lma -lpeigs -ltcgmsg -llapack -lblas -lU77 -lM -lm
  CDEBUG =
  FDEBUG = -g
  FOPTIONS =  +ppu
# COPTIONS =
  COPTIONS = -Aa -D_HPUX_SOURCE +e
  FOPTIMIZE = -g
  FVECTORIZE = -O +Oaggressive
# COPTIMIZE = -fthread-jumps -fdefer-pop -fdelayed-branch -fomit-frame-pointer -finline-functions -ffast-math
  COPTIMIZE = -O
  RANLIB = echo

  DEFINES = -DHPUX -DEXTNAME -DPARALLEL_DIAG

endif


ifeq ($(TARGET),CONVEX-SPP)
#
# Convex SPP-1200 running SPP-UX 3.2
#

        CPP = /lib/cpp -P
         FC = fc

# -g is not recognized; 
# Convex debug flag -cxdb does not disable optimization 
     CDEBUG = -no
     FDEBUG = -no
   FOPTIONS = -ppu -or none
   COPTIONS = -or none
  FOPTIMIZE = -O1
  COPTIMIZE = -O

    DEFINES = -DCONVEX -DHPUX -DEXTNAME

# &%@~* Convex compiler will preprocess only *.f and *.FORT files !
  EXPLICITF = TRUE
   FCONVERT = $(CPP) $(CPPFLAGS) $< > $*.f
endif



ifeq ($(TARGET),IBM)
#
# IBM AIX
#

    CORE_SUBDIRS_EXTRA = lapack
#blas
         FC = xlf
         CC = xlc
    ARFLAGS = urs
     RANLIB = echo
  MAKEFLAGS = -j 1 --no-print-directory
    INSTALL = @echo $@ is built
        CPP = /usr/lib/cpp -P

   FOPTIONS = -qEXTNAME -qnosave # -qalign=4k 
# -qinitauto=FF
   COPTIONS = 
# -qstrict required with -O3 (according to Edo)
# -qfloat=rsqrt gives faster square roots (off by -qstrict)
# -qfloat=fltint gives faster real-integer conversion (off by -qstrict)
# -qhot seems to break a lot of things so don't ever use it
# -qarch=pwr (for peril) com (for any) , pwr2  or ppc
  FOPTIMIZE = -O3 -qstrict -qfloat=rsqrt:fltint -NQ40000 -NT80000
  COPTIMIZE = -O

    DEFINES = -DIBM -DAIX -DEXTNAME
ifdef USE_ESSL
   DEFINES += -DESSL
endif

       LIBPATH += -L/usr/lib -L/msrc/apps/lib

       CORE_LIBS = -lchemio -lglobal -lma -lutil -lpeigs -llapack -lblas \
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
              -brename:.dsyr2_,.dsyr2

#              -brename:.dtrsv_,.dtrsv \
#              -brename:.ztrsv_,.ztrsv 
#              -brename:.dsymm_,.dsymm \

#              -brename:.dznrm2_,.dznrm2 \
#              -brename:.zaxpy_,.zaxpy \
#              -brename:.zcopy_,.zcopy \
#              -brename:.zdotc_,.zdotc \
#              -brename:.zdscal_,.zdscal \

#              -brename:.zgemv_,.zgemv \
#              -brename:.zgerc_,.zgerc \
#              -brename:.zhemm_,.zhemm \
#              -brename:.zhemv_,.zhemv \
#              -brename:.zher2_,.zher2 \
#              -brename:.zher2k_,.zher2k \
#              -brename:.zherk_,.zherk \
#              -brename:.zscal_,.zscal \
#              -brename:.zswap_,.zswap \
#              -brename:.ztrmm_,.ztrmm \
#              -brename:.ztrmv_,.ztrmv \
#              -brename:.ztrsm_,.ztrsm \
#

##comment out from dtrmm_ inclusive
#ifdef USE_ESSL
#       CORE_LIBS += -lessl
#endif

  EXPLICITF = TRUE
  FCONVERT = $(CPP) $(CPPFLAGS) $< > $*.f
#
endif


ifeq ($(TARGET),SP)
#
     OLD_GA = y 
    CORE_SUBDIRS_EXTRA = lapack blas
         FC = mpxlf -qnohpf
# -F/u/d3g681/xlhpf.cfg:rjhxlf
         CC = mpcc
    ARFLAGS = urs
     RANLIB = echo
  MAKEFLAGS = -j 1 --no-print-directory
    INSTALL = @echo $@ is built
        CPP = /usr/lib/cpp -P

LARGE_FILES = YES

  LDOPTIONS = -lc -lm -qEXTNAME -qnosave -g -bmaxdata:0x20000000 -bloadmap:nwchem_map
ifeq ($(NWCHEM_TARGET_CPU),604)
  LDOPTIONS = -lxlf90 -lm -qEXTNAME -qnosave -g -bmaxdata:0x20000000 -bloadmap:nwchem_map
endif

   LINK.f   = mpxlf -qnohpf $(LDFLAGS)
   FOPTIONS = -qEXTNAME -qnosave
# -qinitauto=7F # note that grad_force breaks with this option
   COPTIONS = 
  FOPTIMIZE = -O3 -qstrict -qfloat=rsqrt:fltint -NQ40000 -NT80000
  COPTIMIZE = -O
ifeq ($(NWCHEM_TARGET_CPU),604)
  FOPTIMIZE += -qarch=604
  COPTIMIZE += -qarch=ppc
endif

ifeq ($(NWCHEM_TARGET_CPU),P2SC)
# These from George from Kent Winchell
  FOPTIMIZE += -qcache=type=d:level=1:size=128:line=256:assoc=4:cost=14 \
        -qcache=type=i:level=1:size=32:line=128
  COPTIMIZE += -qcache=type=d:level=1:size=128:line=256:assoc=4:cost=14 \
        -qcache=type=i:level=1:size=32:line=128
endif

    DEFINES = -DSP1 -DAIX -DEXTNAME -DPARALLEL_DIAG
#
# Prefix LIBPATH with -L/usr/lib for AIX 3.2.x
#
#  LIBPATH += -L/sphome/harrison/peigs2.0

  CORE_LIBS = -lchemio -lglobal -lma -lutil -lpeigs -llapack -lblas


   USE_ESSL = YES
#   USE_BLAS = YES
ifdef USE_ESSL
   DEFINES += -DESSL
# renames not needed for 4.1.  Still are for 3.2.
ifeq ($(NWCHEM_TARGET_CPU),P2SC)
 CORE_LIBS += -lesslp2 -lpesslp2 -lblacsp2
else
CORE_LIBS += -lessl -lpessl -lblacs
endif

#	      -brename:.daxpy_,.daxpy \
#	      -brename:.dgesv_,.dgesv \
#	      -brename:.dcopy_,.dcopy \
#	      -brename:.ddot_,.ddot \
#	      -brename:.dgemm_,.dgemm \
#	      -brename:.dgemv_,.dgemv \
#	      -brename:.dgetrf_,.dgetrf \
#	      -brename:.dgetrs_,.dgetrs \
#	      -brename:.dscal_,.dscal \
#	      -brename:.dspsvx_,.dspsvx \
#	      -brename:.dpotrf_,.dpotrf \
#	      -brename:.dpotri_,.dpotri \
#	      -brename:.idamax_,.idamax 
ifdef USE_BLAS
 CORE_LIBS += -lblas -brename:.xerbla_,.xerbla -brename:.lsame_,.lsame
endif

else
    CORE_SUBDIRS_EXTRA += blas
             CORE_LIBS += -lblas
endif


# IMPORTANT:  These renames are necessary if you try to link against
# a copy of PeIGS built for MPI instead of TCGMSG. (Not recommended, 
# see INSTALL)
# mpipriv is a common block used in MPICH's implementation of MPI.  It
# is critical that this common block is renamed correctly because
# the linker will not detect any problems (there will be separate
# common blocks labeled mpipriv_ and mpipriv) but the program will not
# operate correctly.
#ifdef USE_MPI
#   CORE_LIBS += -brename:.mpi_recv_,.mpi_recv \
#		-brename:.mpi_initialized_,.mpi_initialized \
#		-brename:.mpi_init_,.mpi_init \
#		-brename:.mpi_comm_rank_,.mpi_comm_rank \
#		-brename:.mpi_comm_size_,.mpi_comm_size \
#		-brename:.mpi_finalize_,.mpi_finalize \
#		-brename:.mpi_send_,.mpi_send \
#		-brename:mpipriv_,mpipriv
#endif

 EXPLICITF = TRUE
  FCONVERT = $(CPP) $(CPPFLAGS) $< > $*.f
#
endif


ifeq ($(TARGET),LAPI)
#
    CORE_SUBDIRS_EXTRA = lapack blas
         FC = mpxlf_r -qnohpf
# -F/u/d3g681/xlhpf.cfg:rjhxlf
         CC = mpcc_r
    ARFLAGS = urs
     RANLIB = echo
  MAKEFLAGS = -j 1 --no-print-directory
    INSTALL = @echo $@ is built
        CPP = /usr/lib/cpp -P
     MPILIB = 
LARGE_FILES = YES

ifeq ($(NWCHEM_TARGET_CPU),604)
  LDOPTIONS = -lxlf90_r -lm_r -qEXTNAME -qnosave -g -bmaxdata:0x20000000 -bloadmap:nwchem.lapi_map
   LINK.f   = mpxlf_r   $(LDFLAGS)
else
  LDOPTIONS = -lc_r -lxlf90_r -lm_r -qEXTNAME -qnosave -g -bmaxdata:0x20000000 -bloadmap:nwchem.lapi_map
   LINK.f   = mpcc_r   $(LDFLAGS)
endif
   FOPTIONS = -qEXTNAME -qnosave
# -qinitauto=7F # note that grad_force breaks with this option
   COPTIONS = 
  FOPTIMIZE = -O3 -qstrict -qfloat=rsqrt:fltint -NQ40000 -NT80000
  COPTIMIZE = -O
ifeq ($(NWCHEM_TARGET_CPU),P2SC)
# These from George from Kent Winchell
  FOPTIMIZE += -qcache=type=d:level=1:size=128:line=256:assoc=4:cost=14 \
        -qcache=type=i:level=1:size=32:line=128
  COPTIMIZE += -qcache=type=d:level=1:size=128:line=256:assoc=4:cost=14 \
        -qcache=type=i:level=1:size=32:line=128
endif
ifeq ($(NWCHEM_TARGET_CPU),604)
	FC += -qarch=604 -qtune=604 -qthreaded
	CC += -qarch=ppc -qtune=604
endif


    DEFINES = -DLAPI -DSP1 -DAIX -DEXTNAME -DPARALLEL_DIAG
#
# Prefix LIBPATH with -L/usr/lib for AIX 3.2.x
#
#  LIBPATH += -L/sphome/harrison/peigs2.0

  CORE_LIBS = -lchemio -lglobal -lma -lutil -lpeigs -llapack -lblas


   USE_ESSL = YES
#   USE_BLAS = YES
ifdef USE_ESSL
   DEFINES += -DESSL
# renames not needed for 4.1.  Still are for 3.2.
ifeq ($(NWCHEM_TARGET_CPU),P2SC)
 CORE_LIBS += -lpesslp2_t -lblacsp2_t -lesslp2_r
else
 CORE_LIBS += -lpessl -lblacs -lessl
endif

#	      -brename:.daxpy_,.daxpy \
#	      -brename:.dgesv_,.dgesv \
#	      -brename:.dcopy_,.dcopy \
#	      -brename:.ddot_,.ddot \
#	      -brename:.dgemm_,.dgemm \
#	      -brename:.dgemv_,.dgemv \
#	      -brename:.dgetrf_,.dgetrf \
#	      -brename:.dgetrs_,.dgetrs \
#	      -brename:.dscal_,.dscal \
#	      -brename:.dspsvx_,.dspsvx \
#	      -brename:.dpotrf_,.dpotrf \
#	      -brename:.dpotri_,.dpotri \
#	      -brename:.idamax_,.idamax 
ifdef USE_BLAS
 CORE_LIBS += -lblas -brename:.xerbla_,.xerbla -brename:.lsame_,.lsame
endif

else
    CORE_SUBDIRS_EXTRA += blas
             CORE_LIBS += -lblas
endif


# IMPORTANT:  These renames are necessary if you try to link against
# a copy of PeIGS built for MPI instead of TCGMSG. (Not recommended, 
# see INSTALL)
# mpipriv is a common block used in MPICH's implementation of MPI.  It
# is critical that this common block is renamed correctly because
# the linker will not detect any problems (there will be separate
# common blocks labeled mpipriv_ and mpipriv) but the program will not
# operate correctly.
#ifdef USE_MPI
#   CORE_LIBS += -brename:.mpi_recv_,.mpi_recv \
#		-brename:.mpi_initialized_,.mpi_initialized \
#		-brename:.mpi_init_,.mpi_init \
#		-brename:.mpi_comm_rank_,.mpi_comm_rank \
#		-brename:.mpi_comm_size_,.mpi_comm_size \
#		-brename:.mpi_finalize_,.mpi_finalize \
#		-brename:.mpi_send_,.mpi_send \
#		-brename:mpipriv_,mpipriv
#endif

 EXPLICITF = TRUE
  FCONVERT = $(CPP) $(CPPFLAGS) $< > $*.f
#
endif

ifeq ($(TARGET),DECOSF)
#
# DEC AXP OSF1
#
# JN 96/10/02:
# Replaced -DLongInteger with -DEXT_INT for consistency with GA, DRA, PEIGS ...

    CORE_SUBDIRS_EXTRA = blas lapack
                  NICE = nice
                SHELL := $(NICE) /bin/sh
                    FC = f77
                    AR = ar
                RANLIB = echo
	 	   CPP = /usr/bin/cpp -P -C	

               INSTALL = @echo nwchem is built
             MAKEFLAGS = -j 1 --no-print-directory

# -fpe2 and call to util/dec_fpe.f from nwchem.F necessary to avoid
# braindead alpha undflows inside texas (c6h6 6-31g)

#              FOPTIONS = -i8 -assume noaccuracy_sensitive -align dcommons -math_library fast -fpe2 -check nounderflow
# assume noaccuracy_sensitive was breaking the code in recent versions (EA)
              FOPTIONS = -i8 -align dcommons -math_library fast -fpe2 -check nounderflow
              COPTIONS = 
             FOPTIMIZE = -O 
             COPTIMIZE = -O

               DEFINES = -DDECOSF -DEXT_INT
             CORE_LIBS = -lutil -lchemio -lglobal -lma -lpeigs -llapack -lblas
            EXTRA_LIBS = -laio -lpthreads 
endif

ifeq ($(TARGET),LINUX)
#
# Most Linux distributions are using EGCS
#
  EGCS = YES
#
# Linux running on an x86 using g77
# f2c has not been tested in years and is not supported
#
       NICE = nice -2
      SHELL := $(NICE) /bin/sh
    CORE_SUBDIRS_EXTRA = blas lapack
         CC = gcc
     RANLIB = ranlib
  MAKEFLAGS = -j 2 --no-print-directory
    INSTALL = @echo $@ is built

ifeq ($(BUILDING_PYTHON),python)
   EXTRA_LIBS += -ltk -ltcl -L/usr/X11/lib -lX11 -ldl
   INCPATH += -I/usr/include/python1.5
endif

# defaults are for X86 platforms
         FC  = g77
  FOPTIONS   = -fno-second-underscore 
  FOPTIMIZE  = -g -O2 
  COPTIONS   = -Wall -m486 -malign-double
  COPTIMIZE  = -g -O2
ifdef EGCS
  FOPTIONS  += -fno-globals -Wunused -fno-silent -m486 -malign-double
  FOPTIMIZE += -Wuninitialized -ffast-math -funroll-loops -fstrength-reduce 
  FOPTIMIZE += -fno-move-all-movables -fno-reduce-all-givs -fno-rerun-loop-opt 
  FOPTIMIZE += -fforce-mem -fforce-addr -ffloat-store
endif

ifeq ($(NWCHEM_TARGET_CPU),ALPHA)
  FOPTIONS   = -fno-second-underscore
  FOPTIMIZE  = -g -O2 
  COPTIONS   = -Wall
  COPTIMIZE  = -g -O2
endif
ifeq ($(NWCHEM_TARGET_CPU),POWERPC)
  FOPTIONS   = -fno-second-underscore -fno-globals
  FOPTIMIZE  = -g -O2 
  COPTIONS   = -Wall
  COPTIMIZE  = -g -O2
endif

    DEFINES = -DLINUX

  LDOPTIONS = -g
     LINK.f = g77 $(LDFLAGS)
 EXTRA_LIBS += -lm
ifndef EGCS
 EXTRA_LIBS += -lf2c -lm
endif

  CORE_LIBS = -lutil -lchemio -lglobal -lma -lpeigs -llapack -lblas

        CPP = gcc -E -nostdinc -undef -P
   FCONVERT = (/bin/cp $< /tmp/$$$$.c; \
			$(CPP) $(CPPFLAGS) /tmp/$$$$.c | sed '/^$$/d' > $*.f; \
			/bin/rm -f /tmp/$$$$.c) || exit 1


endif

ifeq ($(TARGET),FUJITSU_VPP)
#
# FUJITSU VX/VPP
#
# HAF Sept. 97

         FC = frt
      CPP = /lib/cpp -P -C
     RANLIB = echo
  MAKEFLAGS = 
    INSTALL = @echo $@ is built
                        
    DEFINES = -DFUJITSU_VPP
    USE_MPI = TRUE

#change DEFINES so that frt understands them and simply add them to FOPTIONS
 FDEFINES_1 = $(strip  $(DEFINES))
   FDEFINES = -Wp,$(subst $(space),$(comma),$(FDEFINES_1))   
   FOPTIONS = -w -Sw $(FDEFINES)
   COPTIONS = 
     FDEBUG = -Ob -g
  FOPTIMIZE = -Kfast -Wv,-s8
  COPTIMIZE = -K4

# removed global, ma, tcgmsg-mpi, as they are part of the native GA
NW_CORE_SUBDIRS = include basis geom inp input chemio ma \
	pstat rtdb task symmetry util $(CORE_SUBDIRS_EXTRA)

       CORE_LIBS = -lutil -lchemio -lpeigs \
                   -L/home/fruechtl/lib -lglobal -lma -ltcgmsg-mpi \
                   -llapackvp -lblasvp
      EXTRA_LIBS = -L /opt/tools/lib/ -lmp -lgen  -lpx -lelf -Wl,-J,-P,-t
endif

ifeq ($(TARGET),PGLINUX)
#
# Linux running on an x86 using g77
# to use f2c/gcc, define environment variable USE_F2C
#
       NICE = nice
      SHELL := $(NICE) /bin/sh
    CORE_SUBDIRS_EXTRA = blas lapack
         CC = gcc
     RANLIB = ranlib
  MAKEFLAGS = -j 1 --no-print-directory
    INSTALL = @echo $@ is built

  FOPTIONS  = -Mdalign -Minform,warn -Mnolist 
#         FC = sleep 2;pgf77
          FC = pgf77
#         FC = sleep 2;pghpf -Mf90

   COPTIONS =  -Wall -m486 -malign-double
ifeq ($(NWCHEM_TARGET_CPU),604)
   COPTIONS =  -Wall
endif
  FOPTIMIZE = -O2
  COPTIMIZE = -g -02

    DEFINES = -DLINUX -DPGLINUX

  LDOPTIONS = -g
     LINK.f = pgf77 $(LDFLAGS)
  CORE_LIBS = -lutil -lchemio -lglobal -lma -lpeigs -llapack -lblas
 EXTRA_LIBS = 

        CPP = gcc -E -nostdinc -undef -P
   FCONVERT = (/bin/cp $< /tmp/$$$$.c; \
			$(CPP) $(CPPFLAGS) /tmp/$$$$.c | sed '/^$$/d' > $*.f; \
			/bin/rm -f /tmp/$$$$.c) || exit 1
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
CORE_LIBS += $(PYTHONHOME)/lib/python$(PYTHONVERSION)/config/libpython$(PYTHONVERSION).a
endif

###################################################################
#  All machine dependent sections should be above here, otherwise #
#  some of the definitions below will be 'lost'                   #
###################################################################
# MPI version requires tcgmsg-mpi library
ifdef USE_MPI
   ifdef MPI_LIB
       LIBPATH += -L$(MPI_LIB)
   endif
   CORE_LIBS += -ltcgmsg-mpi $(MPILIB)
else
    CORE_LIBS += -ltcgmsg
endif

#the new GA uses ARMCI library
ifndef OLD_GA
      CORE_LIBS += -larmci
endif


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
.SUFFIXES:	.o .s .F .f .c

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
	$(FC) -c $(FFLAGS) $(CPPFLAGS) $<
endif

(%.o):	%.f
	$(FC) -c $(FFLAGS) $<

(%.o):	%.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) -o $% $<

(%.o):  %.o

# Preceding line has a tab to make an empty rule

# a .F.f rule is needed for any system target where the default .F.f rule does not work
# AND the EXPLICITF is not already true.  Right now this is only LINUX with g77
ifeq ($(TARGET),LINUX)
.F.f:
	$(FC) -c $(FFLAGS) -E $(CPPFLAGS) $< -o $*.f
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


