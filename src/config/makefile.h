#
# $Id: makefile.h,v 1.354 2001-04-13 22:25:25 edo Exp $
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
# Note that the preprocessor flags for CRAY-T3D and CRAY-T3E are CRAY_T3D and CRAY_T3E respectively
#

ifndef NWCHEM_TARGET
error2:
	@echo You must define NWCHEM_TARGET in your environment to be the name
	@echo of the machine you wish to build for ... for example
	@echo     setenv NWCHEM_TARGET SOLARIS
	@echo Known targets are SOLARIS, SGI_N32, ...
	@echo See the INSTALL instructions for a complete list
	@exit 1
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
    LIBPATH = -L$(SRCDIR)/tools/lib/$(TARGET)

#
# Define INCPATH to be directories to get includes for
# libraries that you are not building now.  These directories
# will be searched AFTER anything you are building now.
#
    INCPATH = 
    INCPATH = -I$(SRCDIR)/tools/include

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

       CORE_LIBS =  -lutil -lpario -lglobal -lma -lpeigs 
#
# Machine specific stuff
#

ifeq ($(TARGET),SOLARIS)
      SHELL := $(NICE) /bin/sh
	CPP = /usr/ccs/lib/cpp
     RANLIB = echo
  MAKEFLAGS = -j 4 --no-print-directory
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
ifdef BLASOPT
CORE_SUBDIRS_EXTRA = blas lapack 
endif

# Note that WS6 does not optimize robustly and if using this you must 
#   - put "-nodpend -xvector=no" on FOPTIONS after -fast
#   - remove "-fsimple=2 -depend -xvector=yes" from FOPTIMIZE.
#   - remove -lmvec from CORELIBS
# Also in WS6 -fast implies -native so NWCHEM_TARGET_CPU does not
# really add anything.
#
# These options are set for WS5

   COPTIONS = 
  COPTIMIZE = -g -O
   FOPTIONS = -stackvar -dalign 
  FOPTIMIZE = -fast -O5 -fsimple=2 -depend -xvector=yes
     FDEBUG = -g -O1 -nodepend

ifeq ($(NWCHEM_TARGET_CPU), ULTRA)
  FOPTIMIZE += -xarch=v8plusa 
endif

    DEFINES = -DSOLARIS  -DNOAIO  -DPARALLEL_DIAG
  LDOPTIONS = -xildoff
  LINK.f = $(FC) $(LDFLAGS) $(FOPTIONS)
ifndef BLASOPT
  CORE_LIBS +=  -xlic_lib=sunperf -lmvec
else
  CORE_LIBS +=  -llapack $(BLASOPT) -lblas  -lmvec
endif

      EXTRA_LIBS = -ldl 
# this creates a static executable
#EXTRA_LIBS = -Bdynamic -ldl -lXext -lnsl  -Bstatic  

ifdef LARGE_FILES
  LDOPTIONS  += $(shell getconf LFS_LDFLAGS)
  EXTRA_LIBS += $(shell getconf LFS_LIBS)
  DEFINES    += -DLARGE_FILES
endif

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
	CPP = /usr/ccs/lib/cpp
         CC = cc
         FC = f77

   COPTIONS = -xarch=v9 -dalign
  COPTIMIZE = -O
     RANLIB = echo
  MAKEFLAGS = -j 28 --no-print-directory
    INSTALL = echo $@ is built

# These options are set for WS6.1

   FOPTIONS = -stackvar -fast -nodepend -xvector=no -xarch=v9a -xtypemap=real:64,double:64,integer:64
  FOPTIMIZE = -g -O5
     FDEBUG = -g -O1

    DEFINES = -DSOLARIS  -DNOAIO -DSOLARIS64 -DEXT_INT -DPARALLEL_DIAG

  LDOPTIONS = -xs -xildoff
  LINK.f = $(FC) $(LDFLAGS) $(FOPTIONS)
       CORE_LIBS += -llapack $(BLASOPT)   -lblas 
      EXTRA_LIBS =  -ldl -lfsu

ifdef LARGE_FILES
  LDOPTIONS  += $(shell getconf LFS_LDFLAGS)
  EXTRA_LIBS += $(shell getconf LFS_LIBS)
  DEFINES    += -DLARGE_FILES
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

            CORE_LIBS += -llapack $(BLASOPT) -lblas 

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

              DEFINES = -DCRAY_T3E -DCRAY_T3D -D__F90__ -DPARALLEL_DIAG

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
 FOPTIMIZE_10K = -O3 -OPT:const_copy_limit=20000:Olimit=4000 -TENV:X=1 -WK,-so=1,-o=1,-r=3,-dr=AKC
FVECTORIZE_10K = -O3 -TENV:X=1 -WK,-dr=AKC

#optimization flags for R12000 (IP30)
 FOPTIMIZE_12K = -O3 -OPT:const_copy_limit=20000:Olimit=4000 -TENV:X=1 -WK,-so=1,-o=1,-r=3,-dr=AKC
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

  DEFINES = -DSGI -DSGITFP -DEXT_INT -DPARALLEL_DIAG
  CORE_LIBS += -llapack $(BLASOPT) -lblas
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

       CORE_LIBS = -llapack $(BLASOPT) -lblas
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
    DEFINES = -DSGI  -DSGI_N32 -DPARALLEL_DIAG

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
endif

ifeq ($(TARGET),HPUX)
#
# HPUX 11.0
#

  CORE_SUBDIRS_EXTRA = blas lapack
  MAKEFLAGS = -j 1 --no-print-directory
  CPP = /lib/cpp -P
  CC = cc
#  FC = f77
  FC = f90
  LDOPTIONS = -Wl,+vallcompatwarnings -g -L/usr/lib +U77
#  LINK.f = fort77   $(LDFLAGS) 
  LINK.f = f90   $(LDFLAGS) 
  CORE_LIBS +=   -llapack $(BLASOPT) -lblas   -lm
  FDEBUG = -g
  FOPTIONS =  +ppu
  COPTIONS = -Aa -D_HPUX_SOURCE +e
  FOPTIMIZE = +O2
  FVECTORIZE = +Oall +Onofltacc
  COPTIMIZE = -O
  RANLIB = echo

 DEFINES = -DHPUX -DEXTNAME -DPARALLEL_DIAG
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
  MAKEFLAGS = -j 1 --no-print-directory
  CPP = /lib/cpp -P
  CC = cc
  FC = f90
  LDOPTIONS = -Wl,+vallcompatwarnings
  CORE_LIBS +=  -llapack $(BLASOPT) -lblas -lm
  CDEBUG =
  FDEBUG = -g
  FOPTIONS =  +ppu +i8 +DA2.0W +U77  
  COPTIONS = -Aa -D_HPUX_SOURCE +e +DA2.0W 
  FOPTIMIZE = +O2
  FVECTORIZE = +Oall +Onofltacc
  COPTIMIZE = -O
  RANLIB = echo

 DEFINES = -DHPUX -DEXTNAME -DPARALLEL_DIAG -DHPUX64 -DEXT_INT
# DEFINES = -DHPUX -DEXTNAME -DEXT_INT

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

   FOPTIONS = -qEXTNAME -qnosave -qalign=4k 
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
              -brename:.ztrmv_,.ztrmv 
#              -brename:.zherk_,.zherk \
#              -brename:.zhemm_,.zhemm \
#              -brename:.ztrsm_,.ztrsm \
#              -brename:.dtrsv_,.dtrsv \
#              -brename:.ztrsv_,.ztrsv 
#              -brename:.dsymm_,.dsymm \
#

##comment out from dtrmm_ inclusive
#ifdef USE_ESSL
#       CORE_LIBS += -lessl
#endif

  EXPLICITF = TRUE
  FCONVERT = $(CPP) $(CPPFLAGS) $< > $*.f
#
endif

ifeq ($(TARGET),IBM64)
# 
# IBM AIX 64-bit
# tested on ecs1 May 10 2000 AIX 4.3.3.10
# does not run on AIX  4.3.2.1 (skunkworks)
#

    CORE_SUBDIRS_EXTRA = lapack blas
         FC = xlf
         CC = xlc
    ARFLAGS = -X 64 urs 
     RANLIB = echo
  MAKEFLAGS = -j 3 --no-print-directory
    INSTALL = @echo $@ is built
        CPP = /usr/lib/cpp -P

   FOPTIONS = -qEXTNAME -qnosave -qalign=4k -q64 -qintsize=8 
   COPTIONS = -q64
  FOPTIMIZE = -O3 -qstrict -qfloat=rsqrt:fltint -NQ40000 -NT80000  -qarch=auto -qtune=auto
  FVECTORIZE = -O5 -qhot -qfloat=rsqrt:fltint:hssngl -NQ40000 -NT80000  -qarch=auto -qtune=auto 
   COPTIMIZE = -O -qmaxmem=8192

    DEFINES = -DIBM -DAIX -DEXTNAME -DEXT_INT -DIBM64
       LIBPATH += -L/usr/lib -L/lib 

       CORE_LIBS += -llapack $(BLASOPT) -lblas

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
  LDOPTIONS = -lxlf90_r -lm_r -qEXTNAME -qnosave -qalign=4k -g -bmaxdata:0x20000000 -bloadmap:nwchem.lapi_map
   LINK.f   = mpxlf_r   $(LDFLAGS)
else
  LDOPTIONS = -lc_r -lxlf90_r -lm_r -qEXTNAME -qnosave -qalign=4k -g -bmaxdata:0x40000000 -bloadmap:nwchem.lapi_map
   LINK.f   = mpcc_r   $(LDFLAGS)
endif
   FOPTIONS = -qEXTNAME -qnosave -qalign=4k 
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
ifeq ($(NWCHEM_TARGET_CPU),PWR3)
	FC += -qarch=pwr3 -qtune=pwr3 -qcache=auto -qthreaded
	CC += -qarch=pwr3 -qtune=pwr3 -qcache=auto
endif


    DEFINES = -DLAPI -DSP1 -DAIX -DEXTNAME -DPARALLEL_DIAG
#
# Prefix LIBPATH with -L/usr/lib for AIX 3.2.x
#
#  LIBPATH += -L/sphome/harrison/peigs2.0


USE_ESSL = YES
ifdef USE_ESSL
   DEFINES += -DESSL
ifeq ($(NWCHEM_TARGET_CPU),P2SC)
 CORE_LIBS += -lpesslp2_t -lblacsp2_t -lesslp2_r
else
 CORE_LIBS += -lpessl -lblacs -lessl
endif
# renames not needed for 4.1.  Still are for 3.2.
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
endif

# Need ESSL before our own BLAS library but still need our
# own stuff for misc. missing routines

CORE_LIBS +=  -llapack -lblas

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

ifeq ($(TARGET),LAPI64)
#
    CORE_SUBDIRS_EXTRA = lapack blas
         FC = mpxlf_r 
         CC = mpcc_r
    ARFLAGS = urs
     RANLIB = echo
  MAKEFLAGS = -j 1 --no-print-directory
    INSTALL = @echo $@ is built
        CPP = /usr/lib/cpp -P
     MPILIB = 
LARGE_FILES = YES

  LDOPTIONS = -lc_r -lxlf90_r -lm_r -qEXTNAME -qnosave -q64 -g -bmaxdata:0x40000000 -bloadmap:nwchem.lapi64_map
   LINK.f   = mpcc_r   $(LDFLAGS)

   FOPTIONS = -qEXTNAME -qnosave -q64 #-qalign=4k 
       AR   = ar -X 64
   COPTIONS = -q64
  FOPTIMIZE = -O3 -qstrict -qfloat=rsqrt:fltint -NQ40000 -NT80000
  FOPTIMIZE += -qarch=auto -qtune=auto
  COPTIMIZE = -O

    DEFINES = -DEXT_INT -DLAPI64  -DIBM64 -DEXTNAME -DLAPI -DSP1 -DAIX -DEXTNAME -DPARALLEL_DIAG
 CORE_LIBS += -lessl_r # need 64-bit essl

CORE_LIBS +=  -llapack -lblas

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
              FOPTIONS = -i8 -align dcommons -math_library fast -fpe2 -check nounderflow -check nopower -check nooverflow  -warn argument_checking -warn unused -automatic

             COPTIONS = 
             LDOPTIONS = -O
             LINK.f = f77 $(LDFLAGS)
             FOPTIMIZE =  -O4  -tune host -arch host  -math_library fast
             FVECTORIZE = -fast -O5 -tune host -arch host
             COPTIMIZE = -O

               DEFINES = -DDECOSF -DEXT_INT -DPARALLEL_DIAG
             CORE_LIBS +=  -llapack $(BLASOPT) -lblas 
            EXTRA_LIBS = -laio 
ifeq ($(BUILDING_PYTHON),python)
      EXTRA_LIBS += -lX11

ifdef LARGE_FILES
  LDOPTIONS  += $(shell getconf LFS_LDFLAGS)
  EXTRA_LIBS += $(shell getconf LFS_LIBS)
  DEFINES    += -DLARGE_FILES
endif


endif
endif

ifeq ($(TARGET),$(findstring $(TARGET),LINUX CYGNUS))
#
#
# Linux or Cygwin under Windows running on an x86 using g77
# f2c has not been tested in years and is not supported
#
       NICE = nice -2
      SHELL := $(NICE) /bin/sh
    CORE_SUBDIRS_EXTRA = blas lapack
         CC = gcc
     RANLIB = ranlib
  MAKEFLAGS = -j 1 --no-print-directory
    INSTALL = @echo $@ is built

ifeq ($(BUILDING_PYTHON),python)
   EXTRA_LIBS += -ltk -ltcl -L/usr/X11R6/lib -lX11 -ldl
# needed if python was built with pthread support
#   EXTRA_LIBS += -lpthread
   INCPATH += -I/usr/include/python1.5
endif

  DEFINES = -DLINUX -DPARALLEL_DIAG
ifeq ($(TARGET),CYGNUS)
    DEFINES += -DCYGNUS
endif

ifeq ($(FC),pgf77)
  DEFINES   += -DPGLINUX
# added -Kieee to get dlamc1 to work on pgf77 3.1-3 EA Jun 8th 2000
  FOPTIONS   = -Mdalign -Minform,warn -Mnolist -Minfo=loop -Munixlogical -Kieee
# for pentium
  FOPTIONS  += -tp p5  
# for Pentium Pro or Pentium II
# FOPTIONS  += -tp p6
  FOPTIMIZE  = -O2 -Mvect
  COPTIONS   = -Wall -m486 -malign-double
  COPTIMIZE  = -g -O2
  MAKEFLAGS += FC=pgf77
else
# defaults are for X86 platforms
         FC  = g77
  FOPTIONS   = -fno-second-underscore 
  FOPTIMIZE  = -g -O2 -Wuninitialized 
  COPTIONS   = -Wall -m486 -malign-double 
  COPTIMIZE  = -g -O2
# Most Linux distributions are using EGCS
#
  EGCS = YES
ifdef EGCS
  FOPTIONS  += -Wno-globals
  FOPTIONS  += -fno-globals -Wunused -fno-silent -m486 -malign-double
  FOPTIMIZE += -Wuninitialized -ffast-math -funroll-loops -fstrength-reduce 
  FOPTIMIZE += -fno-move-all-movables -fno-reduce-all-givs -fno-rerun-loop-opt 
  FOPTIMIZE += -fforce-mem -fforce-addr -ffloat-store
endif

ifeq ($(NWCHEM_TARGET_CPU),POWERPC)
  FOPTIONS   = -fno-second-underscore -fno-globals -Wno-globals
  FOPTIMIZE  = -g -O2 
  COPTIONS   = -Wall
  COPTIMIZE  = -g -O2
endif
endif


ifeq ($(FC),pgf77)
  LDOPTIONS = -g
     LINK.f = pgf77 $(LDFLAGS)
 EXTRA_LIBS += -lm
else
  LDOPTIONS = -g -Xlinker -export-dynamic
     LINK.f = g77 $(LDFLAGS)
 EXTRA_LIBS += -lm
ifndef EGCS
 EXTRA_LIBS += -lf2c -lm
endif
endif

CORE_LIBS += -llapack $(BLASOPT) -lblas

        CPP = gcc -E -nostdinc -undef -P
   FCONVERT = (/bin/cp $< /tmp/$$$$.c; \
			$(CPP) $(CPPFLAGS) /tmp/$$$$.c | sed '/^$$/d' > $*.f; \
			/bin/rm -f /tmp/$$$$.c) || exit 1

# end of Linux, Cygnus
endif

ifeq ($(NWCHEM_TARGET),LINUX64)
# using COMPAQ/DEC compilers (EA 3/13/2000)
    CORE_SUBDIRS_EXTRA = blas lapack
     RANLIB = echo
  FC         = fort
  CC         = ccc      
# this creates a static executable
#  LINK.f = fort $(LDFLAGS)   -Wl,-Bstatic
  LINK.f = fort $(LDFLAGS)  
  DEFINES   +=   -DEXT_INT -DLINUX -DLINUX64 -DPARALLEL_DIAG
  FOPTIONS   = -i8 -assume no2underscore -align dcommons -fpe3 -check nooverflow -assume accuracy_sensitive -check nopower -check nounderflow  -noautomatic
  EXTRA_LIBS = 
  FOPTIMIZE =  -O4  -tune host -arch host  -math_library fast
  FVECTORIZE = -fast -O5 -tune host -arch host
  CORE_LIBS += -llapack $(BLASOPT) -lblas
ifeq ($(BUILDING_PYTHON),python)
#   EXTRA_LIBS += -ltk -ltcl -L/usr/X11R6/lib -lX11 -ldl
   EXTRA_LIBS +=  -ldl
endif
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
  FOPTIMIZE = -Kfast -Wv,-s8
  COPTIMIZE = -K4

# removed global, ma, tcgmsg-mpi, as they are part of the native GA
 NW_CORE_SUBDIRS = include basis geom inp input  \
       pstat rtdb task symmetry util peigs $(CORE_SUBDIRS_EXTRA)

        CORE_LIBS = -lutil -lpeigs \
                    -L$(GA_LIBDIR) -lglobal -lpario -lma \
                    -ltcgmsg-mpi -L/usr/lang/mpi2/lib32 -lmpi -lmp
       EXTRA_LIBS = -llapackvp -lblasvp -lsocket -Wl,-J,-P,-t,-dy
endif

#-do not use# ifeq ($(TARGET),PGLINUX)
#-do not use# #
#-do not use# # Linux running on an x86 using g77
#-do not use# # to use f2c/gcc, define environment variable USE_F2C
#-do not use# #
#-do not use#        NICE = nice
#-do not use#       SHELL := $(NICE) /bin/sh
#-do not use#     CORE_SUBDIRS_EXTRA = blas lapack
#-do not use#          CC = gcc
#-do not use#      RANLIB = ranlib
#-do not use#   MAKEFLAGS = -j 1 --no-print-directory
#-do not use#     INSTALL = @echo $@ is built
#-do not use#
#-do not use#   FOPTIONS  = -Mdalign -Minform,warn -Mnolist 
#-do not use# #         FC = sleep 2;pgf77
#-do not use#           FC = pgf77
#-do not use# #         FC = sleep 2;pghpf -Mf90
#-do not use#
#-do not use#    COPTIONS =  -Wall -m486 -malign-double
#-do not use# ifeq ($(NWCHEM_TARGET_CPU),604)
#-do not use#    COPTIONS =  -Wall
#-do not use# endif
#-do not use#   FOPTIMIZE = -O2
#-do not use#   COPTIMIZE = -g -02
#-do not use#
#-do not use#     DEFINES = -DLINUX -DPGLINUX
#-do not use#
#-do not use#   LDOPTIONS = -g
#-do not use#      LINK.f = pgf77 $(LDFLAGS)
#-do not use#   CORE_LIBS = -lutil -lpario -lglobal -lma -lpeigs -llapack -lblas
#-do not use#  EXTRA_LIBS = 
#-do not use#
#-do not use#         CPP = gcc -E -nostdinc -undef -P
#-do not use#    FCONVERT = (/bin/cp $< /tmp/$$$$.c; \
#-do not use# 			$(CPP) $(CPPFLAGS) /tmp/$$$$.c | sed '/^$$/d' > $*.f; \
#-do not use# 			/bin/rm -f /tmp/$$$$.c) || exit 1
#-do not use# endif


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
ifndef LIBMPI 
  LIBMPI = -lmpi 
endif 
ifdef MPI_LIB 
      CORE_LIBS += -L$(MPI_LIB) 
endif 
  CORE_LIBS += -ltcgmsg-mpi $(LIBMPI) 
else 
  CORE_LIBS += -ltcgmsg 
endif 
# lower level libs used by communication libraries 
 
#the new GA uses ARMCI library
ifndef OLD_GA
      CORE_LIBS += -larmci
endif


ifdef COMM_LIBS 
 CORE_LIBS += $(COMM_LIBS) 
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


