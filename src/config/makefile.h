#
# $Id: makefile.h,v 1.486 2004-10-15 00:14:31 edo Exp $
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
#RELEASE := 4.6

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
	pstat rtdb task symmetry util peigs perfm cons $(CORE_SUBDIRS_EXTRA)

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

       CORE_LIBS =  -lnwcutil -lpario -lglobal -lma -lpeigs -lperfm -lcons
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

    DEFINES = -DSOLARIS  -DNOAIO  -DPARALLEL_DIAG
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
   DEFINES = -DSOLARIS  -DNOAIO -DSOLARIS64 -DPARALLEL_DIAG
ifdef USE_INTEGER4
else
   DEFINES  +=  -DEXT_INT
endif

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
     ifdef USE_INTEGER4
       FOPTIONS += -CcdLL8
     else
       FOPTIONS += -CcdLL8 -CcdII8
     endif
     FOPTIMIZE = -Kfast_GP=2  -KV9FMADD
     FDEBUG=
   else
#  SUN/Solaris f77 options 
     FOPTIONS = -stackvar -fast -nodepend -xvector=no -xarch=v9a
     ifdef USE_INTEGER4
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
              DEFINES = -DCRAY_T3D -DPARALLEL_DIAG -DUSE_FCD

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

              DEFINES = -DCRAY_T3E -DCRAY_T3D -D__F90__ -DPARALLEL_DIAG -DUSE_FCD

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

  DEFINES = -DSGI -DSGITFP -DEXT_INT -DPARALLEL_DIAG
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
  _CPU = $(shell uname -m  )
  MAKEFLAGS = -j 1 --no-print-directory
  CPP = /lib/cpp -P
  CC = cc
  FC = f90
  LDOPTIONS = -Wl,+vallcompatwarnings 
  CORE_LIBS +=  -llapack $(BLASOPT) -lblas  -lm
  CDEBUG =
  FDEBUG = -g
  FOPTIONS =  +ppu  +U77  
  COPTIONS = -Aa -D_HPUX_SOURCE +e 
  ifeq ($(_CPU),ia64)
    FOPTIONS += +DD64 +DSitanium2 +Ofltacc=relaxed +Olibcalls +Onolimit +FPD
    COPTIONS += +DD64
    FOPTIMIZE = +O3 +Oloopblock +Oinline_budget=200 
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

 DEFINES = -DHPUX -DEXTNAME -DPARALLEL_DIAG -DHPUX64 
 ifdef USE_INTEGER4
 else
 DEFINES +=  -DEXT_INT
 FOPTIONS +=  +i8 
 endif

endif




ifeq ($(TARGET),IBM)
#
# IBM AIX
#

    CORE_SUBDIRS_EXTRA = lapack blas
         FC = xlf
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

    DEFINES = -DIBM -DAIX -DEXTNAME -DPARALLEL_DIAG
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
         AR = ar -X 64
     RANLIB = echo
  MAKEFLAGS = -j 11 --no-print-directory
    INSTALL = @echo $@ is built
        CPP = /usr/lib/cpp -P

   FOPTIONS = -qEXTNAME -qnosave -qalign=4k -q64 -qxlf77=leadzero
   COPTIONS = -q64
  FOPTIMIZE = -O3 -qstrict -NQ40000 -NT80000  -qarch=auto -qtune=auto
  ifdef RSQRT
    FOPTIMIZE  += -qfloat=rsqrt:fltint
    FVECTORIZE  += -qfloat=rsqrt:fltint
  endif
   COPTIMIZE = -O -qmaxmem=8192

    DEFINES = -DIBM -DAIX -DEXTNAME -DPARALLEL_DIAG
    DEFINES += -DCHKUNDFLW
       LIBPATH += -L/usr/lib -L/lib 
ifdef USE_INTEGER4
   FOPTIONS += -qintsize=4
   ifdef USE_ESSL
      DEFINES += -DESSL
      CORE_LIBS += -lessl 
   endif
else
   FOPTIONS += -qintsize=8 
   DEFINES  += -DEXT_INT 
endif
  ifdef  USE_GPROF
    FOPTIONS += -pg
    LDOPTIONS += -pg
  endif
  LDOPTIONS += -bloadmap:nwchem.ibm64map -bbigtoc
  LDOPTIONS += -bmaxstack:0x80000000 -bmaxdata:0x80000000 # needed because of bigtoc
   CORE_LIBS += -llapack $(BLASOPT) -lblas


  EXPLICITF = TRUE
  FCONVERT = $(CPP) $(CPPFLAGS) $< > $*.f
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


    DEFINES = -DLAPI -DSP1 -DAIX -DEXTNAME -DPARALLEL_DIAG
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

  LDOPTIONS = -lc_r -lxlf90_r -lm_r -qEXTNAME -qnosave -q64  -bloadmap:nwchem.lapi64_map $(LAPI64LIBS)
   LINK.f   = mpcc_r   $(LDFLAGS)

   FOPTIONS = -qEXTNAME -qnosave -q64 -qalign=4k -qxlf77=leadzero -qthreaded
       AR   = ar -X 64
   COPTIONS = -q64
  FOPTIMIZE = -O3 -qstrict -NQ40000 -NT80000
  FOPTIMIZE += -qarch=auto -qtune=auto -qcache=auto
  ifdef RSQRT
    FOPTIMIZE  += -qfloat=rsqrt:fltint
  endif
  COPTIMIZE = -O

    DEFINES = -DLAPI64 -DEXTNAME -DLAPI -DSP1 -DAIX -DPARALLEL_DIAG
    DEFINES += -DCHKUNDFLW
ifdef USE_INTEGER4
   FOPTIONS += -qintsize=4
 CORE_LIBS += -lessl_r -llapack -lblas # need 64-bit essl
else
   FOPTIONS += -qintsize=8
        DEFINES += -DEXT_INT
  CORE_LIBS +=  -llapack -lblas
endif
  LDOPTIONS += -bloadmap:nwchem.lapi64map -bbigtoc
  LDOPTIONS += -bmaxstack:0x80000000 -bmaxdata:0x80000000 # needed because of bigtoc


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
              FDEBUG = -g -O0
              FOPTIONS = -align dcommons -math_library fast -fpe2 -check nounderflow -check nopower -check nooverflow  -warn argument_checking -warn unused -automatic -math_library fast

             COPTIONS = 
             LDOPTIONS = -O
             LINK.f = f77 $(LDFLAGS)
             FOPTIMIZE =  -O4  -tune host -arch host  
             FVECTORIZE = -fast -O4 -tune host -arch host
             COPTIMIZE = -O

               DEFINES = -DDECOSF -DPARALLEL_DIAG
ifdef USE_INTEGER4
              FOPTIONS +=  -i4
             CORE_LIBS +=   -lcxml_ev6
else
               DEFINES +=  -DEXT_INT 
              FOPTIONS +=  -i8
             CORE_LIBS +=  -llapack $(BLASOPT) -lblas 
endif
            EXTRA_LIBS = -laio 
ifeq ($(BUILDING_PYTHON),python)
      EXTRA_LIBS += -lX11
endif
endif
ifeq ($(TARGET),MACX)
#
# MacOSX 
#
# 

ifdef USE_VECLIB
    CORE_SUBDIRS_EXTRA =  blas
else
    CORE_SUBDIRS_EXTRA =  blas lapack
endif
               _CPU = $(shell machine  )
                    FC = g77
               INSTALL = @echo nwchem is built
               RANLIB = ranlib
             MAKEFLAGS = -j 1 --no-print-directory
             DEFINES =-DMACX -DPARALLEL_DIAG

  ifeq ($(FC),xlf)
    XLFMAC=y
    FOPTIONS = -qextname -qfixed -qnosave  -qalign=4k
    FOPTIONS +=  -NQ40000 -NT80000 -NS2048 -qmaxmem=8192 -qxlf77=leadzero
    FOPTIMIZE= -O3 -qstrict  -qarch=auto -qtune=auto -qcache=auto
    ifdef RSQRT
      FOPTIMIZE  += -qfloat=rsqrt:fltint
    endif
    FOPTIMIZE+=  -qunroll=yes 
#    FDEBUG= -O2 
    DEFINES  +=-DXLFLINUX -DCHKUNDFLW
     FOPTIONS += $(INCLUDES) -WF,"$(DEFINES)" $(shell echo $(LIB_DEFINES) | sed -e "s/-D/-WF,-D/g"   | sed -e 's/\"/\\\"/g'  | sed -e "s/\'/\\\'/g")
  else
#g77, only decent one form Fink http://fink.sf.net
#gcc version 3.4 20031015 (experimental)
    _G77V33= $(shell g77 -v  2>&1|egrep spec|head -1|awk ' /3.3/  {print "Y"}')
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
    ifeq ($(CC),xlc)
      COPTIONS  +=  -qlanglvl=extended
    else
      COPTIONS   = -Wall -no-cpp-precomp
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
    EXTRA_LIBS += -lm -lcc_dynamic
  endif

endif


ifeq ($(TARGET),$(findstring $(TARGET),LINUX CYGNUS CYGWIN INTERIX))
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
        CPP = gcc -E -nostdinc -undef -P
   FCONVERT = (/bin/cp $< /tmp/$$$$.c; \
			$(CPP) $(CPPFLAGS) /tmp/$$$$.c | sed '/^$$/d' > $*.f; \
			/bin/rm -f /tmp/$$$$.c) || exit 1

         LINUXCPU = $(shell uname -m |\
                 awk ' /sparc/ { print "sparc" }; /i*86/ { print "x86" };  /ppc*/ { print "ppc"} ' )

ifeq ($(BUILDING_PYTHON),python)
#   EXTRA_LIBS += -ltk -ltcl -L/usr/X11R6/lib -lX11 
#   EXTRA_LIBS += -L/home/edo/tcltk/lib/LINUX -ltk8.3 -ltcl8.3 -L/usr/X11R6/lib -lX11 -ldl
# needed if python was built with pthread support
   EXTRA_LIBS += -lz  -lreadline -lncurses -lnwcutil  -lpthread -lutil -ldl
endif

  DEFINES = -DLINUX -DPARALLEL_DIAG

ifeq ($(LINUXCPU),x86) 
  ifeq ($(TARGET),CYGNUS)
    DEFINES += -DCYGNUS
  endif
  ifeq ($(TARGET),CYGWIN)
    DEFINES += -DCYGWIN -DCYGNUS
  endif
  
  _CPU = $(shell uname -m  )
  _G77V33= $(shell g77 -v  2>&1|egrep spec|head -1|awk ' /3.3/  {print "Y"}')

# FC  = g77

    ifeq ($(_CPU),i686)
     _GOTSSE2= $(shell cat /proc/cpuinfo | egrep sse2 | tail -1 | awk ' /sse2/  {print "Y"}')
     _PENTIUM_M= $(shell cat /proc/cpuinfo | egrep " M processor" | tail -1 | awk ' /M/  {print "Y"}')
      ifeq ($(_GOTSSE2),Y) 
        _CPU=i786
      endif
    endif

    ifeq ($(_CPU),i786)
      COPTIONS   =  -march=i686 
      ifdef USE_GCC31
        FDEBUG=-O1 -g
        COPTIMIZE +=-march=pentium4 -mcpu=pentium4 -msse2 -mfpmath=sse 
        COPTIMIZE +=-fprefetch-loop-arrays -minline-all-stringops -fexpensive-optimizations
        FOPTIMIZE +=-march=pentium4 -mcpu=pentium4 -msse2 -mfpmath=sse 
        FOPTIMIZE +=-fprefetch-loop-arrays -minline-all-stringops -fexpensive-optimizations
      else
        FOPTIMIZE  +=  -march=i686
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
  FOPTIONS   += -fno-second-underscore   
  FOPTIONS   += -fno-f90  -ffixed-line-length-72 -ffixed-form
  FOPTIMIZE  +=  -O2  -malign-double -finline-functions 
  COPTIONS   += -Wall  -malign-double 
  COPTIMIZE  += -g -O2
    FOPTIONS  +=  -malign-double -fno-globals -Wno-globals  -Wunused  -fno-silent
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
 _FC=noifc
 ifeq ($(FC),ifc)
     _FC=ifc
 endif
 ifeq ($(FC),ifort)
     _FC=ifc
 endif
  ifeq ($(_FC),ifc)
  FOPTIONS   =  -align    -mp1 -w -g -vec_report3
  ifdef  USE_GPROF
    FOPTIONS += -qp
  endif
    _IFCV80= $(shell ifc -v  2>&1|egrep 8|awk ' /8\.0/  {print "Y"}')
    _IFCV8= $(shell ifc -v  2>&1|egrep 8|awk ' /8\./  {print "Y"}')
    ifeq ($(_IFCV8),Y)
      DEFINES+= -DIFCV8
      ifeq ($(FC),ifc)
          FOPTIONS += -quiet
      endif
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
      ifeq ($(_PENTIUM_M),Y)
        ifeq ($(_IFCV8),Y)
          FOPTIMIZE +=  -tpp7 -xB    # this are for Pentium M (aka Centrino)
        else
          FOPTIMIZE +=  -tpp7 -xW    # this are for PentiumIV
        endif
      else
        FOPTIMIZE +=  -tpp7 -xW    # this are for PentiumIV
      endif
    endif
    DEFINES   += -DIFCLINUX
    ifeq ($(_IFCV8),Y)
      FOPTIMIZE += -ansi_alias-
    endif
  endif

  ifeq ($(CC),icc)
    COPTIONS   =   -mp1 -w -g -vec_report3
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
    ifeq ($(FC),xlf)
      FOPTIONS  = -q32  -qextname -qfixed -qnosave -qsmallstack  -qalign=4k
      FOPTIONS +=  -NQ40000 -NT80000 -NS2048 -qmaxmem=8192 -qsigtrap -qxlf77=leadzero
      FOPTIMIZE= -O3 -qstrict  -qarch=auto -qtune=auto
      ifdef RSQRT
        FOPTIMIZE  += -qfloat=rsqrt:fltint
      endif
      FDEBUG= -O2 -g
      EXPLICITF = TRUE
      DEFINES  +=   -DXLFLINUX
      CPP=/usr/bin/cpp  -P -C -traditional
      FCONVERT = $(CPP) $(CPPFLAGS) $< > $*.f
    else
      FOPTIONS   = -fno-second-underscore -fno-globals -Wno-globals
      FOPTIMIZE  = -g -O2
    endif
    ifeq ($(CC),xlc)
      COPTIONS  +=  -q32 -qlanglvl=extended
    else
      COPTIONS   = -Wall
      COPTIMIZE  = -g -O2
    endif
    LDOPTIONS += -v
  endif



      LINK.f = $(FC) $(FOPTIONS) $(LDFLAGS) 
ifeq ($(LINUXCPU),x86)
  ifeq ($(FC),pgf77)
   LDOPTIONS=-g
   EXTRA_LIBS += -lm
  else
    ifeq ($(_FC),ifc)
    ifeq ($(_IFCV8),Y)
      EXTRA_LIBS +=  
    else
      EXTRA_LIBS +=   -Vaxlib  
    endif
      ifeq ($(_CPU),i786)
        EXTRA_LIBS +=  -lsvml
      endif
      EXTRA_LIBS += #-static
    else
  LDOPTIONS = -Xlinker --export-dynamic 
#  LDOPTIONS = --Xlinker -O -Xlinker -static
      EXTRA_LIBS += -lm
    endif
  endif
endif
#EXTRA_LIBS +=-lefence # link against Electricfence
ifeq ($(LINUXCPU),ppc)
  EXTRA_LIBS += -lm
    ifeq ($(FC),xlf)
      LINK.f   = xlf_r -Wl,-Bstatic $(LDFLAGS) 
    endif
endif


CORE_LIBS += -llapack $(BLASOPT) -lblas


# end of Linux, Cygnus
endif

ifeq ($(NWCHEM_TARGET),LINUX64)
   ifeq ($(FC),g77)
      ifndef USE_INTEGER4
      integer4:
	@echo You must define USE_INTEGER4 to compile with g77 
	@echo on 64-bit architectures.
	@echo Please type
	@echo 
	@echo " make FC=g77 USE_INTEGER4=y"
	@echo 
	@exit 1
      endif
   endif
      _FC=noifc
      ifeq ($(FC),g77)
       _FC=g77
      endif
      ifeq ($(FC),g77-ssa)
       _FC=g77
      endif
       ifdef USE_INTEGER4
         ifneq ($(_FC),g77)
           FOPTIONS += -i4 
         endif
       else
         FOPTIONS += -i8
         DEFINES  += -DEXT_INT
       endif
  MAKEFLAGS = -j 1 --no-print-directory
     _CPU = $(shell uname -m  )
     CORE_SUBDIRS_EXTRA = blas lapack
     RANLIB = echo
     DEFINES   +=   -DLINUX -DLINUX64 -DPARALLEL_DIAG
     ifeq ($(_CPU),alpha)
# using COMPAQ/DEC compilers (EA 3/13/2000)
       FC  = fort
       CC  = ccc      
       FOPTIONS   += -assume no2underscore -align dcommons -check nooverflow -assume accuracy_sensitive -check nopower -check nounderflow  -noautomatic
       DEFINES   +=   -DLINUXALPHA
       FOPTIMIZE =  -O4  -tune host -arch host  -math_library fast
       FVECTORIZE = -fast -O5 -tune host -arch host

       ifdef USE_INTEGER4
         FOPTIONS +=  -fpe0
# needed: binutils 2.11 for -taso option with some more hacking on bfd/elf.c
         LINK.f = fort -Wl,-taso $(LDFLAGS)  
         CORE_LIBS += -lcxml
       else
         FOPTIONS += -fpe3
         LINK.f = fort $(LDFLAGS)  
# this creates a static executable
#  LINK.f = fort $(LDFLAGS)   -Wl,-Bstatic
         CORE_LIBS += -llapack $(BLASOPT) -lblas
       endif
     endif

    ifeq ($(_CPU),ia64)
# Itanium  
# only working decent compiler: efc 6.0  release 127
# Version 6.0 Beta, Build 20020215   
# g77 not working 
# i4 not working 
#
      FC=efc
      CC=gcc
      DEFINES   +=   -DLINUXIA64 
      COPTIMIZE = -O1

      ifeq ($(FC),efc)
       _IFCV81= $(shell efc -V  2>&1|egrep -v Inte|egrep -v efc |egrep 8|awk ' /8\.1/  {print "Y"}')
       _IFCV8= $(shell efc -V  2>&1|egrep -v Inte|egrep -v efc |egrep 8|awk ' /8\./  {print "Y"}')
       ifeq ($(_IFCV8),Y)
         DEFINES+= -DIFCV8
         FOPTIONS += -quiet
       endif	
       ifeq ($(_IFCV81),Y)
         DEFINES+= -DIFCV81
       endif	
        ITANIUMNO = $(shell   cat /proc/cpuinfo | egrep family | head -1  2>&1 | awk ' /Itanium 2/ { print "-tpp2"; exit };/Itanium/ { print "-tpp1"}')
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
        FVECTORIZE =  -noalign -O3 -pad  -mP2OPT_hlo_level=2  -tpp1
        FOPTIMIZE =  -O3 -pad -mP2OPT_hlo_level=2  #-mP2OPT_align_array_to_cache_line=TRUE 
       ifeq ($(_IFCV81),Y)
#         FOPTIMIZE+= -IPF_fp_relaxed # breaks nwdft/xc/xc_pw91lda
       endif
        ifeq ($(_IFCV8),Y)
         EXTRA_LIBS += -quiet
         FDEBUG = -g -O2
        else
         EXTRA_LIBS += -Vaxlib 
         FDEBUG = -g -O2
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
        LINK.f = efc   -Qoption,link,-v  $(LDFLAGS)  
      endif
      ifeq ($(FC),g77)
        FOPTIONS  += -Wno-globals
        FOPTIONS  += -fno-globals -Wunused -fno-silent  
        FOPTIONS  += -fno-second-underscore  -g
        LDOPTIONS =   -Wl,--relax  -Wl,-Bstatic  
      endif
      ifeq ($(CC),ecc)
        COPTIONS   =   -ftz
        COPTIMIZE =  -O3 -hlo   -mP2OPT_hlo_level=2  
      endif
      ifeq ($(CC),gcc)
        COPTIONS   =   -O3 -funroll-loops -ffast-math
      endif

     CORE_LIBS +=  $(BLASOPT) -llapack -lblas
endif # end of ia32 bit
    ifeq ($(_CPU),x86_64)
# USE_INTEGER4=y need for g77
# make FC=g77 CC=gcc USE_INTEGER4=y USE_GCC34=y
#
      MAKEFLAGS = -j 2 --no-print-directory
      COPTIMIZE = -O1
      ifeq ($(FC),f77)
        defineFC: 
	@echo 
	@echo please set FC equal to g77, pathf90 or pgf90
	@echo e.g. type
	@echo 
	@echo "    make FC=pathf90"
	@echo 
	@exit 1
      endif
      ifeq ($(FC),pgf90)
        _FC=pgf90
      endif
      ifeq ($(FC),pgf77)
        _FC=pgf90
      endif
      ifeq ($(FC),ifc)
       _FC=ifc
      endif
      ifeq ($(FC),ifort)
       _FC=ifc
      endif
      ifeq ($(_FC),ifc)
       _IFCV81= $(shell ifc -v  2>&1|egrep 8|awk ' /8\.1/  {print "Y"}')
       ifeq ($(_IFCV81),Y)
# to get EM64T
# Intel 8.1 is required
       else
           @echo ifort 8.1 is required for x86_64 CPUs
           @exit 1
       endif
        FOPTIONS += -align -w -g -vec_report3
        DEFINES+= -DIFCV8 -DIFCLINUX
        ifeq ($(FC),ifc)
          FOPTIONS += -quiet
        endif
        FDEBUG= -O2 -g
        FOPTIMIZE = -O3 -prefetch  -unroll 
        FOPTIMIZE +=  -tpp7 -xW -ip
      endif	
      
      USE_LIB64 = y #for python linking

      ifeq ($(_FC),pgf90)
        FOPTIONS   +=    -Mrecursive -Mdalign -Mllalign -Kieee 
        FOPTIONS   +=    -tp k8-64  
#        FOPTIONS   +=    -Ktrap=fp
        FOPTIMIZE   =  -O3 -fastsse -Mnounroll -Minfo=loop -Mipa=fast
        FVECTORIZE   = -fast  -fastsse  -O3   -Mipa=fast
        FDEBUG = -g -O0 
        DEFINES  += -DCHKUNDFLW -DPGLINUX
      endif
      ifeq ($(FC),pathf90)
#pathscale 1.3 compiler
# tested Sep 30 2004 on RH AW3
        FOPTIONS   += -cpp -Wp,-P
        FOPTIONS   += -fno-second-underscore -fixedform
        FOPTIONS   += -align64
        FOPTIMIZE   = -O3 -OPT:Ofast:IEEE_arith=1:IEEE_NaN_inf=ON:Olimit=12000:ro=2:fold_reassociate=OFF#:div_split=OFF:fast_nint=OFF
        FVECTORIZE  = -O3 -OPT:Ofast -fno-math-errno
        DEFINES  += -DCHKUNDFLW -DPSCALE
        FDEBUG = -g -O2
        LDOPTIONS = -Wl,--warn-once   -Wl,--relax
      endif
      ifeq ($(_FC),g77)
        FOPTIONS  +=  -fno-globals -Wno-globals 
        DEFINES  +=   -DBAD_GACCESS
        FOPTIONS  += -Wunused  -fno-silent
        FOPTIONS  += -fno-second-underscore  -fno-f90 
        FOPTIMIZE  += -O3 -ffast-math 
        FOPTIMIZE  += -fstrength-reduce -fno-move-all-movables -fno-reduce-all-givs
        FOPTIMIZE  += -finline-functions -fschedule-insns2 
        FOPTIMIZE  += -mfpmath=sse -fstrict-aliasing  #-fprefetch-loop-arrays
        FOPTIMIZE  +=  -fexpensive-optimizations  
        ifdef USE_GCC34
          FOPTIMIZE  += -march=k8 -mtune=k8
          FOPTIMIZE  +=  -fold-unroll-all-loops
        else
          FOPTIONS  +=  -Wno-globals  
          FOPTIMIZE  +=  -funroll-all-loops
        endif
        FDEBUG = -Os -g
      endif
      COPTIONS   =   -O3 -funroll-loops -ffast-math  
      ifdef USE_GCC34
        COPTIONS  +=   -march=k8 -mtune=k8
      endif
     CORE_LIBS +=  $(BLASOPT) -llapack -lblas
     ifeq ($(BUILDING_PYTHON),python)
     EXTRA_LIBS += -lz  -lreadline -lncurses -lnwcutil  -lpthread -lutil -ldl
     endif
     ifdef  USE_GPROF
       FOPTIONS += -pg
       COPTIONS += -pg
       LDOPTIONS += -pg
     endif
endif

    ifeq ($(_CPU),ppc64)
      FC=xlf
      CC=/opt/cross/bin/powerpc64-linux-gcc
      ifeq ($(FC),xlf)
        FOPTIONS  =  -q64 -qextname -qfixed -qnosave  -qalign=natural
        FOPTIONS +=  -NQ40000 -NT80000 -qmaxmem=8192 -qxlf77=leadzero
        ifdef  USE_GPROF
          FOPTIONS += -pg
          LDOPTIONS += -pg
        endif
        FOPTIMIZE= -O3 -qstrict  -qarch=auto -qtune=auto -qcache=auto
        ifdef RSQRT
          FOPTIMIZE  += -qfloat=rsqrt:fltint
        endif
        FDEBUG= -O2 -g
        EXPLICITF = TRUE
        FCONVERT = $(CPP) $(CPPFLAGS) $< > $*.f
        DEFINES  +=   -DXLFLINUX -DCHKUNDFLW
        CPP=/usr/bin/cpp  -P -C -traditional
#    LDOPTIONS = -v #-F/home/eapra/nwchem/src/xlf.cfg:edo
#ld with SLES 8 broken. Needed snapshot from binutils
     ifdef USE_GPROF
LINK.f = /home/eapra/bin/ld   -L/home/eapra/nwchem/lib/LINUX64 -L/home/eapra/nwchem/src/tools/lib/LINUX64 /opt/cross/lib/gcc-lib/powerpc64-linux/3.2/../../../../powerpc64-linux/lib/gcrt1.o -L/opt/ibmcmp/xlf/8.1/lib64 -L/opt/ibmcmp/xlsmp/1.3/lib64 -L/opt/ibmcmp/xlf/8.1/lib64 -R/opt/ibmcmp/xlf/8.1/../../lib64 -R/opt/ibmcmp/xlf/8.1/../../lib64 -R/opt/ibmcmp/xlsmp/1.3/../../lib -L/opt/cross/lib/gcc-lib/powerpc64-linux/3.2/../../../../powerpc64-linux/lib -L/opt/cross/lib/gcc-lib/powerpc64-linux/3.2 -R/opt/cross/lib/gcc-lib/powerpc64-linux/3.2/../../../../powerpc64-linux/lib -R/opt/cross/lib/gcc-lib/powerpc64-linux/3.2 $(LDFLAGS) 
      else
LINK.f = /home/eapra/bin/ld  -L/home/eapra/nwchem/src/tools/lib/LINUX64 /opt/cross/powerpc64-linux/lib/crt1.o -L/opt/ibmcmp/xlf/8.1/lib64 -L/opt/ibmcmp/xlsmp/1.3/lib64 -L/opt/ibmcmp/xlf/8.1/lib64 -R/opt/ibmcmp/xlf/8.1/../../lib64 -R/opt/ibmcmp/xlf/8.1/../../lib64 -R/opt/ibmcmp/xlsmp/1.3/../../lib -L/opt/cross/powerpc64-linux/lib -L/opt/cross/lib/gcc-lib/powerpc64-linux/3.2 -R/opt/cross/powerpc64-linux/lib -R/opt/cross/lib/gcc-lib/powerpc64-linux/3.2 $(LDFLAGS)  -lpthread
      endif
        ifdef USE_INTEGER4
          FOPTIONS += -qintsize=4
        else
          FOPTIONS += -qintsize=8
        endif
      endif
      ifeq ($(CC),xlc)
        COPTIONS  +=  -q64 -qlanglvl=extended
      endif
     CORE_LIBS +=  $(BLASOPT) -llapack -lblas
     EXTRA_LIBS +=  -dynamic-linker /lib64/ld64.so.1 -melf64ppc -lxlf90_r -lxlopt -lxlomp_ser -lxl -lxlfmath -ldl -lm -lc -lgcc -lm
    endif

ifeq ($(BUILDING_PYTHON),python)
#   EXTRA_LIBS += -ltk -ltcl -L/usr/X11R6/lib -lX11 -ldl
   EXTRA_LIBS += -lpthread -ldl
endif

endif


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

        CORE_LIBS = -lnwcutil \
                    -L$(GA_LIBDIR) -lglobal -lpario -lma -lpeigs \
                    -ltcgmsg-mpi -L/usr/lang/mpi2/lib32 -lmpi -lmp
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

        CORE_LIBS = -lnwcutil \
                    -L$(GA_LIBDIR) -lglobal -lpeigs -lpario -lma \
                    -ltcgmsg-mpi -L/usr/lang/mpi2/lib64 -lmpi -lmp
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
   FOPTIMIZE = -O scalar3,aggress,unroll2,vector3
      FDEBUG = -O scalar1,vector1
   COPTIMIZE = -O -h inline2  -h aggress


       EXTRA_LIBS = -lsci64 -llapack  -lblas  # need make dbl_to_sngl for this
#       EXTRA_LIBS =  -llapack  -lblas 

#      EXPLICITF     = TRUE
      FCONVERT      = $(CPP) $(CPPFLAGS)  $< | sed '/^\#/D'  > $*.f
#end of sv2
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
#-do not use#   CORE_LIBS = -lnwcutil -lpario -lglobal -lma -lpeigs -llapack -lblas
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
ifdef USE_LIB64
CORE_LIBS += -L$(PYTHONHOME)/lib64/python$(PYTHONVERSION)/config -lpython$(PYTHONVERSION)
else
CORE_LIBS += $(PYTHONHOME)/lib/python$(PYTHONVERSION)/config/libpython$(PYTHONVERSION).a
endif
endif
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


###################################################################
#  All machine dependent sections should be above here, otherwise #
#  some of the definitions below will be 'lost'                   #
###################################################################
#the new GA uses ARMCI library
      CORE_LIBS += -larmci

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

(%.o):  %.o

# Preceding line has a tab to make an empty rule

# a .F.f rule is needed for any system target where the default .F.f rule does not work
# AND the EXPLICITF is not already true.  Right now this is only LINUX with g77
      ifneq ($(FC),xlf)
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


