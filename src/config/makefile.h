
# $Id: makefile.h,v 1.71 1994-10-11 19:35:44 d3g681 Exp $

# Common definitions for all makefiles ... these can be overridden
# either in each makefile by putting additional definitions below the
# include statement, or on the command line


#
# TOPDIR points to your top-level directory that contains
# src, lib, config, ... (SRCDIR, etc., are derived from TOPDIR)
# Do a setenv for NWCHEM_TOP to be the top level directory
#

ifndef NWCHEM_TOP
error1:
	@echo You must define NWCHEM_TOP in your environment to be the path
	@echo of the top level nwchem directory ... something like
	@echo     setenv NWCHEM_TOP /msrc/home/bill_clinton/nwchem
	@exit 1
endif

     TOPDIR = $(NWCHEM_TOP)
     SRCDIR = $(TOPDIR)/src
     LIBDIR = $(TOPDIR)/lib
     BINDIR = $(TOPDIR)/bin
     INCDIR = $(TOPDIR)/src/include
     CNFDIR = $(TOPDIR)/src/config

#
# Do a setenv for NWCHEM_TARGET to be the machine you wish to build for
# (one of SUN, DELTA, IBM, KSR, PARAGON)
#

ifndef NWCHEM_TARGET
error2:
	@echo You must define NWCHEM_TARGET in your environment to be the name
	@echo of the machine you wish to build for ... for example
	@echo     setenv NWCHEM_TARGET SUN
	@echo Known targets are SUN, DELTA, KSR, PARAGON, IBM
	@exit 1
endif

     TARGET = $(NWCHEM_TARGET)

#
# Define NWSUBDIRS to be list of subdirectories of SRC to be made
#
# The include directory should be first so that the include
# files are all present and correct before any compilation
#

    NWSUBDIRS = include ddscf NWints develop global db rtdb basis inp util \
                moints atomscf geom input ma tcgmsg gradients rimp2 riscf \
                stepper pstat $(SUBDIRS_EXTRA)

#
# Define LIBPATH to be paths for libraries that you are linking in
# from precompiled sources and are not building now. These libraries
# will be searched AFTER anything you are building now.
# e.g. LIBPATH = -L/msrc/proj/mss/lib
#
    LIBPATH = 

#
# Define INCPATH to be directories to get includes for
# libraries that you are not building now.  These directories
# will be searched AFTER anything you are building now.
#
    INCPATH = 


##########################################################
#                                                        #
# Should NOT need to modify below here unless porting to #
# a new machine or changing compiler options             #
#                                                        #
##########################################################

# Only the SUN, KSR, PARAGON, DELTA and IBM versions are up to date

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
#              for parallel builds. -s says be quiet about changing directory.
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
#  LDOPTIONS = additional options to be passed to the linker (LDFLAGS is
#              built from this and the library path info).  LDOPTIONS is
#              the best way to add to the link command.
#
#    ARFLAGS = options for AR (ru)
#
#    DEFINES = C preprocessor defines for both C and Fortran
#
#       LIBS = List of libraries and paths for libraries in addition 
#              to the LIBDIR and LIBPATH options.
#
# SUBDIRS_EXTRA = List of additional directories (e.g., BLAS) that
#                 are needed on this machine.
#
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


#
# Establish some required defaults which may need overriding
# for some machines

      SHELL = /bin/sh
    ARFLAGS = ru
     FDEBUG = -g
     CDEBUG = -g

#
# Machine specific stuff
#

ifeq ($(TARGET),SUN)
#
# Sun running SunOS
#

    SUBDIRS_EXTRA = blas lapack
         CC = gcc
     RANLIB = ranlib
  MAKEFLAGS = -j2
    INSTALL = @echo $@ is built

   FOPTIONS = -Nl99
   COPTIONS = -Wall
  FOPTIMIZE = -O3
  COPTIMIZE = -g -O2

    DEFINES = -DSUN

       LIBS = -ltest -lddscf -lriscf -lrimp2 -lgradients -lnwints \
              -lstepper -lmoints \
              -lguess -lglobal -lutil \
	      -ltcgmsg -llapack -lblas

  EXPLICITF = FALSE
endif

ifeq ($(TARGET),KSR)
#
# KSR running OSF
#
    SUBDIRS_EXTRA = blas
     RANLIB = @echo
  MAKEFLAGS = -j20
    INSTALL = @echo $@ is built

   FOPTIONS = -r8
   COPTIONS = 
  FOPTIMIZE = -xfpu3 -O1
  COPTIMIZE = -xfpu3 -O1

    DEFINES = -DKSR -DPARALLEL_DIAG -DLongInteger

       LIBS = -L/home/d3g681/TCGMSG_DISTRIB \
              -ltest -lddscf -lriscf -lrimp2 -lgradients -lnwints\
              -lstepper -lmoints \
              -lguess -lglobal -lutil \
	      -lpeigs \
              -lksrlapk -lksrblas -llapack2 -lblas2  -ltcgmsg -para -lrpc

  EXPLICITF = FALSE
endif

ifeq ($(TARGET),PARAGON)
#
# Intel Paragon running OSF
# (native build, or cross-compiled -- the latter is faster in most cases)
#
    SUBDIRS_EXTRA = lapack

         FC = if77
         CC = icc
         AR = ar860
     RANLIB = @echo

  MAKEFLAGS = -j4 
    INSTALL = @echo $@ is built

  FOPTIONS = -Knoieee
  COPTIONS = -Knoieee
 FOPTIMIZE = -O2 -Minline=1000
FVECTORIZE = -O2 -Minline=1000 -Mvect
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
       LIBS = -L/home/delilah11/gifann/lib \
              -ltest -lddscf -lriscf -lrimp2 -lgradients -lnwints \
	      -lstepper -lmoints \
              -lguess -lglobal -lutil \
	      -lpeigs_paragon -ltcgmsg -llapack $(LIBDIR)/liblapack.a \
              -lkmath -nx

  EXPLICITF = FALSE
endif

ifeq ($(TARGET),DELTA)
#
# DELTA/IPSC running NX
#
    SUBDIRS_EXTRA = lapack
# blas
        FC = if77
        CC = icc
       CPP = /usr/lib/cpp
        AR = ar860
    RANLIB = @echo

   INSTALL = "strip860 nwchem; rcp nwchem delta1: &"
 MAKEFLAGS = -j2 

  FOPTIONS = -Knoieee
  COPTIONS = -Knoieee
 FOPTIMIZE = -O2 -Minline=1000
FVECTORIZE = -O2 -Miniline=1000 -Mvect
 COPTIMIZE = -O2

   DEFINES = -DNX -DDELTA -DIPSC -DNO_BCOPY  -D__IPSC__ -DPARALLEL_DIAG
      LIBS = -L/home/delilah11/gifann/lib \
             -ltest -lddscf -lriscf -lrimp2 -lgradients -lnwints \
	     -lstepper -lmoints \
             -lguess -lglobal -lutil \
             -lglobal -lpeigs_delta -ltcgmsg $(LIBDIR)/liblapack.a -llapack -lkmath -node

 EXPLICITF = FALSE
endif


ifeq ($(TARGET),SGITFP)
#
# SGI power challenge
#
# SUBDIRS_EXTRA are those machine specific libraries required 

    SUBDIRS_EXTRA = blas lapack
         FC = f77
         CC = cc
         AR = ar
     RANLIB = echo

    INSTALL = @echo nwchem is built
  MAKEFLAGS = -j 4

  FOPTIONS = -d8 -i8 -mips4 -64 -r8 -G 0 -OPT:roundoff=3:IEEE_arithmetic=3
  COPTIONS = -fullwarn -mips4 
 FOPTIMIZE = -O3 -OPT:fold_arith_limit=4000 -TENV:X=3
FVECTORIZE = -O3 -OPT:fold_arith_limit=4000 -TENV:X=3 -WK,-so=1,-o=1

 COPTIMIZE = -O

    DEFINES = -DSGITFP -DSGI -DLongInteger
       LIBS = -L$(LIBDIR) $(LIBPATH) \
              -ltest -lddscf -lriscf -lrimp2 -lgradients -lnwints \
	      -lstepper -lmoints \
              -lguess -lutil -lglobal\
	      -ltcgmsg -llapack -lblas

  EXPLICITF = FALSE
endif


ifeq ($(TARGET),SGI)
#
# SGI normal
#
# SUBDIRS_EXTRA are those machine specific libraries required 

    SUBDIRS_EXTRA = blas lapack
         FC = f77
         CC = cc
         AR = ar
     RANLIB = echo

    INSTALL = @echo nwchem is built
  MAKEFLAGS = -j 4

  FOPTIONS = -mips2
  COPTIONS = -mips2 -fullwarn
 FOPTIMIZE = -O2
 COPTIMIZE = -O2

    DEFINES = -DSGI 
       LIBS = -L$(LIBDIR) $(LIBPATH) \
              -ltest -lddscf -lriscf -lrimp2 -lgradients -lnwints \
	      -lstepper -lmoints \
              -lguess -lutil -lglobal\
	      -ltcgmsg -llapack -lblas -lmalloc 

  EXPLICITF = FALSE
endif

ifeq ($(TARGET),IBM)
#
# IBM AIX . tested rak 4/94
# note: using only source blas.
#

    SUBDIRS_EXTRA = lapack blas
         FC = xlf
    ARFLAGS = urs
     RANLIB = echo
  MAKEFLAGS = -j2
    INSTALL = @echo $@ is built
        CPP = /usr/lib/cpp -P

   FOPTIONS = -qEXTNAME
   COPTIONS =
  FOPTIMIZE = -O3
  COPTIMIZE = -O

    DEFINES = -DIBM -DEXTNAM

       LIBS = -L$(LIBDIR) $(LIBPATH) \
              -ltest -lddscf -lriscf -lrimp2 -lgradients -lnwints \
	      -lstepper -lmoints \
              -lguess -lglobal -lutil \
	      -ltcgmsg -llapack -lblas

 EXPLICITF = TRUE
#
endif



###################################################################
#  All machine dependent sections should be above here, otherwise #
#  some of the definitions below will be 'lost'                   #
###################################################################

ifdef OPTIMIZE
    FFLAGS = $(FOPTIMIZE) $(FOPTIONS) 
    CFLAGS = $(COPTIMIZE) $(COPTIONS) 
else
    FFLAGS = $(FDEBUG) $(FOPTIONS) 
    CFLAGS = $(CDEBUG) $(COPTIONS) 
endif
  INCLUDES = -I. $(LIB_INCLUDES) -I$(INCDIR) $(INCPATH)
  CPPFLAGS = $(INCLUDES) $(DEFINES) $(LIB_DEFINES)
   LDFLAGS = $(LDOPTIONS) -L$(LIBDIR) $(LIBPATH)


#
# Define known suffixes mostly so that .p files don't cause pc to be invoked
#

.SUFFIXES:	
.SUFFIXES:	.o .s .F .f .c

ifeq ($(EXPLICITF),TRUE)
#
# Needed on machines where FCC does not preprocess .F files
# with CPP to get .f files
#
.SUFFIXES:	
.SUFFIXES:	.o .s .F .f .c

.F.o:	
	$(MAKE) $*.f
	$(FC) -c $(FFLAGS) $*.f
	/bin/rm -f $*.f

.F.f:	
ifeq ($(TARGET),IBM)
	$(CPP) $(CPPFLAGS) $*.F > $*.f
else
	$(CPP) $(CPPFLAGS) < $*.F | sed '/^#/D' | sed '/^[a-zA-Z].*:$/D' > $*.f
endif

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $*.c
endif




