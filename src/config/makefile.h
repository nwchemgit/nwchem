
# $Id: makefile.h,v 1.9 1994-04-07 18:44:10 d3e129 Exp $

# Common definitions for all makefiles ... these can be overridden
# either in each makefile by putting additional definitions below the
# include statement, or on the command line


#
# TOPDIR points to your top-level directory that contains
# src, lib, config, ... (SRCDIR, etc., are derived from TOPDIR)
#
# Either do a setenv for NWCHEMTOP or define NWCHEMTOP here
# ... it is preferable to do the setenv then this file is indep of
# who is using it.
#
# NWCHEMTOP = /msrc/home/d3g681
     TOPDIR = $(NWCHEMTOP)
     SRCDIR = $(TOPDIR)/src
     LIBDIR = $(TOPDIR)/lib
     BINDIR = $(TOPDIR)/bin
     INCDIR = $(TOPDIR)/src/include

#
# Define TARGET to be the machine you wish to build for
# (one of SUN, IPSC, KSR)
#
# Either do a setenv for NWCHEM_TARGET or define it here
# ... it is preferable to do the setenv then this file is independent
# ... of who is using it!!!
#
# NWCHEM_TARGET = SUN
     TARGET = $(NWCHEM_TARGET)

#
# Define SUBDIRS to be list of subdirectories of SRC to be made
#
# The include directory should be first so that the include
# files are all present and correct before any compilation
#
    SUBDIRS = include develop global db NWints rtdb basis inp util \
              geom input ma tcgmsg
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

# !!!! Only the SUN version is up to date !!!!!

ifeq ($(TARGET),SUN)
#
# Sun running SunOS
#
         FC = f77
         CC = gcc
         AR = ar
     RANLIB = ranlib
      SHELL = /bin/sh
       MAKE = make
  MAKEFLAGS = -j 2
    INSTALL = echo $@ is built

       FOPT = -g -u -Nl99
   FOPT_REN = $(FOPT)
       COPT = -g
     FLDOPT = $(FOPT)
     CLDOPT = $(COPT)
  INCLUDES =  -I. $(LIB_INCLUDES) -I$(INCDIR) $(INCPATH)
   WARNINGS = -Wall
#-Wshadow -Wcast-qual -Wwrite-strings -Wpointer-arith
    DEFINES = -DSUN $(LIB_DEFINES) 
     FFLAGS = $(FOPT) $(INCLUDES) $(DEFINES)
     CFLAGS = $(COPT) $(INCLUDES) $(DEFINES) $(WARNINGS)
    ARFLAGS = rcv

       LIBS = -L$(LIBDIR) $(LIBPATH) \
              -ltest -lnwints \
              -linput -lgeom -lbasis -lutil -lglobal -lrtdb -ldb -linp \
	      -lutil -lma -ltcgmsg

  EXPLICITF = FALSE
endif

ifeq ($(TARGET),IPSC)
#
# DELTA/IPSC running NX
#
        FC = if77
        CC = icc
       CPP = /usr/lib/cpp
        AR = ar860

    RANLIB = echo
     SHELL = /bin/sh
   INSTALL = rcp $@ delta2:
      FOPT = -O2 -Knoieee -Mquad -node -Minline=100
  FOPT_REN = -O2 -Knoieee -Mquad -Mreentrant -Mrecursive -node
      COPT = -O2 -Knoieee -Mreentrant -node
  INCLUDES =  -I. -I$(SRCDIR)/rtdb -I$(SRCDIR)/global -I$(SRCDIR)/tcgmsg -I$(SRCDIR)/NWints \
              -I$(SRCDIR)/util -I$(SRCDIR)/ma -I$(SRCDIR)/db -I$(SRCDIR)/tcgmsg/ipcv4.0
   DEFINES = -DNX -DIPSC -DNO_BCOPY  $(LIB_DEFINES)
#  -DGA_TRACE
    FFLAGS = $(FOPT)
    CFLAGS = $(COPT) $(INCLUDES) $(DEFINES)
 MAKEFLAGS = -j 2
    FLDOPT = $(FOPT) -node
    CLDOPT = $(COPT) -node
   ARFLAGS = rcv
      LIBS = $(SRCDIR)/input/libinput.a \
             $(SRCDIR)/ddscf/libddscf.a \
             $(SRCDIR)/NWints/libnwints.a  \
             $(SRCDIR)/rtdb/librtdb.a \
             $(SRCDIR)/db/libdb.a \
             $(SRCDIR)/global/libglobal.a \
             $(SRCDIR)/trace/libtrace.a \
             $(SRCDIR)/tcgmsg/ipcv4.0/libtcgmsg.a \
             $(SRCDIR)/util/libutil.a \
             $(SRCDIR)/ma/libma.a \
             $(SRCDIR)/peigs1.0/libpeigs.a \
             $(SRCDIR)/peigs1.0/liblapack.a \
             -lkmath 

 EXPLICITF = TRUE
endif


ifeq ($(TARGET),IBM)
#
# IBM AIX .... NOT YET TESTED !!!!!
#
#        FC = xlf
#        CC = xlc
#        AR = ar
#    RANLIB = ranlib
#   INSTALL = echo
#     SHELL = /bin/sh
#      FOPT = -g
#      COPT = -g
#  INCLUDES = -I. -I../ma
#   DEFINES = -DTCGMSG
#    FFLAGS = -qEXTNAME $(FOPT)
#    FLDOPT = $(FOPT) -b rename:.exit_,.exit
#    CFLAGS = $(COPT) $(INCLUDES) $(DEFINES)
#    CLDOPT = $(COPT)
#   ARFLAGS = rcv
#      LIBS = ../tcgmsg/ipcv4.0/libtcgmsg.a ../ma/libma.a -lc
# EXPLICITF = TRUE
#
endif

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
	$(CPP) $(INCLUDES) $(DEFINES) < $*.F | sed '/^#/D' > $*.f

.c.o:
	$(CC) $(CFLAGS) -c $*.c
endif
