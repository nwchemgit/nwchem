##############################################################################
#
# include makefile for nmake under Windows NT
#
##############################################################################

#
#	$Id$
#

#
# If NWCHEM_TOP_WIN32 is set, take it.  Otherwise, if NWCHEM_TOP is set,
# use that.  One of the two must be set, however.
#
!IFNDEF NWCHEM_TOP_WIN32
!IFDEF NWCHEM_TOP
NWCHEM_TOP_WIN32 = $(NWCHEM_TOP)
!ELSE
!MESSAGE You must define NWCHEM_TOP in your environment to be the path
!MESSAGE of the top level nwchem directory in DOS format ... e.g.
!MESSAGE     NWCHEM_TOP="D:\PNNL\nwchem"
!MESSAGE
!ERROR NWCHEM_TOP not set.
!ENDIF
!ENDIF

TOPDIR = $(NWCHEM_TOP_WIN32)
SRCDIR = $(TOPDIR)\src
# Set LIB_DIR externally to override library dir name under $(TOPDIR)\lib
!IFDEF LIB_DIR
LIBDIR = $(TOPDIR)\lib\$(LIB_DIR)
!ELSE
LIBDIR = $(TOPDIR)\lib\win32
!ENDIF
# !!! This is called LIB_DISTRIB in prev NT makefiles
LIB_DISTRIB = $(LIBDIR)
BINDIR = $(TOPDIR)\bin\win32
INCDIR = $(TOPDIR)\src\include
CNFDIR = $(TOPDIR)\src\config
OBJDIR = .\obj

GLOB_DEFINES =-DWIN32 -DUSE_FCD
DEFINES = $(GLOB_DEFINES) $(LIB_DEFINES)
GLOB_INCLUDES= -I"$(SRCDIR)\include" -I"$(SRCDIR)\tools\include"
INCLUDES = $(GLOB_INCLUDES) $(LIB_INCLUDES)

#AR = link.exe -lib -nologo
AR = lib -nologo
ARFLAGS = /out:$(LIBRARY_PATH)

CC = cl -nologo
FC = f90 -nologo

!IFDEF NWDEBUG
COPT = -Z7
FOPT = /debug:full /nooptimize
!ELSE
COPT = 
FOPT = /fast /optimize:5 /noinline /nofltconsistency
# Added /noinline since it breaks LAPACK dlamach routines
!ENDIF

CFLAGS = -W3 $(COPT) $(INCLUDES) $(DEFINES) -Fo"$(OBJDIR)/" -c
FFLAGS = $(FOPT) $(INCLUDES) $(DEFINES)  /automatic /check:none /traceback /fpscomp=nogeneral /warn:argument_checking /warn:nofileopt /warn:nouncalled /object:"$(OBJDIR)/" /fpp:"/c /m" /nodefine /nokeep -c

.SUFFIXES:
.SUFFIXES:      .obj .s .F .c

.c{$(OBJDIR)}.obj:
	$(CC) $(CFLAGS) $<

.F{$(OBJDIR)}.obj:
	$(FC) $(FFLAGS) $<
