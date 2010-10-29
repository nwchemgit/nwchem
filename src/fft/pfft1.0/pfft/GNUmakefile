#
# $Id$
#
#======================================================================
#
# DISCLAIMER
#
# This material was prepared as an account of work sponsored by an
# agency of the United States Government.  Neither the United States
# Government nor the United States Department of Energy, nor Battelle,
# nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
# ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
# COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
# SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
# INFRINGE PRIVATELY OWNED RIGHTS.
#
# ACKNOWLEDGMENT
#
# This software and its documentation were produced with Government
# support under Contract Number DE-AC06-76RLO-1830 awarded by the United
# States Department of Energy.  The Government retains a paid-up
# non-exclusive, irrevocable worldwide license to reproduce, prepare
# derivative works, perform publicly and display publicly by or for the
# Government, including the right to distribute to other Government
# contractors.
#
#======================================================================
#
#  -- PFFT routine (version 1.0) --
#     Pacific Northwest Laboratory
#     April 5, 1995
#
#======================================================================
#
# set FFTLIB = the library to link with for single processor 1-D FFT routines
# set FFTNAM = the name to define in CPP for compiling *.cpp files.
#              Used to generate correct calls to FFTLIB routines.

LIBRARY   = libpfft.a

IO_STYLE = FILE_IO
LIB_DEFINES = -D${IO_STYLE} -DIBM ${GENDEFS} -DTCGMSG -DFFTPACK

#CPP		= /usr/lib/cpp -P -C -D${COMM_PKG} -D${IO_STYLE} -DIBM ${GENDEFS}

###################################		dcopy.o 
OBJ_OPTIMIZE=	blktrn.o \
		drain.o \
		fillup.o \
		mapfil.o \
		mapchk.o \
		prlst2.o \
		xstop.o \
		pairup.o \
		gexit.o \
		igdiff.o \
		lstcnt.o \
		lstind.o \
		pfft.o \
		pfft3d.o \
		prfft3.o \
		prfft.o \
		prffti.o \
		chk3d.o \
		mfftc1.o \
		mfftr1.o \
		mcfftr.o \
		mrfftr.o \
		pfft1i.o \
		pfft1.o \
		wrkspc2.o

USES_BLAS = drain.F fillup.F
include ../../../config/makefile.h
include ../../../config/makelib.h

