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
#  -- PFFT wrapper routine Makefile --
#     Pacific Northwest Laboratory
#     Jan, 1997
#
#======================================================================

# fft routine: EMSLFFT or ESSLFFT
# parallel layer: TCGPROCS or MPIPROCS
PROCTYPE  = TCGPROCS
#PROCTYPE = MPIPROCS

FFTTYPE   = EMSLFFT
#FFTTYPE  = ESSLFFT
#-------------------------------------------

ifeq ($(PROCTYPE),TCGPROCS)
  ifeq ($(FFTTYPE),EMSLFFT)
     LIB_DEFINES = -DEMSLFFT -DTCGPROCS -DFFTTIME
  else
     LIB_DEFINES = -DESSLFFT -DTCGPROCS -DFFTTIME
  endif
  OBJ = testfft3.o 
else
  ifeq ($(FFTTYPE),EMSLFFT)
     LIB_DEFINES = -DEMSLFFT -DMPIPROCS -DFFTTIME
  else
     LIB_DEFINES = -DESSLFFT -DMPIPROCS -DFFTTIME
  endif
  OBJ = testfft3.o o mympi.o mympic.o  
endif


DO_ESSL = YES
ifeq ($(NWCHEM_TARGET),SP)
  LIB_INC = -I../../include -I/usr/lpp/ppe.poe/include
  LIBDIRS = -L/usr/lpp/ppe.poe/lib -L/usr/lib -L$(HOME)/pfft1.0/lib  \
          -L../../../lib/SP
  ifeq ($(FFTTYPE),EMSLFFT)
     DO_ESSL = NO
  endif
else
  LIB_INC = -I../../include 
  LIBDIRS = -L../../../lib/$(TARGET)
  ifeq ($(FFTTYPE),EMSLFFT)
     DO_ESSL = NO
  endif
endif


include ../../config/makefile.h

WRAP = ../../fft/pfftwrap

test:	$(OBJ) 
	touch $(WRAP)/nwfft3d.F
	touch $(WRAP)/iessl.F
	(cd $(WRAP); $(MAKE) USE_ESSL='$(DO_ESSL)';)
	touch testfft3.F
	$(MAKE) testfft3.o
	$(LINK.f) -o testfft3 $(OBJ) -lnwctask $(LIBS)	


