# Makefile for movecs <--> ASCII translators
# $Id: Makefile,v 1.2 1996-12-09 20:03:14 gg502 Exp $

# Tested on KSR, Cray-T3D
# Requires NWCHEM_TARGET to be set.
# Requires NWChem includes and library for MA (in libutil.a)

ifeq ($(NWCHEM_TARGET),CRAY-T3D) # (must run on MPP not host)
   FFLAGS = -Wf"-dp" -Wp"-F" -Ccray-t3d 
endif

ifeq ($(NWCHEM_TARGET),KSR)
   FFLAGS = -r8
endif

ifeq ($(NWCHEM_TARGET),IBM)
   FFLAGS = -qEXTNAME
endif

ifeq ($(NWCHEM_TARGET),IBM)
   CPPFLAGS = -WF,"-I../../src/include"			# Needs MA includes
else
   CPPFLAGS = -I../../src/include			# Needs MA includes
endif
LDFLAGS = -L../../lib/${NWCHEM_TARGET} -lutil	# Needs MA library

all:    mov2asc asc2mov

mov2asc:	mov2asc.o
	${FC} ${FFLAGS} ${CPPFLAGS} -o $@ $^ ${LDFLAGS} 

asc2mov:	asc2mov.o
	${FC} ${FFLAGS} ${CPPFLAGS} -o $@ $^ ${LDFLAGS} 


clean realclean:
	$(RM) *.o mov2asc asc2mov
