# Makefile for movecs <--> ASCII translators
# $Id: Makefile,v 1.5 2002-03-13 18:38:01 windus Exp $

# Tested on KSR, Cray-T3D
# Requires NWCHEM_TARGET to be set.
# Requires NWChem includes and library for MA (in libutil.a)

ifeq ($(NWCHEM_TARGET),CRAY-T3D) # (must run on MPP not host)
   FFLAGS = -Wf"-dp" -Wp"-F" -Ccray-t3d 
endif

ifeq ($(NWCHEM_TARGET),CRAY-T3E)
   FFLAGS = -O 1 -d p -F
   FC     = f90
endif

ifeq ($(NWCHEM_TARGET),LINUX)
   FFLAGS = -fno-second-underscore
endif

ifeq ($(NWCHEM_TARGET),KSR)
   FFLAGS = -r8
endif

ifeq ($(NWCHEM_TARGET),IBM)
   FFLAGS = -qEXTNAME
endif

ifeq ($(NWCHEM_TARGET),IBM)
   CPPFLAGS = -WF,"-I../../src/tools/include"		# Needs MA includes
else
   CPPFLAGS = -I../../src/tools/include			# Needs MA includes
endif
LDFLAGS = -L../../src/tools/lib/${NWCHEM_TARGET} -lma	# Needs MA library

all:    mov2asc asc2mov

mov2asc:	mov2asc.o
	${FC} ${FFLAGS} ${CPPFLAGS} -o $@ $^ ${LDFLAGS} 

asc2mov:	asc2mov.o
	${FC} ${FFLAGS} ${CPPFLAGS} -o $@ $^ ${LDFLAGS} 


clean realclean:
	$(RM) *.o mov2asc asc2mov



%.o:	%.F
	$(FC) -c $(FFLAGS) $(CPPFLAGS) $<

%.o:	%.f
	$(FC) -c $(FFLAGS) $<

%.o:	%.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $<


