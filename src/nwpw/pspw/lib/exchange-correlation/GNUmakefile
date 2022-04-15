#
# $Id$
#
SUBDIRS = vdw-DF

 OBJ_OPTIMIZE = \
            	v_bwexc_all.o \
		v_exc.o v_dirac.o \
		v_mexc.o \
            	v_bwexc.o \
            	v_bwexc_euv.o \
		nwpw_meta_gga.o \
            	pbe96.o \
            	pbesol.o \
            	revpbe.o \
            	blyp.o \
            	b3lyp.o \
		bnl2007_x.o\
		hsepbe.o \
		beef_gga.o \
		pw92c.o \
		vs98.o \
		tpss03.o scan.o r2scan.o pw91c.o  pbe96x.o pbe96c.o \
		m06.o  m06-l.o 

#            vdirac.o \

	  OBJ = 

      LIBRARY = libpspw.a

    USES_BLAS = v_bwexc_euv.F v_bwexc.F nwpw_meta_gga.F v_mexc.F


  LIB_DEFINES = 


include ../../../../config/makefile.h
include ../../../../config/makelib.h

