# $Id: GNUmakefile,v 1.7 1995-01-29 09:34:53 d3e129 Exp $
       HEADERS = 
       OBJ_OPTIMIZE stpr_sjacobi.o
           OBJ = stpr_cheker.o     stpr_chkhess.o    stpr_clnstp.o  \
                 stpr_cneed.o      stpr_cross.o      stpr_driver.o  \
                 stpr_fndmde.o     stpr_formup.o     stpr_gmass.o   \
                 stpr_hesupd.o     stpr_hmass.o      stpr_inihes.o  \
                 stpr_face.o       stpr_locate.o     stpr_matout.o  \
                 stpr_mgeom.o      stpr_mkstep.o     stpr_mktri.o   \
                 stpr_output.o     stpr_partit.o     stpr_place.o   \
                 stpr_pmat.o       stpr_prntpd.o     stpr_projec.o  \
                 stpr_putbck.o     stpr_readsq.o     stpr_save.o    \
                 stpr_setmass.o    stpr_stepcor.o \
                 stpr_stepper.o    stpr_sumstc.o     stpr_tprntd.o  \
                 stpr_uprcase.o    stpr_wlkchk.o     stpr_wrtsq.o   \
                 stpr_walk.o       scf_optimize.o    stpr_datestamp.o \
                 stpr_walk_dft.o distang.o scf_optimize.o  dft_optimize.o   
       LIBRARY = libstepper.a
   LIB_TARGETS = 

     USES_BLAS = 

include ../config/makefile.h
include ../config/makelib.h

stpr_cheker.F:	chow.fh
stpr_driver.F:	chow.fh
stpr_hesupd.F:	chow.fh
stpr_inihes.F:	chow.fh
stpr_face.F:	chow.fh
stpr_mkstep.F:	chow.fh
stpr_partit.F:	chow.fh
stpr_save.F:	chow.fh
stpr_sumstc.F:	chow.fh




