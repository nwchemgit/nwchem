# $Id: GNUmakefile,v 1.26 1996-07-22 06:04:43 d3e129 Exp $
       HEADERS = 
  OBJ_OPTIMIZE = stpr_sjacobi.o stpr_partit.o 
           OBJ = stpr_cheker.o     stpr_chkhess.o    stpr_clnstp.o  \
                 stpr_cneed.o      stpr_cross.o      stpr_driver.o  \
                 stpr_fndmde.o     stpr_formup.o     stpr_gmass.o   \
                 stpr_hesupd.o     stpr_hmass.o      stpr_inihes.o  \
                 stpr_face.o       stpr_locate.o     stpr_matout.o  \
                 stpr_mgeom.o      stpr_mkstep.o     stpr_mktri.o   \
                 stpr_output.o     stpr_place.o   \
                 stpr_pmat.o       stpr_prntpd.o     stpr_projec.o  \
                 stpr_putbck.o     stpr_readsq.o     stpr_save.o    \
                 stpr_setmass.o    stpr_stepcor.o \
                 stpr_stepper.o    stpr_sumstc.o     stpr_tprntd.o  \
                 stpr_uprcase.o    stpr_wlkchk.o     stpr_wrtsq.o   \
                 stpr_walk.o       scf_optimize.o    stpr_datestp.o \
                 distang.o         dft_optimize.o    \
                 stpr_gah_at.o \
                 stpr_gh_at.o \
                 stpr_gah_opt.o \
                 stpr_gen_aux.o    stpr_filenames.o eopt.o

       LIBRARY = libstepper.a
   LIB_TARGETS = 

   USES_BLAS = eopt.F stpr_mkstep.F stpr_partit.F

include ../config/makefile.h
include ../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(stpr_stepcor.o):	$(INCDIR)/util.fh ./cstprfiles.fh 
$(LIBRARY_PATH)(stpr_place.o):	./cstprfiles.fh 
$(LIBRARY_PATH)(stpr_filenames.o):	./cstprfiles.fh 
$(LIBRARY_PATH)(stpr_gen_aux.o):	$(INCDIR)/tcgmsg.fh ./stpr_optimize.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgtypesf.h $(INCDIR)/msgids.fh $(INCDIR)/global.fh ./cstprfiles.fh $(INCDIR)/context.fh 
$(LIBRARY_PATH)(stpr_fndmde.o):	./cstprfiles.fh 
$(LIBRARY_PATH)(scf_optimize.o):	./stpr_optimize.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cstprfiles.fh 
$(LIBRARY_PATH)(eopt.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh 
$(LIBRARY_PATH)(stpr_gh_at.o):	$(INCDIR)/tcgmsg.fh ./stpr_optimize.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgtypesf.h $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cstprfiles.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(dft_optimize.o):	./stpr_optimize.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cstprfiles.fh 
$(LIBRARY_PATH)(stpr_formup.o):	./cstprfiles.fh 
$(LIBRARY_PATH)(stpr_save.o):	./cstprfiles.fh 
$(LIBRARY_PATH)(stpr_gah_at.o):	$(INCDIR)/tcgmsg.fh ./stpr_optimize.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgtypesf.h $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cstprfiles.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(stpr_walk.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh ./stpr_optimize.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/msgtypesf.h $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geomP.fh $(INCDIR)/geom.fh ./cstprfiles.fh $(INCDIR)/context.fh 
$(LIBRARY_PATH)(stpr_driver.o):	./chow.fh 
$(LIBRARY_PATH)(stpr_wlkchk.o):	$(INCDIR)/util.fh 
$(LIBRARY_PATH)(stpr_sumstc.o):	$(INCDIR)/util.fh ./chow.fh 
$(LIBRARY_PATH)(stpr_partit.o):	$(INCDIR)/util.fh ./chow.fh 
$(LIBRARY_PATH)(stpr_output.o):	$(INCDIR)/util.fh 
$(LIBRARY_PATH)(stpr_mkstep.o):	$(INCDIR)/util.fh ./stpr_clean.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh ./chow.fh 
$(LIBRARY_PATH)(stpr_inihes.o):	$(INCDIR)/util.fh ./stpr_clean.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh ./chow.fh 
$(LIBRARY_PATH)(stpr_hesupd.o):	./chow.fh 
$(LIBRARY_PATH)(stpr_face.o):	$(INCDIR)/util.fh ./stpr_inP.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh ./chow.fh 
$(LIBRARY_PATH)(stpr_cheker.o):	./chow.fh 
