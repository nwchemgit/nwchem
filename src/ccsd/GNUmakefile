	  OBJ = aoccsd.o \
		ccsd.o \
		ccsd_aux.o \
		ccsd_driver.o \
		ccsd_directbf.o \
		ccsd_fsig1.o \
		ccsd_fsig1ao.o \
		ccsd_fsig2.o \
		ccsd_get_param.o \
		ccsd_intmeda.o \
		ccsd_intmedb.o \
		ccsd_iterdrv.o \
		ccsd_mkcinv.o \
		ccsd_mktask.o \
		ccsd_mktau.o \
		ccsd_mkt2.o \
		ccsd_pampt.o \
		ccsd_pdiis.o \
		ccsd_sortmo.o \
		ccsd_rdtrpo.o \
		ccsd_tengy.o \
		ccsd_trpdrv.o \
		moints_trfidx.o \
		moints_trp.o \
		qstat.o \
		cputime.o \
		rddcom.o

    LIBRARY = libccsd.a

  USES_BLAS = \
	      ccsd_directbf.F \
	      ccsd_fsig1.f \
	      ccsd_fsig1ao.F \
	      ccsd_fsig2.f \
	      ccsd_intmeda.f \
	      ccsd_intmedb.f \
	      ccsd_iterdrv.f \
	      ccsd_pampt.F \
	      ccsd_pdiis.f \
	      ccsd_sortmo.F \
	      ccsd_trpdrv.f \
	      moints_trp.F



LIB_DEFINES = 

#
#


include ../config/makefile.h
include ../config/makelib.h
# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
aoccsd.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
ccsd.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
ccsd_directbf.o :	 $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/bas.fh 
ccsd_driver.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
ccsd_get_param.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
ccsd_mkcinv.o :	 $(INCDIR)/bas.fh 
ccsd_mktask.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/bas.fh 
ccsd_pampt.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/bas.fh 
moints_trp.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/schwarz.fh $(INCDIR)/msgids.fh $(INCDIR)/moints_stats.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
rddcom.o :	 $(INCDIR)/global.fh 
