#$Id: GNUmakefile,v 1.37 2001-10-09 23:51:47 windus Exp $

include ../config/makefile.h

       LIBRARY = liblapack.a
       PRECISION = double

#
# Determine if single or double precision is required
#

       VPATH = double:single

ifeq ($(TARGET),CRAY-T3D)
       PRECISION = single
endif
ifeq ($(TARGET),CRAY-T3E)
       PRECISION = single
endif
ifeq ($(TARGET),KSR)
       PRECISION = single
endif

ifndef PRECISION
error:
	@echo Lapack makefile does not have PRECISION defined
	@exit 1
endif

ifeq ($(PRECISION),double)
OBJ_DLAMCH = dlamch.o
OBJ_OPTIMIZE = dbdsqr.o      dgebd2.o      dgebrd.o       dgelq2.o     \
               dgelqf.o      dgeqr2.o      dgeqrf.o       dgesv.o      \
               dgesvd.o      dgetf2.o      dgetrf.o       dgetrs.o     \
               dlabrd.o      dlacon.o      dlacpy.o       dladiv.o     \
               dlae2.o       dlaev2.o      dlagtf.o       dlagts.o     \
                             dlange.o      dlansp.o       dtrti2.o     \
               dlanst.o      dlansy.o      dlapy2.o       dlapy3.o     \
               dlarf.o       dlarfb.o      dlarfg.o       dlarft.o     \
               dlarnv.o      dlartg.o      dlaruv.o       dlas2.o      \
               dlascl.o      dlaset.o      dlasr.o        dlasrt.o     \
               dlassq.o      dlasv2.o      dlaswp.o       dlatrd.o     \
               dlauu2.o      dlauum.o      dorg2l.o       dorg2r.o     \
               dorgbr.o      dorgl2.o      dorglq.o       dorgql.o     \
               dorgqr.o      dorgtr.o      dorm2r.o       dormbr.o     \
               dorml2.o      dormlq.o      dormqr.o       dpotf2.o     \
               dpotrf.o      dpotri.o      dspcon.o       dsprfs.o     \
               dspsvx.o      dsptrf.o      dsptrs.o       dsteqr.o     \
               dsterf.o      dsyev.o       dsygs2.o       dsygst.o     \
               dsygv.o       dsytd2.o      dsytrd.o       dtrtri.o     \
               dlasq1.o      dlasq2.o      dlasq3.o       dlasq5.o     \
               dlasq4.o      dlasq6.o      dposv.o \
               ilaenv.o      lsamen.o                                  \
               zheev.o       zhegs2.o      zhegst.o       zhegv.o      \
               zhetd2.o      zhetrd.o      zlacgv.o       zladiv.o     \
               zlanhe.o      zlarf.o       zlarfb.o       zlarfg.o     \
               zlarft.o      zlascl.o      zlaset.o       zlasr.o      \
               zlassq.o      zlatrd.o      zpotf2.o       zpotrf.o     \
               zsteqr.o      zung2l.o      zung2r.o       zungql.o     \
               zspcon.o      zspmv.o       zspr.o         zsprfs.o     \
               zspsvx.o      zsptrf.o      zsptrs.o       zgeru.o      \
               zlacon.o      zlacpy.o      zlansp.o \
               zungqr.o      zungtr.o

ifeq ($(FC),pgf77)
#pgf77 on Linux (and possibly others) puts dlamch routines into infinite loop
        OBJ = $(OBJ_DLAMCH)
else
	OBJ_OPTIMIZE += $(OBJ_DLAMCH)
endif
else
########################################################################
OBJ_OPTIMIZE = sbdsqr.o      sgebd2.o      sgebrd.o       sgelq2.o     \
               sgelqf.o      sgeqr2.o      sgeqrf.o       sgesv.o      \
               sgesvd.o      sgetf2.o      sgetrf.o       sgetrs.o     \
               slabrd.o      slacon.o      slacpy.o       sladiv.o     \
               slae2.o       slaev2.o      slagtf.o       slagts.o     \
               slamch.o      slange.o      slansp.o       strti2.o     \
               slanst.o      slansy.o      slapy2.o       slapy3.o     \
               slarf.o       slarfb.o      slarfg.o       slarft.o     \
               slarnv.o      slartg.o      slaruv.o       slas2.o      \
               slascl.o      slaset.o      slasr.o        slasrt.o     \
               slassq.o      slasv2.o      slaswp.o       slatrd.o     \
               slauu2.o      slauum.o      sorg2l.o       sorg2r.o     \
               sorgbr.o      sorgl2.o      sorglq.o       sorgql.o     \
               sorgqr.o      sorgtr.o      sorm2r.o       sormbr.o     \
               sorml2.o      sormlq.o      sormqr.o       spotf2.o     \
               spotrf.o      spotri.o      sspcon.o       ssprfs.o     \
               sspsvx.o      ssptrf.o      ssptrs.o       ssteqr.o     \
               ssterf.o      ssyev.o       ssygs2.o       ssygst.o     \
               ssygv.o       ssytd2.o      ssytrd.o       strtri.o     \
               slasq1.o      slasq2.o      slasq3.o       slasq5.o     \
               slasq4.o      sposv.o \
               ilaenv.o      lsamen.o                                  \
               cheev.o       chegs2.o      chegst.o       chegv.o      \
               chetd2.o      chetrd.o      clacgv.o       cladiv.o     \
               clanhe.o      clarf.o       clarfb.o       clarfg.o     \
               clarft.o      clascl.o      claset.o       clasr.o      \
               classq.o      clatrd.o      cpotf2.o       cpotrf.o     \
               csteqr.o      cung2l.o      cung2r.o       cungql.o     \
               cungqr.o      cungtr.o
endif

ifeq ($(TARGET),CRAY-T3D)
OBJ_OPTIMIZE = slagtf.o slagts.o lsamen.o slaset.o \
               ssyev.o slansy.o ssteqr.o ssterf.o slae2.o \
               slaev2.o slasr.o slanst.o slascl.o slasrt.o \
               ssygv.o ssygst.o ssygs2.o slas2.o
endif

ifeq ($(TARGET),CRAY-T3E)
OBJ_OPTIMIZE = slagtf.o slagts.o lsamen.o slaset.o \
               ssyev.o slansy.o ssteqr.o ssterf.o slae2.o \
               slaev2.o slasr.o slanst.o slascl.o slasrt.o \
               ssygv.o ssygst.o ssygs2.o slas2.o clanhe.o clascl.o
endif



include ../config/makelib.h



