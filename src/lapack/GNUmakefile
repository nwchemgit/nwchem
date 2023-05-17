#$Id$
ifeq ($(BLAS_SIZE),4)
USE_I4FLAGS=y
endif
include ../config/makefile.h

       LIBRARY = libnwclapack.a
       PRECISION = double

#
# Determine if single or double precision is required
#

       VPATH = double:single

ifeq ($(TARGET),KSR)
       PRECISION = single
endif

ifeq ($(USE_PRECISION),single)
       PRECISION = single
endif
ifeq ($(USE_PRECISION),double)
       PRECISION = double
endif

ifndef PRECISION
error:
	@echo Lapack makefile does not have PRECISION defined
	@exit 1
endif

ifeq ($(PRECISION),double)
OBJ_DLAMCH = dlamch.o dlaisnan.o 
OBJ_OPTIMIZE = dbdsqr.o      dgebd2.o      dgebrd.o       dgelq2.o     \
               dgelqf.o      dgeqr2.o      dgeqrf.o       dgesv.o      \
               dgesvd.o      dgetf2.o      dgetrf.o       dgetrs.o     \
               disnan.o dlahr2.o dlaqr0.o \
               dlabrd.o      dlacon.o      dlacpy.o       dladiv.o     \
               dlae2.o       dlaev2.o  dlaexc.o    dlagtf.o       dlagts.o     \
		dgels.o       dlange.o      dlansp.o       dtrti2.o     \
               dlanst.o      dlansy.o      dlapy2.o       dlapy3.o     \
               dlarf.o       dlarfb.o      dlarfg.o       dlarft.o     \
               dlarnv.o      dlartg.o      dlaruv.o       dlas2.o  dlasy2.o    \
               dlacn2.o dlaqr1.o dlaqr2.o dlaqr3.o dlaqr4.o dlaqr5.o  \
               dlascl.o      dlaset.o      dlasr.o        dlasrt.o     \
               dlassq.o      dlasv2.o      dlaswp.o       dlatrd.o     \
               dlauu2.o      dlauum.o      dorg2l.o       dorg2r.o     \
               dormhr.o      dgtsv.o \
               dorgbr.o      dorgl2.o      dorglq.o       dorgql.o     \
               dorgqr.o      dorgtr.o      dorm2r.o       dormbr.o     \
               dorml2.o      dormlq.o      dormqr.o       dpotf2.o     \
               dpotrf.o      dpotri.o      dspcon.o       dsprfs.o     \
               dspsvx.o      dsptrf.o      dsptrs.o       dsteqr.o     \
               dsterf.o      dsyev.o       dsygs2.o       dsygst.o     dsysv.o \
               dsygv.o       dsytd2.o      dsytrd.o       dtrexc.o dtrtri.o     \
               dlasq1.o      dlasq2.o      dlasq3.o       dlasq5.o     \
               dlasq4.o      dlasq6.o      dposv.o        dtrtrs.o     \
               dpotrs.o      dlaebz.o      dstedc.o       \
               dlaed0.o dlaed1.o dlaed2.o dlaed3.o dlaed4.o dlaed5.o dlaed6.o dlaed7.o dlaed8.o dlaed9.o \
               dlamrg.o      dlaeda.o dsyevd.o dsyevx.o  dstebz.o dstein.o \
               ilaenv.o ieeeck.o      lsamen.o  iladlr.o iladlc.o ilazlc.o ilazlr.o iparmq.o \
	       zbdsqr.o  zgebd2.o zgebrd.o zgelq2.o zgelqf.o zgeqr2.o \
               zgeqrf.o zgesvd.o \
               zheev.o       zhegs2.o      zhegst.o       zhegv.o      \
               zhetd2.o      zhetrd.o      zlacgv.o       zladiv.o     \
               zlanhe.o      zlarf.o       zlarfb.o       zlarfg.o     \
               zlabrd.o \
               zlarft.o      zlascl.o      zlaset.o       zlasr.o      \
               zlassq.o      zlatrd.o      zpotf2.o       zpotrf.o     \
               zsteqr.o      zung2l.o      zung2r.o       zungql.o     \
               zspcon.o      zspmv.o       zspr.o         zsprfs.o     \
               zspsvx.o      zsptrf.o      zsptrs.o            \
               zlacon.o      zlacpy.o      zlansp.o  zungbr.o zungl2.o zunglq.o \
               zungqr.o      zungtr.o zunm2r.o zunmbr.o zunml2.o zunmlq.o \
               zunmqr.o \
               zlarfx.o zunghr.o zlatrs.o zgebak.o zgebal.o zgeev.o \
               zgehd2.o zgehrd.o zgesv.o zgetf2.o zgetrf.o zgetrs.o \
               zhseqr.o zlahqr.o zlahrd.o zlange.o zlanhs.o dlabad.o \
               ztrevc.o zlaswp.o zrot.o zlaqr3.o ztrexc.o zunmhr.o zlahr2.o zlaqr0.o\
               zlaqr4.o zlartg.o zlaqr1.o xerbla.o zlaqr5.o zlaqr2.o lsame.o \
               dgeev.o dgeevx.o dgebak.o dgebal.o dgehrd.o dtrsna.o dlaqtr.o \
               dhseqr.o dorghr.o dtrevc.o dlaln2.o dlahqr.o dgehd2.o \
               dlanhs.o dlarfx.o dlahrd.o dlanv2.o \
               dlarra.o dlarrb.o dlarrc.o dlarrd.o dlarrf.o dlarrk.o dlarrv.o\
               dlar1v.o dlaneg.o dgelss.o dormtr.o dormql.o dorm2l.o \
               dsfrk.o dsytrf.o  dgetri.o dtrsm.o dtrmv.o dtrmm.o\
	       dpftrf.o dpftri.o dtfsm.o dtftri.o zlacp2.o \
	       zsysv.o zsytrf.o zsytf2.o zsytrs.o zsytrs2.o zsyr.o \
	       zsyconv.o zlasyf.o

#pgf77 on Linux, hpux (and possibly others) puts dlamch routines into infinite loop
        OBJ = $(OBJ_DLAMCH)
else
########################################################################
OBJ_OPTIMIZE = sbdsqr.o      sgebd2.o      sgebrd.o       sgelq2.o     \
               sgelqf.o      sgeqr2.o      sgeqrf.o       sgesv.o      \
               sgesvd.o      sgetf2.o      sgetrf.o       sgetrs.o     \
               sisnan.o      slahr2.o      slaqr0.o       slaqr1.o     \
               slaqr2.o      slaqr3.o      slaqr4.o       slaqr5.o     \
               slasy2.o      sgehd2.o \
               slabrd.o      slacon.o      slacpy.o       sladiv.o     \
               slae2.o       slaev2.o      slagtf.o       slagts.o     \
               slamch.o      slange.o      slansp.o       strti2.o     \
               slanv2.o      slahqr.o      slacn2.o       slarfx.o     \
               slanst.o      slansy.o      slapy2.o       slapy3.o     \
               slarf.o       slarfb.o      slarfg.o       slarft.o     \
               slarnv.o      slartg.o      slaruv.o       slas2.o      \
               slascl.o      slaset.o      slasr.o        slasrt.o     \
               slassq.o      slasv2.o      slaswp.o       slatrd.o     \
               slauu2.o      slauum.o      sorg2l.o       sorg2r.o     \
               sorgbr.o      sorgl2.o      sorglq.o       sorgql.o     \
               sorgqr.o      sorgtr.o      sorm2r.o       sormbr.o     \
               sormhr.o      sgeev.o       slaisnan.o     slaln2.o     \
               sorml2.o      sormlq.o      sormqr.o       spotf2.o     \
               spotrf.o      spotri.o      sspcon.o       ssprfs.o     \
               sspsvx.o      ssptrf.o      ssptrs.o       ssteqr.o     \
               ssterf.o      ssyev.o       ssygs2.o       ssygst.o     \
               ssygv.o       ssytd2.o      ssytrd.o       strtri.o     \
               slasq1.o      slasq2.o      slasq3.o       slasq5.o     \
               slasq4.o      slasq6.o      sposv.o        strtrs.o\
               spotrs.o      slaebz.o      sstedc.o       \
               slaed0.o slaed1.o slaed2.o slaed3.o slaed4.o slaed5.o slaed6.o slaed7.o slaed8.o slaed9.o \
               slamrg.o      slaeda.o      strexc.o       slabrd.o     \
               slaexc.o      ilaslr.o      ilaslc.o       ilaclc.o     \
               ilaclr.o      iparmq.o \
               ilaenv.o    ieeeck.o   lsamen.o                                  \
               cheev.o       chegs2.o      chegst.o       chegv.o      \
               chetd2.o      chetrd.o      clacgv.o       cladiv.o     \
               clanhe.o      clarf.o       clarfb.o       clarfg.o     \
               clarft.o      clascl.o      claset.o       clasr.o      \
               classq.o      clatrd.o      cpotf2.o       cpotrf.o     \
               csteqr.o      cung2l.o      cung2r.o       cungql.o     \
               cungqr.o      cungtr.o      cungbr.o       cunm2r.o     \
               cunmbr.o      cunml2.o      cunmlq.o       cungl2.o     \
               cgehd2.o cgehrd.o clahrd.o cgebak.o cgebal.o cgeev.o \
               chseqr.o clacpy.o clahqr.o clange.o clanhs.o clarfx.o \
               clatrs.o ctrevc.o cunghr.o slabad.o claswp.o \
               sgebak.o sgebal.o sgehrd.o shseqr.o sorghr.o strevc.o \
               ssfrk.o ssytrf.o  \
	       spftrf.o spftri.o stfsm.o stftri.o clacp2.o

endif

include ../config/makelib.h

ifeq ($(_FC),ifc)
  ifeq ($(_IFCV81),Y)
    FOPTIMIZE +=-ip_no_inlining
  endif
endif
#this fixes FPE in lapack/dlagts for dft_feco5
#this fixes hangs in dlamch
ifeq ($(_FC),g77)
    FDEBUG +=-ffloat-store -fno-fast-math
endif
ifeq ($(_FC),gfortran)
    FDEBUG +=-ffloat-store -fno-fast-math
endif

