#$Id: GNUmakefile,v 1.24 1996-07-13 00:21:06 gutowski Exp $

include ../config/makefile.h

       LIBRARY = liblapack.a

#
# Determine if single or double precision is required
#

       VPATH = double:single

ifeq ($(TARGET),CRAY-T3D)
       PRECISION = single
endif
ifeq ($(TARGET),SUN)
       PRECISION = double
endif
ifeq ($(TARGET),SOLARIS)
       PRECISION = double
endif
ifeq ($(TARGET),SGITFP)
       PRECISION = double
endif
ifeq ($(TARGET),SGI)
       PRECISION = double
endif
ifeq ($(TARGET),IBM)
       PRECISION = double
endif
ifeq ($(TARGET),SP1)
       PRECISION = double
endif
ifeq ($(TARGET),KSR)
       PRECISION = double
endif
ifeq ($(TARGET),CRAY-T3D)
       PRECISION = single
endif
ifeq ($(TARGET),PARAGON)
       PRECISION = double
endif
ifeq ($(TARGET),DELTA)
       PRECISION = double
endif
ifeq ($(TARGET),IPSC)
       PRECISION = double
endif
ifeq ($(TARGET),DECOSF)
       PRECISION = double
endif
ifeq ($(TARGET),LINUX)
	PRECISION = double
endif

ifndef PRECISION
error:
	@echo Lapack makefile does not have PRECISION defined
	@exit 1
endif

ifeq ($(PRECISION),double)
OBJ_OPTIMIZE = dlabrd.o   dlarfb.o   dlasv2.o   dormlq.o \
               dbdsqr.o   dlacon.o   dlarfg.o   dorg2r.o   dormqr.o \
               dgebd2.o   dlacpy.o   dlarft.o   dorgbr.o   dspcon.o \
               dgebrd.o   dlaev2.o   dlartg.o   dorgl2.o   dsprfs.o \
               dgelq2.o   dlamch.o   dlas2.o    dorglq.o   dspsvx.o \
               dgelqf.o   dlange.o   dlascl.o   dorgqr.o   dsptrf.o \
               dgeqr2.o   dlansp.o   dlaset.o   dorm2r.o   dsptrs.o \
               dgeqrf.o   dlapy2.o   dlasr.o    dormbr.o   ilaenv.o \
               dgesvd.o   dlarf.o    dlassq.o   dorml2.o   lsamen.o \
               dgesv.o    dgetrs.o   dgetrf.o   dlaswp.o \
               dgetf2.o \
               dlagtf.o   dlagts.o   dlarnv.o   dlaruv.o   dpotrf.o \
               dpotri.o   dtrtri.o   dlauum.o   dlauu2.o   dpotf2.o \
               dsyev.o    dsygst.o   dsygv.o    dsytrd.o   dorgtr.o \
               dsteqr.o   dlansy.o   dsygs2.o   dsterf.o   dsytd2.o \
               dorgql.o   dlatrd.o   dlae2.o    dlazro.o   dorg2l.o \
               dladiv.o   dlae2.o    dlamc1.o   dlamc2.o   dlamc3.o \
               dlamc4.o   dlamc5.o   dlanst.o   dlapy3.o   dlasrt.o \
               dsterf.o   zheev.o    zhegs2.o   zhegst.o   zhegv.o \
               zhetd2.o   zhetrd.o   zlacgv.o   zladiv.o   zlanhe.o \
               zlarf.o    zlarfb.o   zlarfg.o   zlarft.o   zlascl.o \
               zlaset.o   zlasr.o    zlassq.o   zlatrd.o   zpotf2.o \
               zpotrf.o   zsteqr.o   zung2l.o   zung2r.o   zungql.o \
               zungqr.o   zungtr.o
else
OBJ_OPTIMIZE = sgesvd.o   slarf.o    slassq.o   sorml2.o   lsamen.o \
               ilaenv.o   slabrd.o   slarfb.o   slasv2.o   sormlq.o \
               sbdsqr.o   slacon.o   slarfg.o   sorg2r.o   sormqr.o \
               sgebd2.o   slacpy.o   slarft.o   sorgbr.o   sspcon.o \
               sgebrd.o   slaev2.o   slartg.o   sorgl2.o   ssprfs.o \
               sgelq2.o   slamch.o   slas2.o    sorglq.o   sspsvx.o \
               sgelqf.o   slange.o   slascl.o   sorgqr.o   ssptrf.o \
               sgeqr2.o   slansp.o   slaset.o   sorm2r.o   ssptrs.o \
               sgeqrf.o   slapy2.o   slasr.o    sormbr.o \
               sgesv.o    sgetrs.o   sgetrf.o   slaswp.o \
               sgetf2.o \
               slagtf.o   slagts.o   slarnv.o   slaruv.o   spotrf.o \
               spotri.o   strtri.o   slauum.o   slauu2.o   spotf2.o \
               ssyev.o    ssygst.o   ssygv.o    ssytrd.o   sorgtr.o \
               ssteqr.o   slansy.o   ssygs2.o   ssterf.o   ssytd2.o \
               sorgql.o   slatrd.o   slae2.o    slazro.o   sorg2l.o \

endif

ifeq ($(TARGET),CRAY-T3D)
OBJ_OPTIMIZE = slagtf.o slagts.o lsamen.o slaset.o \
               ssyev.o slansy.o ssteqr.o ssterf.o slae2.o \
               slaev2.o slasr.o slazro.o \
               ssygv.o ssygst.o ssygs2.o
endif


include ../config/makelib.h



# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
