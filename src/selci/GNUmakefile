
# $Id$
#


LIBRARY = libselci.a
LIB_TARGETS = conf htoa

BLAS     = mxma.o

UOBJ = \
iprtri.o	prtri.o  \
	selci_icopy.o		selci_yacobi.o \
tred1.o		tqlrat.o	loadw.o		rdhint.o 	rdhcon.o \
rdhwmt.o	rdconf.o	upkcon.o	prtcon.o	wrtcon.o \
iodiff.o	bserch.o	readh.o  \
selci_renorm.o	selci_sread.o		selci_swrite.o	pkcon.o  \
wthcon.o	wtconf.o	vadd.o		idmin.o	\
isum.o		rd1.o		fcidim.o        \
axb.o		axbt.o		gettim.o	cputm.o \
minit.o		mxmb.o  selci_check.o

DVOBJ = davids.o sigma.o 

CIOBJ = \
conf.o \
ci.o            initci.o	rdhead.o 	load.o   	loada.o  \
mkinta.o	pconls.o	makeh.o		mkindb.o	eij.o    \
eijeji.o	makehd.o	makef.o 	loadg.o \
getij.o		makehs.o	eirerj.o 	loadab.o        couple.o \
gtijkl.o	eijkj.o		makhdb.o	eijil.o         mknsneed.o \
eijkl.o		putinb.o	inicij.o	ointer.o	oeq.o 

SELOBJ = \
selci_select.o        selci.o         ptidy.o         anal.o          mkindx.o \
slect.o         tester.o        mkdele.o        mkvc.o          e0mp2.o \
mpdele.o	stool.o 	restk.o		restin.o        orbanal.o \
selci_hpp.o

OBJ_OPTIMIZE = $(UOBJ) $(DVOBJ) $(CIOBJ) $(SELOBJ) $(BLAS)
#OBJ = $(UOBJ) $(DVOBJ) $(CIOBJ) $(SELOBJ) $(BLAS)

USES_BLAS = davids.F fcidim.F makeh.F makehs.F makhdb.F tester.F selci_hpp.F

include ../config/makefile.h
include ../config/makelib.h


ifeq ($(_IFCV81),Y)
   FOPTIMIZE=-O3
   COPTIONS=-O0 -g
endif
ifeq ($(TARGET),LINUX64)
  ifeq ($(_CPU),ia64)
   COPTIONS=-O0 -g
  endif
endif

 LIBS = $(LIBRARY_PATH) $(LIBDIR)/libtcgmsg.a

conf:	conf.o $(LIBRARY_PATH)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

htoa:	htoa.o
	$(FC) -o $@ $^

