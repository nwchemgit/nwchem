

LIBRARY = libselci.a
LIB_TARGETS = conf htoa

BLAS     = mxma.o

UOBJ = \
iprtri.o	prtri.o  \
input.o		icopy.o		output.o	yacobi.o \
tred1.o		tqlrat.o	loadw.o		rdhint.o 	rdhcon.o \
rdhwmt.o	rdconf.o	upkcon.o	prtcon.o	wrtcon.o \
iodiff.o	bserch.o	popcnt32.o	readh.o  \
renorm.o	sread.o		swrite.o	pkcon.o  \
wthcon.o	wtconf.o	vadd.o		idmin.o	\
isum.o		rd1.o		fcidim.o        \
axb.o		axbt.o		gettim.o	cputm.o \
minit.o		mxmb.o  check.o

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
select.o        selci.o         ptidy.o         anal.o          mkindx.o \
slect.o         tester.o        mkdele.o        mkvc.o          e0mp2.o \
mpdele.o	stool.o 	restk.o		restin.o        orbanal.o

OBJ_OPTIMIZE = $(UOBJ) $(DVOBJ) $(CIOBJ) $(SELOBJ) $(BLAS)

USES_BLAS = axbt.f davids.F e0mp2.F fcidim.F makeh.F makehs.F makhdb.F mkvc.F tester.F

include ../config/makefile.h
include ../config/makelib.h

 LIBS = $(LIBRARY_PATH) $(LIBDIR)/libtcgmsg.a

conf:	conf.o $(LIBRARY_PATH)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

htoa:	htoa.o
	$(FC) -o $@ $^

