

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

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(upkcon.o):	$(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(tester.o):	./mptr.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(stool.o):	$(INCDIR)/global.fh ./cselcifiles.fh 
$(LIBRARY_PATH)(slect.o):	./mptr.fh ./implicit.fh $(INCDIR)/global.fh ./ciinfo.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(sigma.o):	$(INCDIR)/global.fh ./cbuff.fh 
$(LIBRARY_PATH)(select.o):	./mptr.fh ./implicit.fh $(INCDIR)/global.fh ./cselcifiles.fh ./ciinfo.fh 
$(LIBRARY_PATH)(selci.o):	$(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cselcifiles.fh 
$(LIBRARY_PATH)(restk.o):	./cselcifiles.fh 
$(LIBRARY_PATH)(restin.o):	$(INCDIR)/global.fh ./cselcifiles.fh 
$(LIBRARY_PATH)(readh.o):	./implicit.fh 
$(LIBRARY_PATH)(rdhint.o):	./implicit.fh 
$(LIBRARY_PATH)(rdhead.o):	./implicit.fh ./cselcifiles.fh ./ciinfo.fh 
$(LIBRARY_PATH)(putinb.o):	./implicit.fh ./cbuff.fh 
$(LIBRARY_PATH)(ptidy.o):	./implicit.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(prtri.o):	./implicit.fh 
$(LIBRARY_PATH)(popcnt32.o):	./cbits.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(pkcon.o):	$(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(pconls.o):	./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(ointer.parthack.o):	$(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(ointer.mod.o):	$(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(ointer.o):	$(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(mpdele.o):	./implicit.fh 
$(LIBRARY_PATH)(mkvc.o):	./mptr.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(mknsneed.o):	./implicit.fh 
$(LIBRARY_PATH)(mkinta.o):	./implicit.fh ./ciinfo.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(mkindx.o):	./implicit.fh ./ciinfo.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(mkdele.o):	./implicit.fh 
$(LIBRARY_PATH)(makhdb.o):	./stfunc.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(makehs.o):	./stfunc.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(makehd.o):	./stfunc.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(makeh.o):	./implicit.fh $(INCDIR)/global.fh ./ciinfo.fh ./cbuff.fh 
$(LIBRARY_PATH)(makef.o):	./stfunc.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(loadw.o):	./implicit.fh 
$(LIBRARY_PATH)(loadg.o):	./implicit.fh ./ciinfo.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(loadab.o):	./implicit.fh ./ciinfo.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(loada.o):	./implicit.fh ./ciinfo.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(load.o):	./implicit.fh ./ciinfo.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(iprtri.o):	./implicit.fh 
$(LIBRARY_PATH)(iodiff.o):	$(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(input.o):	./workc.fh ./work.fh 
$(LIBRARY_PATH)(initci.o):	./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(inicij.o):	./implicit.fh ./ciinfo.fh ./ceij.fh 
$(LIBRARY_PATH)(gtijkl.o):	$(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(getij.o):	$(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(fcidim.o):	./implicit.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(eirerj.o):	./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(eijkl.o):	./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(eijkj.o):	./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(eijil.o):	./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(eijeji.o):	./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(eij.o):	./implicit.fh ./ciinfo.fh ./ceij.fh 
$(LIBRARY_PATH)(e0mp2.o):	./implicit.fh ./cmp2.fh ./ciinfo.fh 
$(LIBRARY_PATH)(davids.o):	./mptr.fh ./implicit.fh $(INCDIR)/global.fh ./cselcifiles.fh ./cbuff.fh 
$(LIBRARY_PATH)(conf.o):	$(INCDIR)/rtdb.fh ./mptr.fh $(INCDIR)/mafdecls.fh ./implicit.fh $(INCDIR)/global.fh ./cselcifiles.fh ./cconf.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(ci.o):	./mptr.fh ./implicit.fh $(INCDIR)/global.fh ./cselcifiles.fh ./ciinfo.fh 
$(LIBRARY_PATH)(anal.o):	./implicit.fh ./ciinfo.fh 
