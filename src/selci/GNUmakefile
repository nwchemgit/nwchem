

LIBRARY = libselci.a
LIB_TARGETS = ci select davids conf

BLAS     = mxma.o ddot.o dnrm2.o dcopy.o dscal.o daxpy.o

UOBJ = \
error.o		ifill.o		iprtri.o	prtri.o  \
input.o		icopy.o		output.o	yacobi.o \
tred1.o		tqlrat.o	loadw.o		rdhint.o 	rdhcon.o \
rdhwmt.o	rdconf.o	upkcon.o	prtcon.o	wrtcon.o \
iodiff.o	bserch.o	popcnt32.o	readh.o  \
renorm.o	sread.o		swrite.o	getnam.o	pkcon.o  \
wthcon.o	wtconf.o	vadd.o		idmin.o	\
isum.o		rd1.o		fcidim.o        pfname.o \
axb.o		axbt.o		gettim.o	cputm.o \
minit.o		ieeetrap.o \
idamax.o	\
dfill.o         dasum.o         mxmb.o 

DVOBJ = sigma.o sigsav.o

CIOBJ = \
initci.o	rdhead.o 	load.o   	loada.o  \
mkinta.o	pconls.o	makeh.o		mkindb.o	eij.o    \
eijeji.o	makehd.o	makef.o 	loadg.o \
getij.o		makehs.o	eirerj.o 	loadab.o        couple.o \
gtijkl.o	eijkj.o		makhdb.o	eijil.o         mknsneed.o \
eijkl.o		putinb.o	inicij.o	ointer.o	oeq.o 

SELOBJ = \
mselct.o        ptidy.o         anal.o          mkindx.o \
slect.o         tester.o        mkdele.o        mkvc.o          e0mp2.o \
mpdele.o	stool.o 	restk.o		restin.o

OBJ_OPTIMIZE = $(UOBJ) $(DVOBJ) $(CIOBJ) $(SELOBJ) $(BLAS)

include ../config/makefile.h
include ../config/makelib.h

 LIBS = $(LIBRARY_PATH) $(LIBDIR)/libtcgmsg.a

all:	select conf ci davids
	@echo Finished making all.

select:	select.o $(LIBRARY_PATH)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

davids:	davids.o $(LIBRARY_PATH)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

ci:	ci.o $(LIBRARY_PATH)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

conf:	conf.o $(LIBRARY_PATH)
	$(FC) $(LDFLAGS) -o $@ $^ $(LIBS)

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(anal.o):	 ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(ci.o):	 ./mptr.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(conf.o):	 ./mptr.fh ./implicit.fh ./cconf.fh ./bitops.fh 
$(LIBRARY_PATH)(davids.o):	 ./mptr.fh ./implicit.fh ./cbuff.fh 
$(LIBRARY_PATH)(e0mp2.o):	 ./implicit.fh ./cmp2.fh ./ciinfo.fh 
$(LIBRARY_PATH)(eij.o):	 ./implicit.fh ./ciinfo.fh ./ceij.fh 
$(LIBRARY_PATH)(eijeji.o):	 ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(eijil.o):	 ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(eijkj.o):	 ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(eijkl.o):	 ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(eirerj.o):	 ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(fcidim.o):	 ./implicit.fh ./bitops.fh 
$(LIBRARY_PATH)(getij.o):	 ./bitops.fh 
$(LIBRARY_PATH)(gtijkl.o):	 ./bitops.fh 
$(LIBRARY_PATH)(inicij.o):	 ./implicit.fh ./ciinfo.fh ./ceij.fh 
$(LIBRARY_PATH)(initci.o):	 ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(input.o):	 ./workc.fh ./work.fh 
$(LIBRARY_PATH)(iodiff.o):	 ./bitops.fh 
$(LIBRARY_PATH)(iprtri.o):	 ./implicit.fh 
$(LIBRARY_PATH)(load.o):	 ./implicit.fh ./ciinfo.fh ./bitops.fh 
$(LIBRARY_PATH)(loada.o):	 ./implicit.fh ./ciinfo.fh ./bitops.fh 
$(LIBRARY_PATH)(loadab.o):	 ./implicit.fh ./ciinfo.fh ./bitops.fh 
$(LIBRARY_PATH)(loadg.o):	 ./implicit.fh ./ciinfo.fh ./bitops.fh 
$(LIBRARY_PATH)(loadw.o):	 ./implicit.fh 
$(LIBRARY_PATH)(makef.o):	 ./stfunc.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(makeh.o):	 ./implicit.fh ./ciinfo.fh ./cbuff.fh 
$(LIBRARY_PATH)(makehd.o):	 ./stfunc.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(makehs.o):	 ./stfunc.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(makhdb.o):	 ./stfunc.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(mkdele.o):	 ./implicit.fh 
$(LIBRARY_PATH)(mkindx.o):	 ./implicit.fh ./ciinfo.fh ./bitops.fh 
$(LIBRARY_PATH)(mkinta.o):	 ./implicit.fh ./ciinfo.fh ./bitops.fh 
$(LIBRARY_PATH)(mknsneed.o):	 ./implicit.fh 
$(LIBRARY_PATH)(mkvc.o):	 ./mptr.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(mpdele.o):	 ./implicit.fh 
$(LIBRARY_PATH)(mselct.o):	 ./mptr.fh ./implicit.fh ./cstool.fh ./ciinfo.fh 
$(LIBRARY_PATH)(ointer.o):	 ./bitops.fh 
$(LIBRARY_PATH)(ointer.mod.o):	 ./bitops.fh 
$(LIBRARY_PATH)(pconls.o):	 ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(pkcon.o):	 ./bitops.fh 
$(LIBRARY_PATH)(popcnt32.o):	 ./cbits.fh ./bitops.fh 
$(LIBRARY_PATH)(prtri.o):	 ./implicit.fh 
$(LIBRARY_PATH)(ptidy.o):	 ./implicit.fh 
$(LIBRARY_PATH)(putinb.o):	 ./implicit.fh ./cbuff.fh 
$(LIBRARY_PATH)(rdhead.o):	 ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(rdhint.o):	 ./implicit.fh 
$(LIBRARY_PATH)(readh.o):	 ./implicit.fh 
$(LIBRARY_PATH)(restin.o):	 ./cstool.fh 
$(LIBRARY_PATH)(restk.o):	 ./cstool.fh 
$(LIBRARY_PATH)(select.o):	 ./mptr.fh ./implicit.fh 
$(LIBRARY_PATH)(sigma.o):	 ./cbuff.fh 
$(LIBRARY_PATH)(sigsav.o):	 ./cbuff.fh 
$(LIBRARY_PATH)(slect.o):	 ./mptr.fh ./implicit.fh ./ciinfo.fh ./bitops.fh 
$(LIBRARY_PATH)(stool.o):	 ./cstool.fh 
$(LIBRARY_PATH)(tester.o):	 ./mptr.fh ./implicit.fh ./ciinfo.fh 
$(LIBRARY_PATH)(upkcon.o):	 ./bitops.fh 
