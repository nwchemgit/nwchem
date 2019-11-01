## $Id$

include ../../../config/makefile.h
include ../../DEFS

COBJ  = \
	clustrfix.o \
	mgscs.o \
	pscale.o \
	clustrxx5.o \
	clustrxx4.o \
	pstein4.o \
	pstein5.o \
	peigs_dlasq1.o \
	b_ortho.o \
         chol_bcst.o \
         choleski9.o \
         ci_entry.o \
         conjug22.o \
         de_sym.o \
         exit.o \
         exit2.o \
         forLL.o \
         forLU.o \
         gmax.o \
         inverse.o \
         lu_mxm2.o \
         mapdif.o \
         mapdif1.o \
         matmax.o \
         mdif2b.o \
        mdiff1.o \
         mdiff2.o \
         memreq.o \
         memreq_f.o \
	mgs5.o \
         mxm25.o \
         mxm5x.o \
         mxm88.o \
         mxm_ll1.o \
         onenorm.o \
         ortho.o \
         pdspev_c.o \
         pdspevx.o \
         pdspgv_c.o \
         pdspgvx.o \
         pdsptri.o \
         peigs_cmod.o \
         pipe_fut.o \
         pmmLSL.o \
         pmmLUL.o \
         pmmlsl2.o \
         pxerbla.o \
         qsort.o \
         resid.o \
         residual.o \
         sfnorm.o \
         shellsort.o \
         soluf.o \
         sonenorm.o \
         soort.o \
         sort.o \
         sorteig.o \
         tred22.o \
         treesort3.o \
         tresid.o \
         util.o \
         xerbla.o \
         mapchk.o \
         pdiff.o \
         pgexit.o \
         reducemap.o\
	pstebz10.o \
	pstebz11.o \
	 inv_it4.o \
	 inv_it5.o \
	inv_it3.o \
	peigs_tldlfact.o \
	clustrf5.o \
	r_ritz.o \
	gdot.o 

#	clustrf4.o \
#
#	inv_it3.o
#	 inv_it4.o 
#	inv_it6.o
#      	 pstebz10.o
#	 clustr31.o
#	pdcomplex.o \


all:		$(COBJ)
		$(peigs_AR) $(peigs_LIB) $?

#.EXPORT_ALL_VARIABLES:
.PHONY:	clean
clean:		
	-rm -f *.o maxx.f F*.f
	-rm -f *.c.*

realclean:	
	clean

maxx.o :	maxx.cpp
		-rm -f maxx.f
		$(peigs_CPP) < maxx.cpp > maxx.f
		$(peigs_FC) $(OPT) -c $(F77NFLAGS) maxx.f

soort.o:	soort.c
		$(peigs_CC) $(OPTC) -c soort.c

lower_mxm.o:	lower_mxm.c
		$(peigs_CC) -c lower_mxm.c

shellsort.o:	shellsort.c
		$(peigs_CC) $(OPTC) -c shellsort.c

treesort.o:	treesort.c
		$(peigs_CC) -c treesort.c

gmax.o:		gmax.c
		$(peigs_CC) $(OPTC) -c gmax.c

forLL.o:	forLL.c
		$(peigs_CC) $(OPTC) -c forLL.c

forLU.o:	forLU.c
		$(peigs_CC) $(OPTC) -c forLU.c

ci_entry.o:	ci_entry.c
		$(peigs_CC) $(OPTC) -c ci_entry.c

util.o:		util.c
		$(peigs_CC) $(OPTC) -c util.c

pdspgv_c.o:	pdspgv_c.c
		$(peigs_CC) $(OPTC) -c pdspgv_c.c

pdspev_c.o:	pdspev_c.c
		$(peigs_CC) $(OPTC) -c pdspev_c.c

pdspevx.o:	pdspevx.c
		$(peigs_CC) $(OPTC) -c pdspevx.c

pdspgvx.o:	pdspgvx.c
		$(peigs_CC) $(OPTC) -c pdspgvx.c

exit.o:		exit.c
		$(peigs_CC) $(OPTC) -c exit.c

exit2.o:	exit2.c
		$(peigs_CC) $(OPTC) -c exit2.c

pxerbla.o:	pxerbla.c
		$(peigs_CC) $(OPTC) -c pxerbla.c

chol_bcst.o:	chol_bcst.c
		$(peigs_CC) $(OPTC) -c chol_bcst.c

upperxfull.o:	upperxfull.c
		$(peigs_CC) $(OPTC) -c upperxfull.c

pipe_fut.o:	pipe_fut.c
		$(peigs_CC) $(OPTC) -c pipe_fut.c

choleski9.o:	choleski9.c
		$(peigs_CC) $(OPTC) -c choleski9.c

tred22.o:	tred22.c
		$(peigs_CC) $(OPTC) -c tred22.c

pstebz9.o:	pstebz9.c
		$(peigs_CC) $(OPTC) -c pstebz9.c


pstein4.o:	pstein4.c
		$(peigs_CC) $(OPTC) -c pstein4.c

pstein5.o:	pstein5.c
		$(peigs_CC) $(OPTC) -c pstein5.c

inverse.o:	inverse.c
		$(peigs_CC) $(OPTC)  -c inverse.c

conjug.o:	conjug.c
		$(peigs_CC) $(OPTC) -c conjug.c 

pmmLUL.o:	pmmLUL.c
		$(peigs_CC) $(OPTC) -c pmmLUL.c

pmmLSL.o:	pmmLSL.c
		$(peigs_CC) $(OPTC) -c pmmLSL.c 

soluf.o:	soluf.c
		$(peigs_CC) $(OPTC) -c soluf.c  

mxm.o:		mxm.c
		$(peigs_CC) $(OPTC) -c mxm.c    

mxm2.o:		mxm2.c
		$(peigs_CC) $(OPTC) -c mxm2.c    

mxm5x.o:	mxm5x.c
		$(peigs_CC) $(OPTC) -c mxm5x.c    

mxm3.o:		mxm3.c
		$(peigs_CC) $(OPTC) -c mxm3.c

mxm4.o:		mxm4.c
		$(peigs_CC) $(OPTC) -c mxm4.c   

mxm5.o:		mxm5.c
		$(peigs_CC) $(OPTC) -c mxm5.c   

mxm8.o:		mxm8.c
		$(peigs_CC) $(OPTC) -c mxm8.c   

reducelst.o:	reducelst.c
		$(peigs_CC) $(OPTC) -c reducelst.c

sort.o:		sort.c
		$(peigs_CC) $(OPTC) -c sort.c

clustrf.o:	clustrf.c
		$(peigs_CC)  -c clustrf.c

inv_it2.o:	inv_it2.c
		$(peigs_CC) $(OPTC) -c inv_it2.c

inv_it5.o:	inv_it5.c
		$(peigs_CC) $(OPTC) -c inv_it5.c

mgs2.o:		mgs2.c
		$(peigs_CC) $(OPTC) -c mgs2.c

mgs1b.o:	mgs1b.c
		$(peigs_CC) $(OPTC) -c mgs1b.c

qsort.o:	qsort.c
		$(peigs_CC) $(OPTC) -c qsort.c

qsort1.o:	qsort1.c
		$(peigs_CC) $(OPTC) -c qsort1.c

memreq.o:	memreq.c
		$(peigs_CC) $(OPTC) -c memreq.c

memreq_f.o:	memreq_f.c
		$(peigs_CC) $(OPTC) -c memreq_f.c

dsteinsch.o:	dsteinsch.c
		$(peigs_CC) $(OPTC) -c  dsteinsch.c

mapdif.o:	mapdif.c
		$(peigs_CC) $(OPTC) -c mapdif.c    

mapdif1.o:	mapdif1.c
		$(peigs_CC) $(OPTC) -c mapdif1.c    

peigs_cmod.o:	peigs_cmod.f
		$(peigs_FC) -c peigs_cmod.f

clustrxx5.o:	clustrxx5.c
		$(peigs_CC) -c clustrxx5.c

clustrxx4.o:	clustrxx4.c
		$(peigs_CC) -c clustrxx4.c

peigs_tldlfact.o:	peigs_tldlfact.c
		$(peigs_CC) -c peigs_tldlfact.c

     CNFDIR := $(NWCHEM_TOP)/src/config
.PHONY:	 64_to_32 32_to_64
64_to_32:
	$(CNFDIR)/64_to_32 "*c *f *F"

32_to_64:
	$(CNFDIR)/32_to_64 "*c *f *F"

.c.o:;	$(peigs_CC) $(OPTC) -c $<
.f.o:;	$(peigs_FC) $(OPT) -c $<
.F.o:;  $(peigs_FC) $(OPT) -c $<
