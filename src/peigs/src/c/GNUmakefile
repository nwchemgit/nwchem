include ../../../config/makefile.h
include ../../DEFS

COBJ  =  b_ortho.o \
         chol_bcst.o \
         choleski9.o \
         ci_entry.o \
         clustrf3.o \
         clustrxx3.o \
         conjug22.o \
         de_sym.o \
         exit.o \
         exit2.o \
         forLL.o \
         forLU.o \
         gmax.o \
         inv_it3.o \
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
         pstebz9.o \
         pstein3.o \
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
         reducemap.o


all:		$(COBJ)

.EXPORT_ALL_VARIABLES:

.PHONY :	clean

clean:		
		-rm -f *.o
		-rm -f *.c.*

realclean:	
		clean

soort.o:	soort.c
		$(peigs_CC) $(peigs_COPT) -c soort.c

lower_mxm.o:	lower_mxm.c
		$(peigs_CC) -c lower_mxm.c

shellsort.o:	shellsort.c
		$(peigs_CC) $(peigs_COPT) -c shellsort.c

treesort.o:	treesort.c
		$(peigs_CC) -c treesort.c

gmax.o:		gmax.c
		$(peigs_CC) $(peigs_COPT) -c gmax.c

forLL.o:	forLL.c
		$(peigs_CC) $(peigs_COPT) -c forLL.c

forLU.o:	forLU.c
		$(peigs_CC) $(peigs_COPT) -c forLU.c

ci_entry.o:	ci_entry.c
		$(peigs_CC) $(peigs_COPT) -c ci_entry.c

util.o:		util.c
		$(peigs_CC) $(peigs_COPT) -c util.c

pdspgv_c.o:	pdspgv_c.c
		$(peigs_CC) $(peigs_COPT) -c pdspgv_c.c

pdspev_c.o:	pdspev_c.c
		$(peigs_CC) $(peigs_COPT) -c pdspev_c.c

pdspevx.o:	pdspevx.c
		$(peigs_CC) $(peigs_COPT) -c pdspevx.c

pdspgvx.o:	pdspgvx.c
		$(peigs_CC) $(peigs_COPT) -c pdspgvx.c

exit.o:		exit.c
		$(peigs_CC) $(peigs_COPT) -c exit.c

exit2.o:	exit2.c
		$(peigs_CC) $(peigs_COPT) -c exit2.c

pxerbla.o:	pxerbla.c
		$(peigs_CC) $(peigs_COPT) -c pxerbla.c

chol_bcst.o:	chol_bcst.c
		$(peigs_CC) $(peigs_COPT) -c chol_bcst.c

upperxfull.o:	upperxfull.c
		$(peigs_CC) $(peigs_COPT) -c upperxfull.c

pipe_fut.o:	pipe_fut.c
		$(peigs_CC) $(peigs_COPT) -c pipe_fut.c

choleski9.o:	choleski9.c
		$(peigs_CC) $(peigs_COPT) -c choleski9.c

tred22.o:	tred22.c
		$(peigs_CC) $(peigs_COPT) -c tred22.c

pstebz9.o:	pstebz9.c
		$(peigs_CC) $(peigs_COPT) -c pstebz9.c

inverse.o:	inverse.c
		$(peigs_CC) $(peigs_COPT)  -c inverse.c

conjug.o:	conjug.c
		$(peigs_CC) $(peigs_COPT) -c conjug.c 

pmmLUL.o:	pmmLUL.c
		$(peigs_CC) $(peigs_COPT) -c pmmLUL.c 

pmmLSL.o:	pmmLSL.c
		$(peigs_CC) $(peigs_COPT) -c pmmLSL.c 

soluf.o:	soluf.c
		$(peigs_CC) $(peigs_COPT) -c soluf.c  

mxm.o:		mxm.c
		$(peigs_CC) $(peigs_COPT) -c mxm.c    

mxm2.o:		mxm2.c
		$(peigs_CC) $(peigs_COPT) -c mxm2.c    

mxm5x.o:	mxm5x.c
		$(peigs_CC) $(peigs_COPT) -c mxm5x.c    

mxm3.o:		mxm3.c
		$(peigs_CC) $(peigs_COPT) -c mxm3.c

mxm4.o:		mxm4.c
		$(peigs_CC) $(peigs_COPT) -c mxm4.c   

mxm5.o:		mxm5.c
		$(peigs_CC) $(peigs_COPT) -c mxm5.c   

mxm8.o:		mxm8.c
		$(peigs_CC) $(peigs_COPT) -c mxm8.c   

reducelst.o:	reducelst.c
		$(peigs_CC) $(peigs_COPT) -c reducelst.c

sort.o:		sort.c
		$(peigs_CC) $(peigs_COPT) -c sort.c

clustrf.o:	clustrf.c
		$(peigs_CC)  -c clustrf.c

clustrxx.o:	clustrxx.c
		$(peigs_CC) $(peigs_COPT) -c clustrxx.c

inv_it2.o:	inv_it2.c
		$(peigs_CC) $(peigs_COPT) -c inv_it2.c

mgs2.o:		mgs2.c
		$(peigs_CC) $(peigs_COPT) -c mgs2.c

mgs1b.o:	mgs1b.c
		$(peigs_CC) $(peigs_COPT) -c mgs1b.c

qsort.o:	qsort.c
		$(peigs_CC) $(peigs_COPT) -c qsort.c

qsort1.o:	qsort1.c
		$(peigs_CC) $(peigs_COPT) -c qsort1.c

memreq.o:	memreq.c
		$(peigs_CC) $(peigs_COPT) -c memreq.c

memreq_f.o:	memreq_f.c
		$(peigs_CC) $(peigs_COPT) -c memreq_f.c

pstein.o:	pstein.c
		$(peigs_CC) $(peigs_COPT) -c pstein.c

dsteinsch.o:	dsteinsch.c
		$(peigs_CC) $(peigs_COPT) -c  dsteinsch.c

mapdif.o:	mapdif.c
		$(peigs_CC) $(peigs_COPT) -c mapdif.c    

mapdif1.o:	mapdif1.c
		$(peigs_CC) $(peigs_COPT) -c mapdif1.c    

peigs_cmod.o:	peigs_cmod.f
		$(peigs_FC) -c peigs_cmod.f


.c.o:;	$(peigs_CC) $(peigs_COPT) -c $<
.f.o:;	$(peigs_FC) $(peigs_FOPT) -c $<

