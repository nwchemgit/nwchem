# $Id: GNUmakefile,v 1.17 1996-09-24 21:41:39 d3e129 Exp $

# optimization breaks on the SGITFP ... there's no point
# optimizing this stuff anyway.

# well actually in large basis sets the atomic guess can take
# a long time ... optimize except on the SGI for now.

#ifeq ($(TARGET),SGITFP)
  XXX = OBJ
#else
  XXX = OBJ_OPTIMIZE
#endif

       $(XXX) = orderd.o \
                atcond.o \
                denat.o \
		hamild.o \
		outpud.o \
		teigd.o \
		atomd.o \
		denmad.o \
		jacod.o \
		pdfded.o \
		tracd.o \
		cmergd.o \
		densid.o \
		densid_z.o \
		oeigd.o \
		shalfd.o \
		trafsd.o \
		creded.o \
		square.o \
		tramad.o \
		datoms.o \
		guess_mem.o  \
		starcd.o \
		writel.o \
		printbasis.o \
                atscf_atscf.o


  USES_BLAS = atomd.F \
              atscf_atscf.F


  LIBRARY = libguess.a


include ../config/makefile.h
include ../config/makelib.h


# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(teigd.o):	./cguess.fh 
$(LIBRARY_PATH)(printbasis.o):	./cguess.fh 
$(LIBRARY_PATH)(outpud.o):	./cguess.fh 
$(LIBRARY_PATH)(osatod.o):	$(INCDIR)/geom.fh ./cguess.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(oeigd.o):	./cguess.fh 
$(LIBRARY_PATH)(hamild.o):	./cguess.fh 
$(LIBRARY_PATH)(guess_mem.o):	./cguess.fh 
$(LIBRARY_PATH)(densid_z.o):	$(INCDIR)/inp.fh ./atomz.fh 
$(LIBRARY_PATH)(denat.o):	$(INCDIR)/util.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cguess.fh $(INCDIR)/bas.fh ./atomz.fh 
$(LIBRARY_PATH)(datoms.o):	$(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cguess.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(creded.o):	$(INCDIR)/global.fh ./cguess.fh 
$(LIBRARY_PATH)(atscf_atscf.o):	$(INCDIR)/util.fh $(INCDIR)/stdio.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cguess.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(atomd.o):	$(INCDIR)/inp.fh ./cguess.fh 
