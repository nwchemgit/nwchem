#
# $Id: GNUmakefile,v 1.16 2004-10-14 21:58:46 edo Exp $
#
include ../../config/makefile.h
include ../DEFS
NOCASENS = CYGNUS MACX CYGWIN INTERIX
ifneq ($(peigs_TARGET),$(findstring $(peigs_TARGET), $(NOCASENS)))
	REMOVESF=TRUE
endif

MXOBJ	=	mxbrod.o \
		mxcombv1.o \
		mxsubs.o \
		sum.o \
		maxx.o \
		peigs_dgop.o

all:		$(MXOBJ)

.EXPORT_ALL_VARIABLES:

.PHONY :	clean

#
# Note that under Cygwin the rm commands are case-insensitive, so the
# commands below intended to remove intermediate .f files actually remove
# the original .F files!  Since intermediate files are not generated for
# Cygwin, this is avoided by commenting out the rm statements.
#
clean:		
		-rm -f *.o
ifdef REMOVESF
		-rm -f mxbrod.f
		-rm -f mxcombv1.f
		-rm -f mxsubs.f
		-rm -f sum.f peigs_dgop.f
endif

realclean:
		clean

.SUFFIXES:	.o .s .F .f .c

.c.o:;	$(peigs_CC) $(peigs_COPT) -c $<
.f.o:;	$(peigs_FC) $(peigs_FOPT) -c $<
.F.o:;	$(peigs_FC) $(peigs_FOPT) -c $<


