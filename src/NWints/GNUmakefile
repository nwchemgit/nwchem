
# $Id: GNUmakefile,v 1.1 1994-04-01 22:50:47 d3e129 Exp $

LIBRARY = libnwint.a

include ../config/makefile.h

API_OBJ    = api/int_init.o

INT_OBJ = int/defNxyz.o int/getNxyz.o int/hf1.o int/hf1mke3.o \
          int/hf1mkr.o int/hf1set3.o int/hf2.o int/hf2mkr.o \
          int/hf2oi.o int/hf3OIs.o int/hfabc.o int/hfctr3.o \
          int/hferi.o int/hfkei.o int/hfmke.o int/hfmkr.o \
          int/hfnai.o int/hfset.o int/igamma.o 

DERIV_OBJ  = 

AUX_OBJ = 

OBJ = $(API_OBJ)  $(INT_OBJ) $(DERIV_OBJ) $(AUX_OBJ)

INT_SUBDIRS = api int

$(LIBRARY):	intlibs
	/bin/rm -f $@
	$(AR) $(ARFLAGS) $@ $(OBJ)
	$(RANLIB) $@
	cp -p $(LIBRARY) $(LIBDIR)

ifdef HEADERS
include_stamp:	$(HEADERS)
	cp -p $(HEADERS) $(INCDIR)
	touch include_stamp
else
include_stamp:
	touch include_stamp
endif

intlibs:
	for intdir in $(INT_SUBDIRS) ; do \
		echo Making $$intdir objects ; \
		(cd $$intdir; $(MAKE)) || exit 1 ; \
	done

clean:
	/bin/rm -f $(LIBRARY) $(OBJ) core include_stamp $(LIB_TARGETS)


realclean:	clean
	/bin/rm -f *~ \#*\#
