
# $Id: GNUmakefile,v 1.2 1994-04-04 20:33:43 d3e129 Exp $

LIBRARY = libnwints.a

include ../config/makefile.h

API_OBJ = api/int_init.o api/exact_mem.o api/int_mem.o api/int_terminate.o api/int_chk_sh.o \
          api/int_2e4c.o api/int_l2e4c.o api/int_canon.o api/int_chk_init.o api/int_1eall.o \
          api/int_l1eall.o

INT_OBJ = int/defNxyz.o int/getNxyz.o int/hf1.o int/hf1mke3.o \
          int/hf1mkr.o int/hf1set3.o int/hf2.o int/hf2mkr.o \
          int/hf2oi.o int/hf3OIs.o int/hfabc.o int/hfctr3.o \
          int/hferi.o int/hfkei.o int/hfmke.o int/hfmkr.o \
          int/hfnai.o int/hfset.o int/igamma.o 

DERIV_OBJ  = 

AUX_OBJ = aux/int_normalize.o aux/nmcoef.o

OBJ = $(API_OBJ)  $(INT_OBJ) $(DERIV_OBJ) $(AUX_OBJ)

INT_SUBDIRS = api int aux

$(LIBRARY):	$(OBJ)
	/bin/rm -f $@
	$(MAKE) intlibs
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
