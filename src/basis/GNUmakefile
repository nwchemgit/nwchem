#$Id$

include ../config/makefile.h

       LIBRARY = libnwcutil.a
       HEADERS = bas.fh basP.fh bas_staticP.fh basdeclsP.fh geobasmapP.fh  \
                 bas_exndcf.fh bas_exndcf_dec.fh bas_exndcf_sfn.fh \
                 bas_ibs.fh bas_ibs_dec.fh bas_ibs_sfn.fh ../NWints/rel/rel_nwc.fh
           OBJ = basis.o bas_input.o bas_contrib.o bas_checksum.o basisP.o bas_blas.o bas_blasP.o bas_vec_info.o

   LIB_TARGETS = testbasis testbasis.o libcheck

     USES_BLAS = bas_blas.F bas_blasP.F basis.F

include ../config/makelib.h

compare:	compare.o
	$(LINK.f) -o $@ compare.o  $(LIBS) 
libcheck:	libcheck.o
	$(LINK.f) -o $@ libcheck.o $(LIBS) 
testbasis:	testbasis.o
	$(LINK.f) -o $@ testbasis.o  $(LIBS) ../stubs.o -lnwctask
libraries_bse:	getlibr.py 
	mkdir -p libraries.bse
	cd libraries.bse ; ../getlibr.py

basP.fh:	bas_staticP.fh
	@touch basP.fh

