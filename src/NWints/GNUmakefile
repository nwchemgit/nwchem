
# $Id: GNUmakefile,v 1.17 1994-08-18 04:52:22 d3e129 Exp $

HEADERS = api/int_tim.fh api/apiP.fh
LIBRARY = libnwints.a

include ../config/makefile.h

API_OBJ = ./api/int_init.o ./api/exact_mem.o ./api/int_mem.o \
          ./api/int_terminate.o ./api/int_chk_sh.o ./api/int_2e4c.o \
          ./api/int_l2e4c.o ./api/int_canon.o ./api/int_chk_init.o \
          ./api/int_1eall.o ./api/int_l1eall.o ./api/int_1eh1.o \
          ./api/int_lgen1e.o ./api/int_1eov.o ./api/cando_sp.o \
          ./api/int_1eke.o ./api/int_1epe.o ./api/int_pgen1e.o \
          ./api/int_2e3c.o ./api/int_2e2c.o ./api/int_l2e3c.o \
          ./api/int_l2e2c.o ./api/int_1e3ov.o ./api/int_l1e3ov.o \
          ./api/int_hf1sp.o ./api/int_nint.o \
          ./api/intd_1eov.o ./api/intd_1eh1.o ./api/intd_2e4c.o \
          ./api/intd_init.o ./api/exactd_mem.o ./api/int_acc.o

INT_OBJ = ./int/defNxyz.o ./int/getNxyz.o ./int/hf1.o ./int/hf1mke3.o \
          ./int/hf1mkr.o ./int/hf1set3.o ./int/hf2.o ./int/hf2mkr.o \
          ./int/hf2oi.o ./int/hf3OIs.o ./int/hfabc.o ./int/hfctr3.o \
          ./int/hferi.o ./int/hfkei.o ./int/hfmke.o ./int/hfmkr.o \
          ./int/hfnai.o ./int/hfset.o ./int/igamma.o ./int/hf_print.o \
          ./int/hfmkr_1.o ./int/hfmkr_2.o ./int/hfmkr_3.o \
          ./int/hfmkr_4.o ./int/hfmkr_5.o ./int/hfmkr_6.o \
          ./int/hfmkr_7.o ./int/hfmkr_8.o ./int/hfmkr_9.o

INTSP_OBJ = ./ints_sp/genr70.o ./ints_sp/init70.o ./ints_sp/qout70.o \
            ./ints_sp/sinfo.o ./ints_sp/sp0000.o ./ints_sp/sp0001.o \
            ./ints_sp/sp0011.o ./ints_sp/sp0101.o ./ints_sp/sp0111.o \
            ./ints_sp/sp1111.o ./ints_sp/tabgen.o 

USES_BLAS = ./api/int_hf1sp.F   ./api/int_l1e3ov.F   ./api/int_l1eall.F \
            ./api/int_l2e2c.F   ./api/int_l2e3c.F    ./api/int_l2e4c.F \
            ./api/int_lgen1e.F  ./api/int_pgen1e.F   ./api/intd_1eh1.F \
            ./api/intd_2e4c.F   ./deriv/hf1d.F       ./deriv/hf2d.F \
            ./ints_sp/sp1111.F

DERIV_OBJ  = ./deriv/hf1PRtoAB.o ./deriv/hf1d.o   ./deriv/hf2PRtoAB.o \
             ./deriv/hf2QStoCD.o ./deriv/hf2d.o   ./deriv/hfd2oi.o \
             ./deriv/hfderi.o    ./deriv/hfdkei.o ./deriv/hfdnai.o \
             ./deriv/hfefi.o

AUX_OBJ = ./aux/int_normalize.o ./aux/nmcoef.o ./aux/int_chk_sum.o

TEXAS_OBJ = 

OBJ = $(API_OBJ)  $(INT_OBJ) $(INTSP_OBJ)	$(DERIV_OBJ) $(AUX_OBJ) $(TEXAS_OBJ)

INT_SUBDIRS = api int aux ints_sp deriv texas

nwints_stamp:	include_stamp
	/bin/rm -f $(LIBRARY) 
	$(MAKE) intlibs
	$(AR) $(ARFLAGS) $(LIBRARY)  $(OBJ)
	$(RANLIB) $(LIBRARY) 
	cp -p $(LIBRARY) $(LIBDIR)
	touch nwints_stamp

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

ifdef USES_BLAS
sngl_to_dbl:
	$(CNFDIR)/sngl_to_dbl $(USES_BLAS)
dbl_to_sngl:
	$(CNFDIR)/dbl_to_sngl $(USES_BLAS)
else
sngl_to_dbl dbl_to_sngl:
	echo $@ : no conversion necessary
endif

redo rebuild:
	/bin/rm -f nwints_stamp
	$(MAKE) nwints_stamp
clean:
	/bin/rm -f $(LIBRARY) $(OBJ) core include_stamp nwints_stamp */*_stamp $(LIB_TARGETS)

realclean:	clean
	/bin/rm -f *~ \#*\#

debug:
	for intdir in $(INT_SUBDIRS) ; do \
		echo Making $$intdir objects ; \
		(cd $$intdir; $(MAKE) 'DEFINES = -DSUN -DINT_TIMING') || exit 1 ; \
	done
debugopt:
	for intdir in $(INT_SUBDIRS) ; do \
		echo Making $$intdir objects ; \
		(cd $$intdir; 	$(MAKE) 'FOPT=-O3 -Nl99' 'DEFINES = -DSUN -DINT_TIMING' ) || exit 1 ; \
	done

opt:
	for intdir in $(INT_SUBDIRS) ; do \
		echo Making $$intdir objects ; \
		(cd $$intdir; $(MAKE) 'FOPT=-O3 -Nl99' 'DEFINES = -DSUN') || exit 1 ; \
	done

