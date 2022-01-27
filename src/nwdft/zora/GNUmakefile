# $Id$
#
	OBJ_OPTIMIZE = \
	dft_zora_utils.o\
	dft_zora_utils_so.o\
	dft_zora_rhos.o\
	dft_zora_rhos_so.o\
	gridNuclearPotential.o\
	gridHartreePotential.o\
	gridQpqPotential.o\
        grid_moments.o\
	zora_getv_EFGZ4_SR.o\
	calc_zora_EFGZ4_SR.o\
	zora_getv_EFGZ4_SO.o\
	calc_zora_EFGZ4_SO.o\
	dft_zora_NMR.o\
	zora_getv_NMRCS_SR.o\
	calc_zora_NMRCS_SR.o\
	dft_zora_EPR.o\
	dft_zora_EPR-NMR_tools.o\
	zora_getv_EPR.o\
	calc_zora_EPR.o\
        dft_zora_Hyperfine.o\
        zora_getv_HFine.o\
        calc_zora_HFine.o\
	zora_getv_sf.o\
	zora_getv_so.o\
	calc_EFG.o\
	calc_zora_sf.o\
	calc_zora_so.o\
	calc_zora_facs.o\
	get_ints_zora.o zora_recalc.o

     LIBRARY = libnwdft.a

 LIB_DEFINES =

LIB_INCLUDES = -I../include -I../grid

   USES_BLAS = dft_zora_NMR.F dft_zora_EPR-NMR_tools.F zora_tabcd.F dft_zora_EPR.F dft_zora_Hyperfine.F zora_tabcd_EFG.F dft_zora_utils_so.F zora_getv_HFine.F


HEADERS = 

include ../../config/makefile.h
include ../../config/makelib.h
