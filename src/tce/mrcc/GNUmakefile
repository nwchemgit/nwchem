
SUBDIRS = include ccsd ccsd_sub ccsd_t uss

OBJ = tce_mrcc_input.o

#OBJ_OPTIMIZE = tce_mrcc_energy.o tce_mrcc_create_cas.o tce_mrcc_readref.o tce_mrcc_calcexcit.o tce_mrcc_sortint.o \
#OBJ_OPTIMIZE = tce_mrcc_energy.o tce_mrcc_init.o tce_mrcc_create_cas.o tce_mrcc_readref.o tce_mrcc_calcexcit.o tce_mrcc_sortint.o 
#OBJ_OPTIMIZE =  tce_mrcc_energy.o tce_mrcc_create_cas.o tce_mrcc_readref.o tce_mrcc_calcexcit.o tce_mrcc_sortint.o \

OBJ_OPTIMIZE =  tce_mrcc_energy.o tce_mrcc_readref.o tce_mrcc_calcexcit.o tce_mrcc_sortint.o \
		tce_mrcc_tile.o tce_mrcc_macpy.o tce_mrcc_fillindexarray.o tce_mrcc_uhf.o create_mr_evl_sorted.o \
		pgcreatefile.o tce_mrcc_mo2e.o tce_mrcc_ao2e.o tce_subnxtask.o tce_mrcc_sublock.o eaf_mrcc.o \
		tce_t3ac_offset.o tce_t4ac_offset.o tce_mrcc_refs_sym.o tce_mrcc_force_orben.o \
		tce_mrcc_diis_init.o tce_mrcc_heff.o  \
		tce_mrcc_diis.o tce_mrcc_jacobi_t1.o tce_mrcc_jacobi_t2.o \
		tce_mrcc_gacp.o tce_mrcc_daxfile.o tce_mrcc_diis_tidy.o tce_jacobi_sub_t1.o tce_jacobi_sub_t2.o \
		tce_mrcc_c1.o tce_mrcc_c1_offset.o tce_mrcc_c2.o tce_mrcc_c2_offset.o \
		tce_mrcc_mul_1.o tce_mrcc_mul_2.o \
		tce_mrcc_diis_save.o tce_mrcc_diis_new.o  tce_mrcc_create_cas1.o atpart_nalength_mrcc.o \
                new_ga4ind_N5uhf_mrcc_lm2.o 
		
#		tce_mrcc_print.o \

LIB_INCLUDES = -I../include -I./include

LIBRARY = libtce.a

#USES_BLAS = tce_mrcc_heff.F tce_mrcc_diis.F tce_mrcc_diis_new.F tce_mrcc_diis_init.F tce_mrcc_diis_save.F tce_mrcc_diis_tidy.F \
#	    tce_mrcc_daxfile.F tce_mrcc_uhf.F \
#	    tce_mrcc_mul_1.F tce_mrcc_mul_2.F tce_mrcc_c1.F tce_mrcc_c1_offset.F tce_mrcc_c2.F tce_mrcc_c2_offset.F \
#	    tce_mrcc_mo2e.F tce_mrcc_ao2e.F  tce_mrcc_calcexcit.F tce_mrcc_sortint.F \
#	    tce_mrcc_tile.F tce_mrcc_macpy.F tce_mrcc_fillindexarray.F

#USES_BLAS = tce_mrcc_energy.F tce_mrcc_create_cas.F tce_mrcc_readref.F tce_mrcc_calcexcit.F tce_mrcc_sortint.F \

USES_BLAS = tce_mrcc_energy.F tce_mrcc_readref.F tce_mrcc_calcexcit.F tce_mrcc_sortint.F \
	    tce_mrcc_tile.F tce_mrcc_macpy.F tce_mrcc_fillindexarray.F tce_mrcc_uhf.F create_mr_evl_sorted.F \
	    pgcreatefile.F tce_mrcc_mo2e.F tce_mrcc_ao2e.F tce_subnxtask.F tce_mrcc_sublock.F eaf_mrcc.F \
	    tce_t3ac_offset.F tce_t4ac_offset.F tce_mrcc_refs_sym.F tce_mrcc_force_orben.F \
	    tce_mrcc_diis_init.F tce_mrcc_heff.F  \
	    tce_mrcc_diis.F tce_mrcc_jacobi_t1.F tce_mrcc_jacobi_t2.F \
	    tce_mrcc_gacp.F tce_mrcc_daxfile.F tce_mrcc_diis_tidy.F tce_jacobi_sub_t1.F tce_jacobi_sub_t2.F \
	    tce_mrcc_c1.F tce_mrcc_c1_offset.F tce_mrcc_c2.F tce_mrcc_c2_offset.F \
	    tce_mrcc_mul_1.F tce_mrcc_mul_2.F \
	    tce_mrcc_diis_save.F tce_mrcc_diis_new.F tce_mrcc_create_cas1.F \
            new_ga4ind_N5uhf_mrcc_lm2.F 
		
#	    tce_mrcc_print.F \


LIB_DEFINES = 

#              -DDEBUG_PRINT

HEADERS = 

ifdef MRCC_METHODS
#ckbn mrcc-r-1 -3
      #LIB_DEFINES += -DMRCC_METHODS -DMRCC_LOCAL_FOCK -DMRCC_LOCAL_T1
      LIB_DEFINES += -DMRCC_METHODS #-DMRCC_LOCAL_FOCK -DMRCC_LOCAL_T1
#LIB_DEFINES += -DMRCC_LOCAL_FOCK -DMRCC_LOCAL_T1
endif

include ../../config/makefile.h
include ../../config/makelib.h

