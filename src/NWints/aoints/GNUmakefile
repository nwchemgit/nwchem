
# $Id


include ../../config/makefile.h
     OBJ_OPTIMIZE = dft_aoints.o       wrt_dft_aoints.o int_1e_sifs.o \
                    int_so_sifs.o      int_2e_sifs.o    int_2e_sifs_a.o \
                    sifs_2e_task.o     int_2e_sifs_b.o  nadct_trace.o \
                    rdhcid.o           egrad_trace.o    rd1mat.o sif2ga.o \
                    sif2da_2e_sort.o   sif2arr.o        rd_d2bl.o \
                    d2geri_trace.o     dint_block_trc.o sequ.o \
                    sopgrdtrc.o	       asif2ga.o \
                    print_dint_block.o print_soblock.o  trc_soblock.o \
                    int_mom_sifs.o     hdoverlap.o 

     LIBRARY = libnwints.a

LIB_INCLUDES = -I../../nwdft/grid -I../../nwdft/include \
 -I../../ddscf # I../../nwdft/util

include ../../config/makelib.h

# FOPTIONS += -Minform=inform -Mdclchk 

