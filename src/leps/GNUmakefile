#$Id$
   OBJ_OPTIMIZE = leps_energy.o leps_gradient.o mc_driver.o mymc_input.o mc_data.o
   OBJ = leps_energy.o leps_gradient.o md_driver.o md_data.o mymd_input.o mc_driver.o mymc_input.o mc_data.o gasdev.o crossing.o mc_main.o md_set.o mc_init.o urand.o 
   LIBRARY = libleps.a

#
# ACCESS='append' is not standard and fails on IBM XL compilers.
# POSITION='append' is standard syntax but fails with FC=g77
#
# The default is POSITION='append' but if your compiler doesn't like this,
# then add your compiler name behind g77 (with a space).
#
ifeq ($(FC),$(findstring $(FC),g77 pgf77))
      LIB_DEFINES += -DUSE_OPEN_ACCESS
endif


include ../config/makefile.h
include ../config/makelib.h


